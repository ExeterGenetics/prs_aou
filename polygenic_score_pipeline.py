"""Compute a polygenic score from per-chromosome PLINK2 files.

Pipeline steps, per chromosome (1-22 and X):

1. Download the trimmed/merged PLINK2 fileset from a Google Cloud bucket.
2. Extract the SNPs of interest with PLINK2.
3. Delete the downloaded files to save disk space.

Afterwards, the per-chromosome extracts are merged into a single fileset
and a polygenic score is computed from a weights file.

Requirements:
    * ``gsutil`` and ``plink2`` available on the PATH.
    * ``snps_sle.txt`` (SNP IDs to extract) and ``extraction_weight.txt``
      (score weights) in the working directory.
"""

import os
import subprocess

# Fill these in before running.
BUCKET_ID = ""
DIRECTORY = ""
VARIANT_ID = ""

CHROMOSOMES = list(range(1, 23)) + ["X"]
EXTENSIONS = ["pgen", "pvar", "psam"]
SNP_LIST_FILE = "snps_sle.txt"
SCORE_WEIGHTS_FILE = "extraction_weight.txt"
MERGE_LIST_FILE = "merge_list.txt"
MERGED_PREFIX = "merged_genome"
SCORE_PREFIX = "polygenic_score"

# Files that may be left behind by a download and PLINK2 run.
CLEANUP_SUFFIXES = [
    ".pgen", ".pvar", ".psam", ".log", ".pgen_.gstmp", ".pvar_.gstmp",
]


def has_chr_i(chrom, snp_file=SNP_LIST_FILE):
    """Return True if the SNP list contains variants on ``chrom``.

    Assumes variant IDs in the first column of ``snp_file`` start with the
    chromosome, in a ``CHR:POS:REF:ALT`` style (an optional ``chr`` prefix
    is ignored), e.g. ``1:12345:A:G`` or ``chrX:67890:C:T``.

    Args:
        chrom: Chromosome number (int) or ``"X"``.
        snp_file: Path to the SNP list, one variant ID per line.

    Returns:
        bool: Whether the chromosome has SNPs of interest.
    """
    target = str(chrom).lower()

    with open(snp_file, "r") as file:
        for line in file:
            fields = line.split()
            if not fields:
                continue

            line_chrom = fields[0].split(":")[0].lower()
            if line_chrom.startswith("chr"):
                line_chrom = line_chrom[3:]

            if line_chrom == target:
                return True

    return False


def is_extracted(out_prefix):
    """Check whether any extracted PLINK2 file exists for a prefix.

    Args:
        out_prefix: Prefix of the extracted fileset.

    Returns:
        bool: True if any of the pgen/pvar/psam files exist.
    """
    return any(os.path.exists(f"{out_prefix}.{ext}") for ext in EXTENSIONS)


def download_chromosome(prefix):
    """Copy a chromosome's fileset from the bucket if not present locally.

    Args:
        prefix: Prefix of the remote fileset, e.g. ``chr1_trimmed_merged``.

    Raises:
        subprocess.CalledProcessError: If ``gsutil`` fails.
    """
    if os.path.exists(f"{prefix}.pgen"):
        return
    subprocess.run(
        [
            "gsutil", "-m", "cp",
            f"gs://{BUCKET_ID}/{DIRECTORY}/{prefix}*",
            ".",
        ],
        check=True,
    )


def extract_snps(prefix, out_prefix):
    """Extract the SNPs listed in ``SNP_LIST_FILE`` using PLINK2.

    Args:
        prefix: Prefix of the input fileset.
        out_prefix: Prefix for the extracted output fileset.

    Raises:
        subprocess.CalledProcessError: If PLINK2 fails.
    """
    subprocess.run(
        [
            "plink2",
            "--pgen", f"{prefix}.pgen",
            "--pvar", f"{prefix}.pvar",
            "--psam", f"{prefix}.psam",
            "--extract", SNP_LIST_FILE,
            "--make-pgen",
            "--out", out_prefix,
        ],
        check=True,
    )


def remove_downloaded_files(prefix):
    """Delete downloaded chromosome files to conserve disk space.

    Args:
        prefix: Prefix of the downloaded fileset.
    """
    for suffix in CLEANUP_SUFFIXES:
        file_path = f"{prefix}{suffix}"
        if os.path.exists(file_path):
            os.remove(file_path)


def extract_all_chromosomes():
    """Download and extract SNPs for every chromosome with relevant SNPs.

    Chromosomes already extracted are reused rather than reprocessed, and
    chromosomes with no SNPs in the SNP list are skipped.

    Returns:
        list[str]: Prefixes of all extracted per-chromosome filesets.
    """
    extracted_prefixes = []

    for chrom in CHROMOSOMES:
        out_prefix = f"extracted_data_{chrom}_chrom"

        if is_extracted(out_prefix):
            print(f"chr{chrom} already extracted, adding to merge list")
            extracted_prefixes.append(out_prefix)
            continue

        if not has_chr_i(chrom):
            print(f"Skipping over chr{chrom}: not in SNP list")
            continue

        prefix = f"chr{chrom}_trimmed_merged"

        download_chromosome(prefix)
        extract_snps(prefix, out_prefix)
        extracted_prefixes.append(out_prefix)
        remove_downloaded_files(prefix)

    return extracted_prefixes


def merge_chromosomes(extracted_prefixes):
    """Merge per-chromosome filesets into one PLINK2 fileset.

    The first prefix is used as the base fileset; the rest are listed in
    ``MERGE_LIST_FILE`` and merged into it.

    Args:
        extracted_prefixes: Prefixes of the filesets to merge.

    Raises:
        ValueError: If there is nothing to merge.
        subprocess.CalledProcessError: If PLINK2 fails.
    """
    if not extracted_prefixes:
        raise ValueError("No extracted chromosomes available to merge.")

    with open(MERGE_LIST_FILE, "w") as merge_file:
        for prefix in extracted_prefixes[1:]:
            merge_file.write(prefix + "\n")

    subprocess.run(
        [
            "plink2",
            "--pfile", extracted_prefixes[0],
            "--pmerge-list", MERGE_LIST_FILE,
            "--make-pgen",
            "--out", MERGED_PREFIX,
        ],
        check=True,
    )


def compute_polygenic_score():
    """Compute the polygenic score on the merged fileset with PLINK2.

    Uses columns 1, 2 and 3 of ``SCORE_WEIGHTS_FILE`` as the variant ID,
    effect allele and effect weight, respectively.

    Raises:
        subprocess.CalledProcessError: If PLINK2 fails.
    """
    subprocess.run(
        [
            "plink2",
            "--pgen", f"{MERGED_PREFIX}.pgen",
            "--pvar", f"{MERGED_PREFIX}.pvar",
            "--psam", f"{MERGED_PREFIX}.psam",
            "--score", SCORE_WEIGHTS_FILE, "1", "2", "3",
            "--out", SCORE_PREFIX,
        ],
        check=True,
    )


def main():
    """Run the full extract, merge and scoring pipeline."""
    extracted_prefixes = extract_all_chromosomes()
    merge_chromosomes(extracted_prefixes)
    compute_polygenic_score()


if __name__ == "__main__":
    main()
