# Polygenic Score Pipeline

Computes a polygenic score from per-chromosome PLINK2 filesets stored in a
Google Cloud Storage bucket. For each chromosome (1-22 and X) the script:

1. Downloads the chromosome's fileset from the bucket.
2. Extracts the SNPs of interest with PLINK2.
3. Deletes the downloaded files to save disk space.

It then merges the extracted chromosomes into a single fileset and scores it
using a weights file.

## Requirements

- Python 3.6 
- [`plink2`](https://www.cog-genomics.org/plink/2.0/) on your `PATH`
- `gsutil` (part of the [Google Cloud SDK](https://cloud.google.com/sdk))
  on your `PATH`, authenticated with read access to the bucket
- Enough free disk space for one full chromosome fileset at a time, plus the
  extracted outputs

Check that the tools are available:

```bash
plink2 --version
gsutil version
```

## Input files

Place these in the directory you run the script from.

| File | Description |
|------|-------------|
| `snps_sle.txt` | SNP list, one variant ID per line, in `CHR:POS:REF:ALT` style (e.g. `1:12345:A:G`). Used both to decide which chromosomes to process and as PLINK2's `--extract` list. |
| `extraction_weight.txt` | Score weights file. Columns 1, 2 and 3 must be variant ID, effect allele and effect weight, in that order. |

The bucket must contain files named `chr<N>_trimmed_merged.pgen`, `.pvar` and
`.psam` for each chromosome, under the directory you configure below.

## Configuration

Open `polygenic_score_pipeline.py` and fill in the constants near the top:

```python
BUCKET_ID = "my-bucket-name"       # bucket name, without gs://
DIRECTORY = "path/to/genotypes"    # folder inside the bucket
```

Files are fetched from `gs://<BUCKET_ID>/<DIRECTORY>/chr<N>_trimmed_merged*`.

`VARIANT_ID` is defined but not currently used by the script. Other
constants (`SNP_LIST_FILE`, `SCORE_WEIGHTS_FILE`, `MERGED_PREFIX`,
`SCORE_PREFIX`) can be changed if you want different file names.

## Running

From the directory containing your input files:

```bash
python polygenic_score_pipeline.py
```

## Output

| File | Description |
|------|-------------|
| `extracted_data_<chrom>_chrom.pgen/.pvar/.psam` | Extracted SNPs for each chromosome |
| `merge_list.txt` | Prefixes of the filesets merged into the first one |
| `merged_genome.pgen/.pvar/.psam` | All extracted chromosomes merged |
| `polygenic_score.sscore` | Final scores, one row per sample |
| `*.log` | PLINK2 logs for each step |

## Resuming and re-running

The script is safe to re-run. If any `extracted_data_<chrom>_chrom` file
already exists for a chromosome, that chromosome is reused instead of being
downloaded and extracted again. To force a chromosome to be redone, delete
its `extracted_data_<chrom>_chrom.*` files first.

The merge and scoring steps always run again, overwriting the previous
merged and score outputs.

If a run fails partway through a chromosome, delete any partial
`extracted_data_<chrom>_chrom.*` files for it before re-running, since a
partial file would be treated as already extracted.

## Troubleshooting

- **`NotImplementedError` / chromosomes skipped unexpectedly:** the script
  decides which chromosomes to process by reading the first `:`-separated
  field of each ID in `snps_sle.txt`. If your IDs are rsIDs or use another
  format, `has_chr_i` needs to be adapted.
- **`gsutil` access errors:** run `gcloud auth list` to confirm you are
  authenticated and that your account can read the bucket.
- **`ValueError: No extracted chromosomes available to merge`:** no
  chromosome had matching SNPs, or none were extracted. Check
  `snps_sle.txt` and the log output.
- **PLINK2 fails on merge with only one chromosome:** the merge list will be
  empty. Use the single `extracted_data_<chrom>_chrom` fileset directly for
  scoring.
- **Out of disk space:** downloaded files are removed after each chromosome,
  but a single chromosome can be large. Free space or run on a machine with
  more storage.
