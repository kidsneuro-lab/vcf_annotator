# VCF Annotator

CLI tool for annotating VCF files with splice-junction distances and INFO fields from one or more auxiliary VCFs. The tool writes a single-ALT VCF and can also emit a compact TSV summary.

## Installation

This project uses Python 3.10+ and `uv` for dependency management.

```bash
uv sync --extra test
```

Run commands through the managed environment:

```bash
uv run python -m vcf_annotator.cli --help
```

For reference download, genePred conversion, MANE setup, and Docker usage, see
[docs/USAGE.md](docs/USAGE.md).

## CLI Interface

```text
usage: python -m vcf_annotator.cli --input INPUT --output OUTPUT
                                      [--annotate-dist GENE_PRED;PREFIX[;MANE]]
                                      [--annotate-vcf VCF;PREFIX;FIELD[,FIELD...]]
                                      [--tsv TSV] [--normalise] [--verbose]
```

Required arguments:

- `--input`: input VCF path. Plain or bgzipped VCF input is supported by `pysam`.
- `--output`: annotated output VCF path.

Optional arguments:

- `--annotate-dist <gene_pred>;<prefix>[;<mane_transcripts>]`: annotate each ALT allele with splice-junction distances using a genePred transcript file. The prefix becomes the INFO tag name, uppercased.
- `--annotate-vcf <vcf>;<prefix>;<fields>`: copy selected INFO values from a matching auxiliary VCF. Repeat this option to use multiple sources. Output fields are named `<PREFIX>_<FIELD>`.
- `--tsv <path>`: also write tab-delimited output with `CHROM`, `POS`, `ID`, `REF`, `ALT`, and annotation columns.
- `--normalise`: split multiallelic input with `bcftools norm -m -any` through `pysam`.
- `--verbose`: enable debug logging.

Splice annotation values use this format:

```text
Allele|Transcript|Gene|VariantType|DonorDist|DonorRegionType|DonorRegionNo|AcceptorDist|AcceptorRegionType|AcceptorRegionNo[|MANE]
```

## Examples

Annotate splice-junction distances with a MANE transcript list:

```bash
uv run python -m vcf_annotator.cli \
  --input tests/data/input.vcf \
  --output /tmp/annotated.splice.vcf \
  --annotate-dist "tests/data/sample.genePred;SJ;tests/data/mane.tsv" \
  --tsv /tmp/annotated.splice.tsv
```

Annotate from a custom VCF:

```bash
uv run python -m vcf_annotator.cli \
  --input tests/data/input.vcf \
  --output /tmp/annotated.custom.vcf \
  --annotate-vcf "tests/data/custom.vcf;CLN;CLNSIG,REVIEW"
```

Combine splice and custom VCF annotations:

```bash
uv run python -m vcf_annotator.cli \
  --input tests/data/input.vcf \
  --output /tmp/annotated.all.vcf \
  --annotate-dist "tests/data/sample.genePred;SJ" \
  --annotate-vcf "tests/data/custom.vcf;CLN;CLNSIG,REVIEW" \
  --tsv /tmp/annotated.all.tsv \
  --verbose
```

## Development

Run the test suite:

```bash
uv run pytest
```

Run only the CLI integration tests:

```bash
uv run pytest tests/test_cli.py
```

Build the development container image:

```bash
docker build -t vcf-annotator .
```

## Test Data

The repository includes small sampled fixtures under `tests/data/`. The notes below describe how selected reference fixtures were generated.

### Sample GenePred

Source:

- <https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2025_08/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz>

Convert annotation GTF to genePred format:

```bash
gtfToGenePred -ignoreGroupsWithoutExons -genePredExt -includeVersion \
  GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
  GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred
```

Obtain a subset of genePred rows for `EMD`, `HBB`, `COL2A1`, and `DYSF`:

```bash
awk -F"\t" '$12~/^EMD$|^HBB$|^COL2A1$|^DYSF$/' \
  GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred \
  > GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred

rm GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred
```

### Sample MANE List

Source:

- <https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz>

```bash
gunzip -c MANE.GRCh38.v1.5.summary.txt.gz | \
  awk -F"\t" 'NR==1 || $4=="COL2A1" || $4=="DYSF" || $4=="EMD" || $3=="HBB"' \
  > MANE.GRCh38.v1.5.summary.sample.txt

rm MANE.GRCh38.v1.5.summary.txt.gz
```
