# Usage

This guide shows how to prepare splice-distance reference inputs and run a
versioned Docker image against an input VCF.

## Prepare a Reference Directory

Create a directory for the downloaded reference files and generated outputs:

```bash
mkdir -p references
```

## Download the GTF and Create a genePred File

Download the GRCh38 RefSeq GTF:

```bash
curl -L \
  -o references/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2025_08/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
```

Convert the GTF to genePred format with the UCSC `gtfToGenePred` tool. The
`gtfToGenePred` executable must be installed and available on your `PATH`.

```bash
gtfToGenePred -ignoreGroupsWithoutExons -genePredExt -includeVersion \
  references/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
  references/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred
```

The generated `.genepred` file is passed to `--annotate-dist` when running the
annotator.

## Download the MANE Annotation

Download the MANE GRCh38 summary file:

```bash
curl -L \
  -o references/MANE.GRCh38.v1.5.summary.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz
```

The annotator can read the gzipped MANE summary file directly, so decompression
is not required.

## Pull a Versioned Docker Image

Set the image version you want to run, then pull that tag:

```bash
VERSION=0.1.0
IMAGE=ghcr.io/kidsneuro-lab/vcf_annotator:${VERSION}

docker pull "${IMAGE}"
```

Replace `0.1.0` with the released version you want to use.

## Run the Docker Image

The container entrypoint is the `vcf-annotator` CLI, so pass CLI arguments after
the image name. This example mounts the current working directory at `/data` and
writes output files next to the input VCF:

```bash
docker run --rm \
  -v "$PWD":/data \
  "${IMAGE}" \
  --input /data/input.vcf \
  --output /data/annotated.vcf \
  --annotate-dist "/data/references/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred;SJ;/data/references/MANE.GRCh38.v1.5.summary.txt.gz" \
  --tsv /data/annotated.tsv
```

The `--annotate-dist` value has three semicolon-separated fields:

- genePred path inside the container.
- INFO field prefix, such as `SJ`.
- optional MANE summary path inside the container.

Use paths under the mounted `/data` directory so the container can read the input
VCF and reference files and write the annotated outputs.

## Interpret SpliceDistance INFO Fields

The splice-distance annotation is written to the INFO tag named by the
`--annotate-dist` prefix. For example, the prefix `SJ` creates an `SJ` INFO
field. The field contains one comma-separated entry per overlapping transcript
for each ALT allele. Each entry is pipe-delimited:

```text
Allele|Transcript|Gene|VariantType|DonorDist|DonorRegionType|DonorRegionNo|AcceptorDist|AcceptorRegionType|AcceptorRegionNo[|MANE]
```

Field meanings:

- `Allele`: ALT allele being annotated.
- `Transcript`: overlapping transcript ID from the genePred file, or `NA` when no transcript overlaps.
- `Gene`: gene symbol/name from the genePred file, or `NA` when no transcript overlaps.
- `VariantType`: variant class reported by the annotator, such as `snp`, `insertion`, `deletion`, or `complex`.
- `DonorDist`: signed distance to the nearest splice donor boundary for that transcript.
- `DonorRegionType`: `exon`, `intron`, or `NA`, describing where the variant anchor used for `DonorDist` falls.
- `DonorRegionNo`: 1-based exon or intron number for `DonorRegionType` in transcript order.
- `AcceptorDist`: signed distance to the nearest splice acceptor boundary for that transcript.
- `AcceptorRegionType`: `exon`, `intron`, or `NA`, describing where the variant anchor used for `AcceptorDist` falls.
- `AcceptorRegionNo`: 1-based exon or intron number for `AcceptorRegionType` in transcript order.
- `MANE`: included only when a MANE file is provided; `1` means the transcript is present in the MANE summary and `0` means it is not.

Distance sign interpretation:

- Donor distances are positive when the variant anchor is intronic after the donor boundary and negative when it is exonic before that donor boundary.
- Acceptor distances are positive when the variant anchor is exonic after the acceptor boundary and negative when it is intronic before that acceptor boundary.
- A distance of `1` or `-1` means the anchor is immediately adjacent to the splice boundary on the indicated side.
- `NA` means the relevant boundary cannot be computed for that transcript, such as a donor in the last exon, an acceptor in the first exon, or no overlapping transcript.

For variants spanning more than one base, the annotator evaluates both the start
and end anchors of the variant interval and reports the closest donor and
acceptor distances. If two anchors are equally close, an intronic assignment is
preferred over an exonic assignment.
