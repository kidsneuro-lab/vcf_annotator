from __future__ import annotations

import argparse
import logging
import sys
import tempfile
from pathlib import Path
from typing import List, Optional

import pysam

from .annotators import CustomVcfAnnotator, SpliceJunctionDistanceAnnotator
from .chromosome import ChromosomeMapper
from .processor import ProcessorConfig, VariantProcessor
from .transcripts import build_transcript_index

LOG = logging.getLogger("vcf_annotator")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Annotate VCF files with splice junction distances and custom VCF fields."
    )
    parser.add_argument("--input", required=True, help="Input VCF file (optionally bgzipped).")
    parser.add_argument("--output", required=True, help="Output VCF file.")
    parser.add_argument("--annotate-dist", help="Splice annotation config: <gene_pred>;<prefix>[;<mane_transcripts>].")
    parser.add_argument(
        "--annotate-vcf",
        action="append",
        dest="annotate_vcf",
        default=[],
        help="Custom VCF annotation config: <vcf.gz>;<prefix>;<info fields comma separated>. Can be supplied multiple times.",
    )
    parser.add_argument("--tsv", help="Optional TSV output path.")
    parser.add_argument("--normalise", action="store_true", help="Normalise input VCF using bcftools norm -m -any.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    tsv_path = Path(args.tsv).expanduser().resolve() if args.tsv else None

    LOG.info("Input VCF: %s", input_path)
    LOG.info("Output VCF: %s", output_path)
    if tsv_path:
        LOG.info("TSV output: %s", tsv_path)
    if args.normalise:
        LOG.info("Normalisation enabled via bcftools norm.")

    effective_input = input_path
    temp_files: List[Path] = []

    try:
        if args.normalise:
            effective_input = _normalise_vcf(input_path)
            temp_files.append(effective_input)
            LOG.info("Normalised VCF written to temporary file %s", effective_input)

        annotators = []

        with pysam.VariantFile(str(effective_input)) as reader:
            chrom_mapper = ChromosomeMapper(reader.header.contigs.keys())

        if args.annotate_dist:
            dist_config = _parse_dist_config(args.annotate_dist)
            transcript_index = build_transcript_index(
                dist_config["gene_pred"],
                chrom_mapper,
                dist_config.get("mane"),
            )
            annotators.append(
                SpliceJunctionDistanceAnnotator(
                    transcript_index,
                    dist_config["prefix"],
                    include_mane=dist_config.get("mane") is not None,
                )
            )
            LOG.info(
                "Loaded transcript model (%d transcripts) using prefix %s",
                sum(len(lst) for lst in transcript_index.by_chrom.values()),
                dist_config["prefix"],
            )

        for entry in args.annotate_vcf:
            custom_config = _parse_custom_config(entry)
            annotators.append(
                CustomVcfAnnotator(
                    custom_config["vcf_path"],
                    custom_config["prefix"],
                    custom_config["fields"],
                )
            )
            LOG.info(
                "Configured custom VCF annotator for %s with prefix %s (fields: %s)",
                custom_config["vcf_path"],
                custom_config["prefix"],
                ",".join(custom_config["fields"]),
            )

        if not annotators:
            LOG.warning(
                "No annotators configured; producing single-allelic output without additional annotations."
            )

        processor = VariantProcessor(
            ProcessorConfig(
                input_path=effective_input,
                output_path=output_path,
                annotators=annotators,
                tsv_path=tsv_path,
            )
        )
        processor.run()

        LOG.info("Annotated VCF written to %s", output_path)
        if tsv_path:
            LOG.info("TSV annotations written to %s", tsv_path)
        return 0
    finally:
        for path in temp_files:
            if path.exists():
                path.unlink(missing_ok=True)


def _parse_dist_config(value: str):
    parts = value.split(";")
    if len(parts) < 2:
        raise ValueError("Invalid --annotate-dist value; expected at least <gene_pred>;<prefix>.")
    gene_pred = Path(parts[0]).expanduser().resolve()
    prefix = parts[1].strip()
    mane = Path(parts[2]).expanduser().resolve() if len(parts) >= 3 and parts[2].strip() else None

    if not gene_pred.exists():
        raise FileNotFoundError(f"Gene pred file not found: {gene_pred}")
    if mane and not mane.exists():
        raise FileNotFoundError(f"MANE transcript file not found: {mane}")

    return {"gene_pred": gene_pred, "prefix": prefix, "mane": mane}


def _parse_custom_config(value: str):
    parts = value.split(";", 2)
    if len(parts) != 3:
        raise ValueError("Invalid --annotate-vcf value; expected <vcf.gz>;<prefix>;<info fields>.")
    vcf_path = Path(parts[0]).expanduser().resolve()
    prefix = parts[1].strip()
    fields = [item.strip() for item in parts[2].split(",") if item.strip()]

    if not vcf_path.exists():
        raise FileNotFoundError(f"Custom VCF not found: {vcf_path}")
    if not fields:
        raise ValueError("No INFO fields provided for custom VCF annotator.")

    return {"vcf_path": vcf_path, "prefix": prefix, "fields": fields}


def _normalise_vcf(path: Path) -> Path:
    if not hasattr(pysam, "bcftools"):
        raise RuntimeError("pysam bcftools wrapper is unavailable; cannot normalise.")

    tmp_fd, tmp_path = tempfile.mkstemp(prefix="vcf_annotator_norm_", suffix=".vcf")
    Path(tmp_path).unlink(missing_ok=True)  # bcftools will create the file

    try:
        pysam.bcftools.norm(
            "--multiallelics",
            "-any",
            "-o",
            tmp_path,
            str(path),
        )
    except Exception as exc:  # pragma: no cover - depends on external tool
        raise RuntimeError(f"bcftools norm failed: {exc}") from exc

    return Path(tmp_path)



if __name__ == "__main__":
    sys.exit(main())
