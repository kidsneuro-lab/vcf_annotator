from __future__ import annotations

from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from vcf_annotator import cli

DATA_DIR = Path(__file__).parent / "data"


def run_cli(args):
    exit_code = cli.main(args)
    assert exit_code == 0


def read_record_map(vcf_path: Path):
    with pysam.VariantFile(str(vcf_path)) as reader:
        return {f"{rec.chrom}:{rec.pos}": rec for rec in reader}


def _info_values(record, key):
    value = record.info.get(key)
    if value is None:
        return []
    if isinstance(value, tuple):
        return list(value)
    return [value]


def test_splice_annotation_with_mane(tmp_path):
    input_vcf = DATA_DIR / "input.vcf"
    output_vcf = tmp_path / "annotated.vcf"
    tsv_path = tmp_path / "annotated.tsv"

    run_cli(
        [
            "--input",
            str(input_vcf),
            "--output",
            str(output_vcf),
            "--annotate-dist",
            f"{DATA_DIR / 'sample.genePred'};SJ;{DATA_DIR / 'mane.tsv'}",
            "--tsv",
            str(tsv_path),
        ]
    )

    records = read_record_map(output_vcf)
    first = records["chr8:18210109"]

    transcripts = _info_values(first, "SJ_TRANSCRIPT")[0]
    assert "NM_001160170.4" in transcripts
    assert _info_values(first, "SJ_VARIANT_TYPE")[0] == "snp"
    assert _info_values(first, "SJ_DDON")[0] != "NA"
    assert _info_values(first, "SJ_DACC")[0] != "NA"

    tsv_content = tsv_path.read_text().splitlines()
    assert tsv_content[0].startswith("CHROM")
    assert any("SJ_DDON" in line for line in tsv_content[1:])


def test_splice_annotation_without_mane(tmp_path):
    input_vcf = DATA_DIR / "input.vcf"
    output_vcf = tmp_path / "annotated_nomane.vcf"

    run_cli(
        [
            "--input",
            str(input_vcf),
            "--output",
            str(output_vcf),
            "--annotate-dist",
            f"{DATA_DIR / 'sample.genePred'};SJ",
        ]
    )

    records = read_record_map(output_vcf)
    minus = records["chr19:58347029"]
    transcripts = _info_values(minus, "SJ_TRANSCRIPT")[0]
    assert "NM_130786.4" in transcripts
    assert _info_values(minus, "SJ_VARIANT_TYPE")[0] == "snp"


def test_custom_vcf_annotation(tmp_path):
    input_vcf = DATA_DIR / "input.vcf"
    custom_vcf = DATA_DIR / "custom.vcf"
    output_vcf = tmp_path / "custom.vcf.out"

    run_cli(
        [
            "--input",
            str(input_vcf),
            "--output",
            str(output_vcf),
            "--annotate-vcf",
            f"{custom_vcf};CLN;CLNSIG,REVIEW",
        ]
    )

    records = read_record_map(output_vcf)
    first = records["chr8:18210109"]

    assert _info_values(first, "CLN_CLNSIG")[0] == "Pathogenic"
    assert _info_values(first, "CLN_REVIEW")[0] == "criteria_provided"

    multi = records["chr19:58347029"]
    clnsig_values = _info_values(multi, "CLN_CLNSIG")
    assert clnsig_values[0] == "Benign"
    assert clnsig_values[1] in {"NA", "Benign"}


def test_no_annotators_copies_input(tmp_path):
    input_vcf = DATA_DIR / "input.vcf"
    output_vcf = tmp_path / "copied.vcf"

    run_cli(
        [
            "--input",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert output_vcf.exists()
    assert output_vcf.read_text() == input_vcf.read_text()
