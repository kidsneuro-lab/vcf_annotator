from __future__ import annotations

import csv
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
        records = {}
        for rec in reader:
            alts = rec.alts or ["."]
            for alt in alts:
                key = f"{rec.chrom}:{rec.pos}:{alt}"
                records.setdefault(key, []).append(rec)
        return records


def read_tsv(path: Path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def info_value(record, key):
    value = record.info.get(key)
    if isinstance(value, tuple):
        if not value:
            return None
        if len(value) == 1:
            return value[0]
    return value
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
    entries = records["chr8:18210109:C"]
    assert len(entries) == 2

    info_by_transcript = {
        info_value(rec, "SJ_TRANSCRIPT"): rec for rec in entries
    }
    assert set(info_by_transcript) == {"NM_001160170.4", "XM_047422397.1"}
    assert "MANE" in entries[0].header.info

    expected_don = {"NM_001160170.4": "-72", "XM_047422397.1": "-72"}
    expected_acc = {"NM_001160170.4": "NA", "XM_047422397.1": "18"}

    for transcript, rec in info_by_transcript.items():
        assert info_value(rec, "SJ_VARIANT_TYPE") == "snp"
        assert info_value(rec, "SJ_DDON") == expected_don[transcript]
        assert info_value(rec, "SJ_DACC") == expected_acc[transcript]

    assert info_value(info_by_transcript["NM_001160170.4"], "SJ_DDON") == "-72"
    assert info_value(info_by_transcript["XM_047422397.1"], "SJ_DACC") == "18"
    assert info_value(info_by_transcript["NM_001160170.4"], "MANE") == 1
    assert info_value(info_by_transcript["XM_047422397.1"], "MANE") == 0

    rows = read_tsv(tsv_path)
    c_rows = [row for row in rows if row["CHROM"] == "chr8" and row["POS"] == "18210109" and row["ALT"] == "C"]
    assert len(c_rows) == 2
    assert {row["SJ_TRANSCRIPT"] for row in c_rows} == {"NM_001160170.4", "XM_047422397.1"}
    assert all(row["SJ_DDON"] != "NA" for row in c_rows)
    mane_values = {row["SJ_TRANSCRIPT"]: row["MANE"] for row in c_rows}
    assert mane_values == {"NM_001160170.4": "1", "XM_047422397.1": "0"}


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
    minus = records["chr19:58347029:T"]
    assert len(minus) == 1
    transcripts = info_value(minus[0], "SJ_TRANSCRIPT")
    assert transcripts == "NM_130786.4"
    assert info_value(minus[0], "SJ_VARIANT_TYPE") == "snp"
    assert "MANE" not in minus[0].header.info


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
    first = records["chr8:18210109:C"][0]

    assert info_value(first, "CLN_CLNSIG") == "Pathogenic"
    assert info_value(first, "CLN_REVIEW") == "criteria_provided"

    benign = records["chr19:58347029:T"][0]
    benign_alt = records["chr19:58347029:G"][0]
    assert info_value(benign, "CLN_CLNSIG") == "Benign"
    assert info_value(benign, "CLN_REVIEW") == "no_assertion"

    assert info_value(benign_alt, "CLN_CLNSIG") in {"NA", "Benign"}


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

    with pysam.VariantFile(str(output_vcf)) as reader:
        records = list(reader)

    assert len(records) == 6
    for rec in records:
        assert len(rec.alts or []) <= 1
