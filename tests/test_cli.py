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
    entries = records["chr12:48006054:A"]
    assert len(entries) == 1

    sj_value = info_value(entries[0], "SJ")
    sj_entries = list(sj_value) if isinstance(sj_value, tuple) else [sj_value]

    expected_entries = [
        "A|XM_017018828.1|COL2A1|snp|1|intron|1|-325|intron|1|0",
        "A|XM_017018829.1|COL2A1|snp|1|intron|1|-325|intron|1|0",
        "A|XM_017018830.1|COL2A1|snp|1|intron|1|-325|intron|1|0",
    ]

    assert sorted(sj_entries) == sorted(expected_entries)
    assert "SJ" in entries[0].header.info
    assert "MANE" not in entries[0].header.info

    rows = read_tsv(tsv_path)
    c_rows = [row for row in rows if row["CHROM"] == "chr12" and row["POS"] == "48006054" and row["ALT"] == "A"]
    assert len(c_rows) == 1
    sj_entries_tsv = c_rows[0]["SJ"].split(",")
    assert sorted(sj_entries_tsv) == sorted(expected_entries)


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
    minus = records["chr12:48006054:A"]
    assert len(minus) == 1
    sj_value = info_value(minus[0], "SJ")
    sj_entries = list(sj_value) if isinstance(sj_value, tuple) else [sj_value]
    assert len(sj_entries) == 3
    fields = sj_entries[0].split("|")
    assert len(fields) == 10
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
    first = records["chr12:48006054:A"][0]

    assert info_value(first, "CLN_CLNSIG") == "Pathogenic"
    assert info_value(first, "CLN_REVIEW") == "criteria_provided"

    benign = records["chr12:47987610:A"][0]
    assert info_value(benign, "CLN_CLNSIG") == "Benign"
    assert info_value(benign, "CLN_REVIEW") == "no_assertion"


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
