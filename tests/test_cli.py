from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from vcf_annotator import cli

DATA_DIR = Path(__file__).parent / "data"


def write_ensembl_style_vcf(source: Path, target: Path) -> None:
    text = source.read_text()

    def replace_chrom(match):
        chrom = match.group(1)
        return "MT" if chrom == "M" else chrom

    target.write_text(re.sub(r"chr([0-9]+|X|Y|M)\b", replace_chrom, text))


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


def test_splice_annotation_with_ensembl_input_and_ucsc_genepred_preserves_input_chrom(tmp_path):
    input_vcf = tmp_path / "input.ensembl.vcf"
    write_ensembl_style_vcf(DATA_DIR / "input.vcf", input_vcf)
    output_vcf = tmp_path / "annotated.ensembl.vcf"
    tsv_path = tmp_path / "annotated.ensembl.tsv"

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
    assert "12:48006054:A" in records
    assert "chr12:48006054:A" not in records

    entry = records["12:48006054:A"][0]
    sj_value = info_value(entry, "SJ")
    sj_entries = list(sj_value) if isinstance(sj_value, tuple) else [sj_value]
    assert "A|XM_017018828.1|COL2A1|snp|1|intron|1|-325|intron|1|0" in sj_entries

    with pysam.VariantFile(str(output_vcf)) as reader:
        assert "12" in reader.header.contigs
        assert "chr12" not in reader.header.contigs

    rows = read_tsv(tsv_path)
    assert any(row["CHROM"] == "12" and row["POS"] == "48006054" for row in rows)
    assert not any(row["CHROM"] == "chr12" for row in rows)


def test_splice_annotation_with_gzipped_mane(tmp_path):
    input_vcf = DATA_DIR / "input.vcf"
    mane_gzip = tmp_path / "mane-compressed.tsv"
    output_vcf = tmp_path / "annotated_gzip_mane.vcf"

    with (
        (DATA_DIR / "mane.tsv").open("rb") as source,
        gzip.open(mane_gzip, "wb") as target,
    ):
        target.write(source.read())

    run_cli(
        [
            "--input",
            str(input_vcf),
            "--output",
            str(output_vcf),
            "--annotate-dist",
            f"{DATA_DIR / 'sample.genePred'};SJ;{mane_gzip}",
        ]
    )

    records = read_record_map(output_vcf)
    entries = records["chr12:48006054:A"]
    sj_value = info_value(entries[0], "SJ")
    sj_entries = list(sj_value) if isinstance(sj_value, tuple) else [sj_value]

    expected_entries = [
        "A|XM_017018828.1|COL2A1|snp|1|intron|1|-325|intron|1|0",
        "A|XM_017018829.1|COL2A1|snp|1|intron|1|-325|intron|1|0",
        "A|XM_017018830.1|COL2A1|snp|1|intron|1|-325|intron|1|0",
    ]

    assert sorted(sj_entries) == sorted(expected_entries)


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

    assert len(records) == 11
    for rec in records:
        assert len(rec.alts or []) <= 1


def test_input_without_contig_header_adds_record_contigs_to_output(tmp_path):
    input_vcf = tmp_path / "no_contigs.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\tvar1\tA\tG\t.\t.\t.\n"
    )
    output_vcf = tmp_path / "annotated.vcf"

    run_cli(
        [
            "--input",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    with pysam.VariantFile(str(output_vcf)) as reader:
        assert "1" in reader.header.contigs
        records = list(reader)

    assert len(records) == 1
    assert records[0].chrom == "1"
    assert records[0].pos == 100
    assert records[0].alts == ("G",)
