from __future__ import annotations

from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from vcf_annotator.annotators.base import AnnotationResult, VariantContext
from vcf_annotator.annotators.custom_vcf import CustomVcfAnnotator
from vcf_annotator.annotators.splice_distance import SpliceJunctionDistanceAnnotator, _format_distance
from vcf_annotator.chromosome import ChromosomeMapper
from vcf_annotator.transcripts import Region, Transcript
from vcf_annotator.variant_utils import VariantCoordinates, compute_variant_bounds


class DummyRecord:
    """Lightweight stand-in for pysam records in unit tests."""

    def __init__(self, chrom: str, pos: int, ref: str, alts: tuple[str, ...]):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.id = "test"
        self.info = {}


class DummyTranscriptIndex:
    def __init__(self, transcripts):
        self.transcripts = transcripts

    def fetch(self, chrom: str, start0: int, end0: int):
        matched = []
        for transcript in self.transcripts:
            if transcript.chrom != chrom:
                continue
            if transcript.overlaps(start0, end0):
                matched.append(transcript)
        return matched


def build_transcript(name: str = "TX1", strand: str = "+"):
    exons = [Region(0, 100, 1, "exon"), Region(200, 300, 2, "exon")]
    introns = [Region(100, 200, 1, "intron")]
    chrom = "chr1"
    if strand == "-":
        exons = [Region(200, 300, 1, "exon"), Region(0, 100, 2, "exon")]
        introns = [Region(100, 200, 1, "intron")]
    return Transcript(
        name=name,
        chrom=chrom,
        strand=strand,
        gene="GENE1",
        tx_start=0,
        tx_end=400,
        exons=exons,
        introns=introns,
        mane=name.endswith("MANE"),
    )


def test_annotation_result_merge_handles_scalars_and_sequences():
    first = AnnotationResult(info={"A": "X"}, tsv_rows=[{"A": "X"}])
    second = AnnotationResult(info={"A": "Y", "B": ["Z", "Q"]}, tsv_rows=[{"B": "Z"}])

    first.merge(second)

    assert first.info["A"] == ("X", "Y")
    assert first.info["B"] == ["Z", "Q"]
    assert first.tsv_rows == [{"A": "X"}, {"B": "Z"}]


def test_variant_context_variant_type_classification():
    ctx = VariantContext(DummyRecord("chr1", 101, "A", ("G",)), "G", 0)
    assert ctx.variant_type() == "snp"

    ctx = VariantContext(DummyRecord("chr1", 101, "A", ("AG",)), "AG", 0)
    assert ctx.variant_type() == "ins"

    ctx = VariantContext(DummyRecord("chr1", 101, "AG", ("A",)), "A", 0)
    assert ctx.variant_type() == "del"

    ctx = VariantContext(DummyRecord("chr1", 101, "AG", ("CT",)), "CT", 0)
    assert ctx.variant_type() == "delins"


def test_compute_variant_bounds_handles_variant_types():
    assert compute_variant_bounds(100, "A", "G") == VariantCoordinates(99, 100)
    assert compute_variant_bounds(100, "A", "GA") == VariantCoordinates(99, 100)
    assert compute_variant_bounds(100, "GA", "G") == VariantCoordinates(99, 101)


def test_chromosome_mapper_translates_between_styles():
    mapper_chr = ChromosomeMapper(["chr1", "chr2"])
    assert mapper_chr.to_vcf("1") == "chr1"
    assert mapper_chr.to_vcf("MT") == "chrM"

    mapper_plain = ChromosomeMapper(["1", "2"])
    assert mapper_plain.to_vcf("chr1") == "1"
    assert mapper_plain.to_vcf("chrM") == "MT"


def test_splice_annotator_reports_distance_and_region():
    transcript = build_transcript()
    index = DummyTranscriptIndex([transcript])
    annotator = SpliceJunctionDistanceAnnotator(index, "sj")

    record = DummyRecord("chr1", 201, "A", ("C",))
    context = VariantContext(record, "C", 0)
    result = annotator.annotate(context)

    assert result.info["SJ_TRANSCRIPT"] == transcript.name
    assert result.info["SJ_GENE"] == transcript.gene
    assert result.info["SJ_VARIANT_TYPE"] == "snp"
    assert result.info["SJ_DACC"] == "0"
    assert result.info["SJ_DDON"] == "100"
    assert result.info["SJ_DDON_REGION_TYPE"] == "exon"
    assert result.info["SJ_DDON_REGION_NO"] == "2"
    assert result.tsv_rows[0]["dacc"] == "0"


def test_splice_annotator_handles_missing_transcripts():
    annotator = SpliceJunctionDistanceAnnotator(DummyTranscriptIndex([]), "sj")
    record = DummyRecord("chr1", 101, "A", ("C",))
    context = VariantContext(record, "C", 0)
    result = annotator.annotate(context)

    assert result.info["SJ_TRANSCRIPT"] == "NA"
    assert result.info["SJ_DDON"] == "NA"


def test_format_distance_helper():
    assert _format_distance(5) == "5"
    assert _format_distance(None) == "NA"


def write_temp_vcf(path: Path, contig: str = "chr1") -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", "VCFv4.2")
    header.contigs.add(contig)
    header.info.add("CLNSIG", 1, "String", " Clinical significance")

    with pysam.VariantFile(str(path), "w", header=header) as handle:
        rec = handle.new_record(contig=contig, start=99, stop=100, alleles=("A", "G"))
        rec.info["CLNSIG"] = "Benign"
        handle.write(rec)


def test_custom_vcf_annotator_fetches_matching_records(tmp_path):
    external = tmp_path / "external.vcf"
    write_temp_vcf(external)

    annotator = CustomVcfAnnotator(external, "cln", ["CLNSIG"])
    context = VariantContext(DummyRecord("1", 100, "A", ("G",)), "G", 0)

    result = annotator.annotate(context)

    assert result.info["CLN_CLNSIG"] == "Benign"
    assert result.tsv_rows == [{"CLN_CLNSIG": "Benign"}]

    annotator.close()


def test_custom_vcf_annotator_returns_na_for_missing(tmp_path):
    external = tmp_path / "external.vcf"
    write_temp_vcf(external, contig="1")

    annotator = CustomVcfAnnotator(external, "cln", ["CLNSIG"])
    context = VariantContext(DummyRecord("chr1", 100, "A", ("T",)), "T", 0)

    result = annotator.annotate(context)

    assert result.info["CLN_CLNSIG"] == "NA"
    assert result.tsv_rows == []

    annotator.close()
