from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

from vcf_annotator.annotators.base import VariantContext
from vcf_annotator.annotators.splice_distance import SpliceJunctionDistanceAnnotator, _format_distance
from vcf_annotator.chromosome import ChromosomeMapper
from vcf_annotator.transcripts import build_transcript_index

DATA_DIR = Path(__file__).parent / "data"


class DummyRecord:
    """Lightweight stand-in for pysam VariantRecord objects."""

    def __init__(self, chrom: str, pos: int, ref: str, alts: tuple[str, ...]):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alts = alts
        self.id = "test"
        self.info = {}


@pytest.fixture(scope="session")
def transcript_index():
    gene_pred = DATA_DIR / "sample.genePred"
    mane = DATA_DIR / "mane.tsv"
    mapper = ChromosomeMapper(["chr8", "chr19"])
    return build_transcript_index(gene_pred, mapper, mane_path=mane)


def _expected_row(
    transcript: str,
    gene: str,
    variant_type: str,
    ddon: str,
    ddon_region_type: str,
    ddon_region_no: str,
    dacc: str,
    dacc_region_type: str,
    dacc_region_no: str,
    mane_flag: int,
) -> dict[str, object]:
    return {
        "SJ_TRANSCRIPT": transcript,
        "SJ_GENE": gene,
        "SJ_VARIANT_TYPE": variant_type,
        "SJ_DDON": ddon,
        "SJ_DDON_REGION_TYPE": ddon_region_type,
        "SJ_DDON_REGION_NO": ddon_region_no,
        "SJ_DACC": dacc,
        "SJ_DACC_REGION_TYPE": dacc_region_type,
        "SJ_DACC_REGION_NO": dacc_region_no,
        "MANE": mane_flag,
    }


@pytest.mark.parametrize(
    ("chrom", "pos", "ref", "alt", "expected_rows"),
    [
        (
            "chr8",
            18210181,
            "G",
            "C",
            [
                _expected_row("NM_001160170.4", "NAT1", "snp", "1", "intron", "1", "-1011", "intron", "1", 1),
                _expected_row("XM_047422397.1", "NAT1", "snp", "1", "intron", "5", "-1011", "intron", "5", 0),
            ],
        )
    ],
)
def test_splice_distance_annotations(transcript_index, chrom, pos, ref, alt, expected_rows):
    annotator = SpliceJunctionDistanceAnnotator(transcript_index, "sj", include_mane=True)
    record = DummyRecord(chrom, pos, ref, (alt,))
    context = VariantContext(record, alt, 0)

    result = annotator.annotate(context)

    expected_map = {row["SJ_TRANSCRIPT"]: row for row in expected_rows}
    observed_map = {row["SJ_TRANSCRIPT"]: dict(row) for row in result.rows}

    assert len(result.tsv_rows) == len(expected_rows)
    
    for tsv_row in result.tsv_rows:
        transcript = tsv_row["SJ_TRANSCRIPT"]
        for key, value in expected_map[transcript].items():
            assert dict(tsv_row)[key] == value, f"Mismatch for {transcript} on {key}"
    
    assert observed_map == expected_map

def test_format_distance_helper():
    assert _format_distance(5) == "5"
    assert _format_distance(None) == "NA"
