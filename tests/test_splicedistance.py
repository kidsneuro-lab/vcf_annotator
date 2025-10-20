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


def _expected_entry(
    allele: str,
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
) -> str:
    return "|".join(
        [
            allele,
            transcript,
            gene,
            variant_type,
            ddon,
            ddon_region_type,
            ddon_region_no,
            dacc,
            dacc_region_type,
            dacc_region_no,
            str(mane_flag),
        ]
    )


@pytest.mark.parametrize(
    ("chrom", "pos", "ref", "alt", "expected_rows"),
    [
        (
            "chr8",
            18210181,
            "G",
            "C",
                [
                    _expected_entry("C", "NM_001160170.4", "NAT1", "snp", "1", "intron", "1", "-1012", "intron", "1", 1),
                    _expected_entry("C", "XM_047422397.1", "NAT1", "snp", "1", "intron", "4", "-1012", "intron", "4", 0),
                ],
        )
    ],
)
def test_splice_distance_annotations(transcript_index, chrom, pos, ref, alt, expected_rows):
    annotator = SpliceJunctionDistanceAnnotator(transcript_index, "sj", include_mane=True)
    record = DummyRecord(chrom, pos, ref, (alt,))
    context = VariantContext(record, alt, 0)

    result = annotator.annotate(context)

    assert len(result.rows) == 1
    assert len(result.tsv_rows) == 1

    observed_entries = result.rows[0]["SJ"].split(",")

    assert sorted(observed_entries) == sorted(expected_rows)
    assert result.tsv_rows[0]["SJ"] == result.rows[0]["SJ"]

def test_format_distance_helper():
    assert _format_distance(5) == "5"
    assert _format_distance(None) == "NA"
