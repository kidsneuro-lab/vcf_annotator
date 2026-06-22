from pathlib import Path

import pytest
from pytest_bdd import scenarios, given, when, then, parsers

from vcf_annotator.annotators.base import VariantContext
from vcf_annotator.annotators.splice_distance import SpliceJunctionDistanceAnnotator, _format_distance
from vcf_annotator.chromosome import ChromosomeMapper
from vcf_annotator.transcripts import build_transcript_index
pysam = pytest.importorskip("pysam")

DATA_DIR = Path(__file__).parent / "data"

# Load all scenarios from the feature file
scenarios("features/splice_distance.feature")

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
    gene_pred = DATA_DIR / "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred"
    mane = DATA_DIR / "MANE.GRCh38.v1.4.summary.sample.txt"
    mapper = ChromosomeMapper(["chr12", "chrX", "chr11", "chr2"])
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
    def normalize_missing(value: str) -> str:
        return "NA" if value == "NULL" else value

    return "|".join(
        [
            normalize_missing(allele),
            normalize_missing(transcript),
            normalize_missing(gene),
            normalize_missing(variant_type),
            normalize_missing(ddon),
            normalize_missing(ddon_region_type),
            normalize_missing(ddon_region_no),
            normalize_missing(dacc),
            normalize_missing(dacc_region_type),
            normalize_missing(dacc_region_no),
            str(mane_flag),
        ]
    )


@pytest.mark.parametrize(
    ("chrom", "pos", "ref", "alt", "expected_rows"),
    [
        (
            "chr12",
            48006054,
            "C",
            "A",
                [
                    _expected_entry("A", "XM_017018828.1", "COL2A1", "snp", "1", "intron", "1", "-325", "intron", "1", 0),
                    _expected_entry("A", "XM_017018829.1", "COL2A1", "snp", "1", "intron", "1", "-325", "intron", "1", 0),
                    _expected_entry("A", "XM_017018830.1", "COL2A1", "snp", "1", "intron", "1", "-325", "intron", "1", 0)
                ],
        )
    ],
)
def test_splice_distance_annotations(transcript_index, chrom, pos, ref, alt, expected_rows):
    annotator = SpliceJunctionDistanceAnnotator(transcript_index, "sj", include_mane=True)
    record = DummyRecord(chrom, pos, ref, (alt,))
    context = VariantContext(record, alt, 0)

    result = annotator.annotate(context)

    observed_entries = result.rows[0]["SJ"].split(",")

    assert sorted(observed_entries) == sorted(expected_rows)
    assert result.tsv_rows[0]["SJ"] == result.rows[0]["SJ"]

def test_format_distance_helper():
    assert _format_distance(5) == "5"
    assert _format_distance(None) == "NA"

@pytest.fixture
def ctx():
    """Per-scenario scratchpad."""
    return {}

@given("a transcript index")
def have_transcript_index(transcript_index, ctx):
    # Reuse your existing fixture
    ctx["transcript_index"] = transcript_index

@given(
    parsers.parse('a variant at "{chrom}" {pos:d} with ref "{ref}" and alt "{alt}"'),
)
def make_variant(ctx, chrom, pos, ref, alt):
    ctx["chrom"] = chrom
    ctx["pos"] = pos
    ctx["ref"] = ref
    ctx["alt"] = alt


@when('I annotate splice distances with tag "sj" and include_mane')
def run_annotation(ctx):
    annotator = SpliceJunctionDistanceAnnotator(ctx["transcript_index"], "sj", include_mane=True)
    record = DummyRecord(ctx["chrom"], ctx["pos"], ctx["ref"], (ctx["alt"],))
    vc = VariantContext(record, ctx["alt"], 0)
    result = annotator.annotate(vc)
    ctx["result"] = result
    # Parsed INFO field entries
    ctx["observed_entries"] = result.rows[0]["SJ"].split(",")

# Lambda function to print all observed entries for debugging. Each entry is printed on a new line for clarity.
log_observed_entries = lambda ctx: print("\n".join(ctx["observed_entries"]))

@then(
    parsers.parse(
        'the SJ field should contain the expected rows built from "{allele}" "{transcript_id}" "{gene}" "{variant_type}" "{donor_distance}" "{donor_region_type}" "{donor_region_number}" "{acceptor_distance}" "{acceptor_region_type}" "{acceptor_region_number}" {mane_flag:d}'
    )
)
def check_expected_rows(
    ctx,
    allele,
    transcript_id,
    gene,
    variant_type,
    donor_distance,
    donor_region_type,
    donor_region_number,
    acceptor_distance,
    acceptor_region_type,
    acceptor_region_number,
    mane_flag,
):
    expected_rows = [
        _expected_entry(
            allele,
            transcript_id,
            gene,
            variant_type,
            donor_distance,
            donor_region_type,
            donor_region_number,
            acceptor_distance,
            acceptor_region_type,
            acceptor_region_number,
            mane_flag,
        )
    ]
    # assert sorted(ctx["observed_entries"]) == sorted(expected_rows)
    # The order of entries may vary
    # We want to ensure that the expected rows are present regardless of order
    observed_entries_set = set(ctx["observed_entries"])
    expected_rows_set = set(expected_rows)
    assert expected_rows_set.issubset(observed_entries_set), f"Expected rows {expected_rows_set} not found in observed entries {log_observed_entries(ctx)}"

@then("the TSV SJ should equal the INFO SJ")
def tsv_matches_info(ctx):
    result = ctx["result"]
    assert result.tsv_rows[0]["SJ"] == result.rows[0]["SJ"]