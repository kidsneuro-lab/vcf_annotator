from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from ..variant_utils import VariantCoordinates, compute_variant_bounds
from .base import AnnotationResult, Annotator, VariantContext


@dataclass
class DistanceResult:
    distance: Optional[int]
    region_type: str
    region_number: str


class SpliceJunctionDistanceAnnotator(Annotator):
    """
    Annotates distance to nearest splice donor/acceptor for overlapping transcripts.
    """

    def __init__(self, transcript_index, prefix: str):
        self.transcripts = transcript_index
        self.prefix = prefix.upper()
        self.info_fields = self._build_info_names()

    def _build_info_names(self) -> Dict[str, str]:
        base = {
            "transcript": "TRANSCRIPT",
            "gene": "GENE",
            "variant_type": "VARIANT_TYPE",
            "ddon": "DDON",
            "ddon_region_type": "DDON_REGION_TYPE",
            "ddon_region_no": "DDON_REGION_NO",
            "dacc": "DACC",
            "dacc_region_type": "DACC_REGION_TYPE",
            "dacc_region_no": "DACC_REGION_NO",
        }
        return {key: f"{self.prefix}_{suffix}" for key, suffix in base.items()}

    def register_fields(self, header) -> None:
        info = self.info_fields
        header.add_line(
            f"##INFO=<ID={info['transcript']},Number=A,Type=String,Description=\"{self.prefix}: Transcript identifiers overlapping each ALT (pipe separated).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['gene']},Number=A,Type=String,Description=\"{self.prefix}: Gene symbols overlapping each ALT (pipe separated).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['variant_type']},Number=A,Type=String,Description=\"{self.prefix}: Variant type per ALT classified as snp/ins/del/delins.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['ddon']},Number=A,Type=String,Description=\"{self.prefix}: Distances to nearest splice donor per ALT (pipe separated per transcript, NA if unavailable).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['ddon_region_type']},Number=A,Type=String,Description=\"{self.prefix}: Region type (exon/intron/NA) used for donor distance per ALT (pipe separated).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['ddon_region_no']},Number=A,Type=String,Description=\"{self.prefix}: Region number used for donor distance per ALT (pipe separated).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['dacc']},Number=A,Type=String,Description=\"{self.prefix}: Distances to nearest splice acceptor per ALT (pipe separated per transcript, NA if unavailable).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['dacc_region_type']},Number=A,Type=String,Description=\"{self.prefix}: Region type (exon/intron/NA) used for acceptor distance per ALT (pipe separated).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['dacc_region_no']},Number=A,Type=String,Description=\"{self.prefix}: Region number used for acceptor distance per ALT (pipe separated).\">"
        )

    def output_fields(self):
        return list(self.info_fields.values())

    def annotate(self, context: VariantContext) -> AnnotationResult:
        variant_type = context.variant_type()
        coords = compute_variant_bounds(context.pos, context.ref, context.alt)

        transcripts = self.transcripts.fetch(context.chrom, coords.start0, coords.end0)
        result = AnnotationResult()

        if not transcripts:
            for key in (
                "transcript",
                "gene",
                "ddon",
                "ddon_region_type",
                "ddon_region_no",
                "dacc",
                "dacc_region_type",
                "dacc_region_no",
            ):
                result.info[self.info_fields[key]] = "NA"
            result.info[self.info_fields["variant_type"]] = variant_type
            return result

        transcript_ids: List[str] = []
        gene_names: List[str] = []
        ddon_values: List[str] = []
        ddon_region_types: List[str] = []
        ddon_region_numbers: List[str] = []
        dacc_values: List[str] = []
        dacc_region_types: List[str] = []
        dacc_region_numbers: List[str] = []

        for transcript in transcripts:
            transcript_ids.append(transcript.name)
            gene_names.append(transcript.gene)

            donor_result = self._distance_to_site(
                transcript, coords, is_donor=True
            )
            acceptor_result = self._distance_to_site(
                transcript, coords, is_donor=False
            )

            ddon_values.append(_format_distance(donor_result.distance))
            ddon_region_types.append(donor_result.region_type)
            ddon_region_numbers.append(donor_result.region_number)
            dacc_values.append(_format_distance(acceptor_result.distance))
            dacc_region_types.append(acceptor_result.region_type)
            dacc_region_numbers.append(acceptor_result.region_number)

            # TSV row
            result.tsv_rows.append(
                {
                    "transcript": transcript.name,
                    "gene": transcript.gene,
                    "variant_type": variant_type,
                    "ddon": ddon_values[-1],
                    "ddon_region_type": ddon_region_types[-1],
                    "ddon_region_no": ddon_region_numbers[-1],
                    "dacc": dacc_values[-1],
                    "dacc_region_type": dacc_region_types[-1],
                    "dacc_region_no": dacc_region_numbers[-1],
                }
            )

        info = self.info_fields
        result.info[info["variant_type"]] = variant_type
        result.info[info["transcript"]] = "|".join(transcript_ids) if transcript_ids else "NA"
        result.info[info["gene"]] = "|".join(gene_names) if gene_names else "NA"
        result.info[info["ddon"]] = "|".join(ddon_values) if ddon_values else "NA"
        result.info[info["ddon_region_type"]] = "|".join(ddon_region_types) if ddon_region_types else "NA"
        result.info[info["ddon_region_no"]] = "|".join(ddon_region_numbers) if ddon_region_numbers else "NA"
        result.info[info["dacc"]] = "|".join(dacc_values) if dacc_values else "NA"
        result.info[info["dacc_region_type"]] = "|".join(dacc_region_types) if dacc_region_types else "NA"
        result.info[info["dacc_region_no"]] = "|".join(dacc_region_numbers) if dacc_region_numbers else "NA"
        return result

    def _distance_to_site(
        self, transcript, coords: VariantCoordinates, is_donor: bool
    ) -> DistanceResult:
        sites: Sequence[int]
        if is_donor:
            sites = tuple(transcript.donors())
        else:
            sites = tuple(transcript.acceptors())

        if not sites:
            return DistanceResult(distance=None, region_type="NA", region_number="NA")

        anchors = self._anchors(transcript, coords)
        best_distance = None
        best_region: Optional[Tuple[str, str]] = None

        for site in sites:
            for anchor_pos, region in anchors:
                if anchor_pos is None:
                    continue
                distance = abs(anchor_pos - site)
                if best_distance is None or distance < best_distance:
                    best_distance = distance
                    best_region = (
                        region.region_type if region else "NA",
                        str(region.number) if region else "NA",
                    )

        if best_distance is None:
            return DistanceResult(distance=None, region_type="NA", region_number="NA")

        region_type, region_number = best_region or ("NA", "NA")
        return DistanceResult(distance=best_distance, region_type=region_type, region_number=region_number)

    def _anchors(self, transcript, coords: VariantCoordinates):
        start_anchor = coords.start0
        end_anchor = max(coords.end0 - 1, coords.start0)

        start_region = transcript.locate(start_anchor)
        end_region = transcript.locate(end_anchor)

        return (
            (start_anchor, start_region),
            (end_anchor, end_region),
        )


def _format_distance(value: Optional[int]) -> str:
    return str(value) if value is not None else "NA"
