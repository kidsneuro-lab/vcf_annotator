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
            f"##INFO=<ID={info['transcript']},Number=A,Type=String,Description=\"{self.prefix}: Transcript identifier for this transcript-specific record.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['gene']},Number=A,Type=String,Description=\"{self.prefix}: Gene symbol for the annotated transcript.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['variant_type']},Number=A,Type=String,Description=\"{self.prefix}: Variant type per ALT classified as snp/ins/del/delins.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['ddon']},Number=A,Type=String,Description=\"{self.prefix}: Distance to the nearest splice donor for this transcript (NA if unavailable).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['ddon_region_type']},Number=A,Type=String,Description=\"{self.prefix}: Region type (exon/intron/NA) used for donor distance.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['ddon_region_no']},Number=A,Type=String,Description=\"{self.prefix}: Region number used for donor distance.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['dacc']},Number=A,Type=String,Description=\"{self.prefix}: Distance to the nearest splice acceptor for this transcript (NA if unavailable).\">"
        )
        header.add_line(
            f"##INFO=<ID={info['dacc_region_type']},Number=A,Type=String,Description=\"{self.prefix}: Region type (exon/intron/NA) used for acceptor distance.\">"
        )
        header.add_line(
            f"##INFO=<ID={info['dacc_region_no']},Number=A,Type=String,Description=\"{self.prefix}: Region number used for acceptor distance.\">"
        )

    def output_fields(self):
        return list(self.info_fields.values())

    def annotate(self, context: VariantContext) -> AnnotationResult:
        variant_type = context.variant_type()
        coords = compute_variant_bounds(context.pos, context.ref, context.alt)

        transcripts = self.transcripts.fetch(context.chrom, coords.start0, coords.end0)
        result = AnnotationResult()

        info = self.info_fields

        if not transcripts:
            row = {
                info["variant_type"]: variant_type,
                info["transcript"]: "NA",
                info["gene"]: "NA",
                info["ddon"]: "NA",
                info["ddon_region_type"]: "NA",
                info["ddon_region_no"]: "NA",
                info["dacc"]: "NA",
                info["dacc_region_type"]: "NA",
                info["dacc_region_no"]: "NA",
            }
            result.rows.append(row)
            result.tsv_rows.append(row.copy())
            return result

        for transcript in transcripts:
            donor_result = self._distance_to_site(
                transcript, coords, is_donor=True
            )
            acceptor_result = self._distance_to_site(
                transcript, coords, is_donor=False
            )

            row = {
                info["variant_type"]: variant_type,
                info["transcript"]: transcript.name,
                info["gene"]: transcript.gene,
                info["ddon"]: _format_distance(donor_result.distance),
                info["ddon_region_type"]: donor_result.region_type,
                info["ddon_region_no"]: donor_result.region_number,
                info["dacc"]: _format_distance(acceptor_result.distance),
                info["dacc_region_type"]: acceptor_result.region_type,
                info["dacc_region_no"]: acceptor_result.region_number,
            }

            result.rows.append(row)
            result.tsv_rows.append(row.copy())

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
