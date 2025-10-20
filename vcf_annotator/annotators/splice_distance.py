from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..variant_utils import VariantCoordinates, compute_variant_bounds
from .base import AnnotationResult, Annotator, VariantContext
from ..transcripts import Region


@dataclass
class DistanceResult:
    distance: Optional[int]
    region_type: str
    region_number: str


class SpliceJunctionDistanceAnnotator(Annotator):
    """
    Annotates distance to nearest splice donor/acceptor for overlapping transcripts.
    """

    def __init__(self, transcript_index, prefix: str, include_mane: bool = False):
        self.transcripts = transcript_index
        self.prefix = prefix.upper()
        self.include_mane = include_mane
        self.info_id = self.prefix

    def register_fields(self, header) -> None:
        description = self._info_description()
        header.add_line(
            f"##INFO=<ID={self.info_id},Number=.,Type=String,Description=\"{description}\">"
        )

    def output_fields(self):
        return [self.info_id]

    def annotate(self, context: VariantContext) -> AnnotationResult:
        variant_type = context.variant_type()
        coords = compute_variant_bounds(context.pos, context.ref, context.alt)

        transcripts = self.transcripts.fetch(context.chrom, coords.start0, coords.end0)
        result = AnnotationResult()
        entries = []

        if not transcripts:
            entries.append(
                self._build_annotation_entry(
                    allele=context.alt,
                    transcript_name="NA",
                    gene="NA",
                    variant_type=variant_type,
                    donor_result=DistanceResult(distance=None, region_type="NA", region_number="NA"),
                    acceptor_result=DistanceResult(distance=None, region_type="NA", region_number="NA"),
                    mane_flag="0" if self.include_mane else None,
                )
            )
        else:
            for transcript in transcripts:
                donor_result = self._distance_to_site(
                    transcript, coords, is_donor=True
                )
                acceptor_result = self._distance_to_site(
                    transcript, coords, is_donor=False
                )

                mane_flag = None
                if self.include_mane:
                    mane_flag = "1" if getattr(transcript, "mane", False) else "0"

                entries.append(
                    self._build_annotation_entry(
                        allele=context.alt,
                        transcript_name=transcript.name,
                        gene=transcript.gene,
                        variant_type=variant_type,
                        donor_result=donor_result,
                        acceptor_result=acceptor_result,
                        mane_flag=mane_flag,
                    )
                )

        value = ",".join(entries) if entries else "NA"
        row = {self.info_id: value}
        result.rows.append(row)
        result.tsv_rows.append(row.copy())
        return result

    def _info_description(self) -> str:
        fields = [
            "Allele",
            "Transcript",
            "Gene",
            "VariantType",
            "DonorDist",
            "DonorRegionType",
            "DonorRegionNo",
            "AcceptorDist",
            "AcceptorRegionType",
            "AcceptorRegionNo",
        ]
        if self.include_mane:
            fields.append("MANE")
        format_description = "|".join(fields)
        return f"Splice annotations per transcript. Format: {format_description}"

    def _build_annotation_entry(
        self,
        *,
        allele: str,
        transcript_name: str,
        gene: str,
        variant_type: str,
        donor_result: DistanceResult,
        acceptor_result: DistanceResult,
        mane_flag: Optional[str],
    ) -> str:
        parts = [
            allele,
            transcript_name,
            gene,
            variant_type,
            _format_distance(donor_result.distance),
            donor_result.region_type,
            donor_result.region_number,
            _format_distance(acceptor_result.distance),
            acceptor_result.region_type,
            acceptor_result.region_number,
        ]
        if mane_flag is not None:
            parts.append(mane_flag)
        return "|".join(parts)

    def _distance_to_site(
        self, transcript, coords: VariantCoordinates, is_donor: bool
    ) -> DistanceResult:
        anchors = self._anchors(transcript, coords)
        best_distance = None
        best_region: Optional[Region] = None

        for anchor_pos, region in anchors:
            if anchor_pos is None or region is None:
                continue

            if is_donor:
                distance = self._donor_distance(transcript, anchor_pos, region)
            else:
                distance = self._acceptor_distance(transcript, anchor_pos, region)

            if distance is None:
                continue

            if best_distance is None or abs(distance) < abs(best_distance):
                best_distance = distance
                best_region = region
            elif best_distance is not None and abs(distance) == abs(best_distance):
                # Prefer intronic assignments when equidistant to align with boundary rules.
                if best_region and best_region.region_type == "exon" and region.region_type == "intron":
                    best_distance = distance
                    best_region = region

        if best_distance is None:
            return DistanceResult(distance=None, region_type="NA", region_number="NA")

        region_type = best_region.region_type if best_region else "NA"
        region_number = str(best_region.number) if best_region else "NA"
        return DistanceResult(distance=best_distance, region_type=region_type, region_number=region_number)

    def _anchors(self, transcript, coords: VariantCoordinates):
        start_anchor = coords.start0 + 1
        end_anchor0 = max(coords.end0 - 1, coords.start0)
        end_anchor = end_anchor0 + 1

        start_region = transcript.locate(start_anchor)
        end_region = transcript.locate(end_anchor)

        return (
            (start_anchor, start_region),
            (end_anchor, end_region),
        )

    def _donor_distance(self, transcript, pos1: int, region: Region) -> Optional[int]:
        if region.region_type == "intron":
            if transcript.strand == "+":
                return pos1 - region.start + 1
            return region.end - pos1 + 1

        # Exonic region
        intron = self._intron_after_exon(transcript, region.number)
        if intron is None:
            return None

        if transcript.strand == "+":
            return -(region.end - pos1 + 1)
        return -(pos1 - region.start + 1)

    def _acceptor_distance(self, transcript, pos1: int, region: Region) -> Optional[int]:
        if region.region_type == "intron":
            if transcript.strand == "+":
                return -(region.end - pos1 + 1)
            return -(pos1 - region.start + 1)

        intron = self._intron_before_exon(transcript, region.number)
        if intron is None:
            return None

        if transcript.strand == "+":
            return pos1 - region.start + 1
        return region.end - pos1 + 1

    @staticmethod
    def _intron_after_exon(transcript, exon_number: int) -> Optional[Region]:
        if exon_number > len(transcript.introns):
            return None
        return transcript.introns[exon_number - 1]

    @staticmethod
    def _intron_before_exon(transcript, exon_number: int) -> Optional[Region]:
        if exon_number <= 1:
            return None
        return transcript.introns[exon_number - 2]


def _format_distance(value: Optional[int]) -> str:
    return str(value) if value is not None else "NA"
