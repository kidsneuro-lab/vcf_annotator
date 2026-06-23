from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Mapping, MutableMapping, Optional

from ..chromosome import to_ucsc


@dataclass
class AnnotationResult:
    """Container for annotator output grouped into per-row dictionaries."""

    rows: List[MutableMapping[str, object]] = field(default_factory=list)
    tsv_rows: List[Mapping[str, object]] = field(default_factory=list)

    def merge(self, other: "AnnotationResult") -> None:
        """Merge another result, broadcasting single rows where needed."""

        other_rows = [dict(row) for row in (other.rows or [{}])]

        if not self.rows:
            self.rows = other_rows
        else:
            if len(other_rows) == len(self.rows):
                for idx, row in enumerate(other_rows):
                    self.rows[idx].update(row)
            elif len(other_rows) == 1:
                filler = other_rows[0]
                for row in self.rows:
                    row.update(filler)
            elif len(self.rows) == 1:
                base = dict(self.rows[0])
                self.rows = [dict(base) for _ in range(len(other_rows))]
                for idx, row in enumerate(other_rows):
                    self.rows[idx].update(row)
            else:
                raise ValueError("Incompatible annotation row counts during merge.")

        if other.tsv_rows:
            self.tsv_rows.extend(other.tsv_rows)


class Annotator(ABC):
    """Base class for annotator implementations."""

    @abstractmethod
    def register_fields(self, header) -> None:
        """
        Update the pysam header with INFO field definitions required by the annotator.

        Args:
            header: pysam.VariantHeader object.
        """

    @abstractmethod
    def annotate(self, context: "VariantContext") -> AnnotationResult:
        """
        Produce annotations for the provided variant context.

        Args:
            context: Single-allele variant context.

        Returns:
            AnnotationResult containing INFO updates and optional TSV rows.
        """

    def close(self) -> None:  # pragma: no cover - hook for subclasses.
        """Hook for releasing resources."""
        return

    def output_fields(self):
        """
        Return INFO field names produced by this annotator.
        """
        return []


class VariantContext:
    """
    Lightweight wrapper representing a single-allele variant.

    Attributes:
        record: pysam.VariantRecord from the input VCF.
        alt: Alternate allele string for this context.
        alt_index: Index of the alternate allele within the original record.
        original_chrom: Chromosome/contig name from the input VCF.
        chrom: Canonical UCSC chromosome/contig name used for annotation lookups.
        ucsc_chrom: Alias for the canonical UCSC chromosome/contig name.
        pos: 1-based genomic position.
        ref: Reference allele string.
    """

    __slots__ = ("record", "alt", "alt_index", "original_chrom", "chrom", "ucsc_chrom", "pos", "ref")

    def __init__(self, record, alt: str, alt_index: int, chrom: Optional[str] = None):
        self.record = record
        self.alt = alt
        self.alt_index = alt_index
        self.original_chrom = record.chrom
        self.ucsc_chrom = chrom or to_ucsc(record.chrom)
        self.chrom = self.ucsc_chrom
        self.pos = record.pos
        self.ref = record.ref

    @property
    def id(self) -> Optional[str]:
        """Return the variant identifier."""
        return self.record.id

    def variant_type(self) -> str:
        """
        Classify the variant based on REF/ALT lengths.

        Returns:
            One of 'snp', 'ins', 'del', 'delins'.
        """
        ref_len = len(self.ref)
        alt_len = len(self.alt)

        if ref_len == 1 and alt_len == 1:
            return "snp"
        if ref_len == 1 and alt_len > 1:
            return "ins"
        if ref_len > 1 and alt_len == 1:
            return "del"
        return "delins"
