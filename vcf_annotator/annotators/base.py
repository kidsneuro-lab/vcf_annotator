from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional


@dataclass
class AnnotationResult:
    """
    Container for results produced by an annotator.

    Attributes:
        info: Mapping of INFO field names to values to be merged into the
            variant record. Values can be scalars or iterables. Lists are
            converted to tuples when applied to pysam records.
        tsv_rows: Optional collection of per-variant rows for TSV output.
            Each row is a mapping of column name to value.
    """

    info: MutableMapping[str, object] = field(default_factory=dict)
    tsv_rows: List[Mapping[str, object]] = field(default_factory=list)

    def merge(self, other: "AnnotationResult") -> None:
        """Merge another result into this one."""
        for key, value in other.info.items():
            if key not in self.info:
                self.info[key] = value
                continue

            existing = self.info[key]
            merged = _merge_values(existing, value)
            self.info[key] = merged

        self.tsv_rows.extend(other.tsv_rows)


def _merge_values(existing: object, new_value: object) -> object:
    """
    Merge two INFO values, preserving order and deduplicating where appropriate.
    """
    if existing is None:
        return new_value
    if new_value is None:
        return existing

    if isinstance(existing, (list, tuple)):
        existing_iter = list(existing)
    else:
        existing_iter = [existing]

    if isinstance(new_value, (list, tuple)):
        new_iter: Iterable[object] = new_value
    else:
        new_iter = [new_value]

    combined = existing_iter + list(new_iter)
    return tuple(combined)


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
        chrom: Chromosome/contig name.
        pos: 1-based genomic position.
        ref: Reference allele string.
    """

    __slots__ = ("record", "alt", "alt_index", "chrom", "pos", "ref")

    def __init__(self, record, alt: str, alt_index: int):
        self.record = record
        self.alt = alt
        self.alt_index = alt_index
        self.chrom = record.chrom
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
