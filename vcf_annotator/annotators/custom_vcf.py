from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pysam

from .base import AnnotationResult, Annotator, VariantContext


class CustomVcfAnnotator(Annotator):
    """
    Annotates variants by looking up matching records in an auxiliary VCF.
    """

    def __init__(self, vcf_path: Path, prefix: str, fields: List[str]):
        self.vcf_path = Path(vcf_path)
        self.prefix = prefix.upper()
        self.fields = [field.strip() for field in fields if field.strip()]

        if not self.fields:
            raise ValueError("At least one INFO field must be provided for custom VCF annotation.")

        self._external = pysam.VariantFile(str(self.vcf_path))
        self._external_has_chr = any(contig.startswith("chr") for contig in self._external.header.contigs)

    def register_fields(self, header) -> None:
        for field in self.fields:
            field_id = self._info_name(field)
            header.add_line(
                f"##INFO=<ID={field_id},Number=.,Type=String,Description=\"{self.prefix}: Values from {self.vcf_path.name} INFO field {field} for matching records.\">"
            )

    def annotate(self, context: VariantContext) -> AnnotationResult:
        chrom = self._map_chrom(context.chrom)

        try:
            records = self._external.fetch(chrom, context.pos - 1, context.pos + len(context.ref))
        except ValueError:
            records = self._iter_all()

        collected: Dict[str, List[str]] = {field: [] for field in self.fields}
        tsv_rows = []

        for record in records:
            if record.pos != context.pos:
                continue
            if record.ref != context.ref:
                continue
            if context.alt not in (record.alts or []):
                continue

            for field in self.fields:
                value = record.info.get(field)
                value_str = _stringify(value)
                if value_str is not None:
                    collected[field].append(value_str)

            row = {self._info_name(field): collected[field][-1] if collected[field] else "NA" for field in self.fields}
            tsv_rows.append(row)

        result = AnnotationResult()
        for field, values in collected.items():
            key = self._info_name(field)
            result.info[key] = "|".join(values) if values else "NA"

        if tsv_rows:
            result.tsv_rows.extend(tsv_rows)

        return result

    def close(self) -> None:
        self._external.close()

    def output_fields(self):
        return [self._info_name(field) for field in self.fields]

    def _map_chrom(self, chrom: str) -> str:
        if self._external_has_chr and not chrom.startswith("chr"):
            return "chrM" if chrom in {"MT", "M"} else f"chr{chrom}"
        if not self._external_has_chr and chrom.startswith("chr"):
            return "MT" if chrom == "chrM" else chrom[3:]
        return chrom

    def _info_name(self, field: str) -> str:
        return f"{self.prefix}_{field}"

    def _iter_all(self):
        with pysam.VariantFile(str(self.vcf_path)) as handle:
            for record in handle:
                yield record


def _stringify(value):
    if value is None:
        return None
    if isinstance(value, tuple):
        return ",".join(str(item) for item in value)
    return str(value)
