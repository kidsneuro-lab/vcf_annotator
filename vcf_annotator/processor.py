from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pysam

from .annotators.base import AnnotationResult, VariantContext

LOG = logging.getLogger(__name__)


@dataclass
class ProcessorConfig:
    input_path: Path
    output_path: Path
    annotators: Iterable
    tsv_path: Optional[Path] = None
    log_every: int = 10_000


class VariantProcessor:
    """
    Coordinates iteration over the input VCF and applies annotators per ALT allele.
    """

    def __init__(self, config: ProcessorConfig):
        self.config = config
        self.annotators = list(config.annotators)

    def run(self) -> None:
        LOG.info("Opening input VCF: %s", self.config.input_path)
        with pysam.VariantFile(str(self.config.input_path)) as reader:
            header = reader.header.copy()
            for annotator in self.annotators:
                annotator.register_fields(header)

            LOG.info("Writing annotated VCF to %s", self.config.output_path)
            with pysam.VariantFile(str(self.config.output_path), "w", header=header) as writer:
                tsv_writer = self._init_tsv_writer(self.config.tsv_path, self.annotators)

                record_count = 0
                variant_count = 0
                for record in reader:
                    record_count += 1
                    alts = record.alts or []
                    if not alts:
                        writer.write(record)
                        continue

                    per_alt_results = self._annotate_record(record, alts)
                    per_field_values = self._collect_per_field(alts, per_alt_results)

                    for key, values in per_field_values.items():
                        record.info[key] = tuple(values)

                    writer.write(record)
                    if tsv_writer is not None:
                        rows = self._build_tsv_rows(record, alts, per_field_values)
                        for row in rows:
                            tsv_writer.writerow(row)
                    variant_count += len(alts)

                    if record_count % self.config.log_every == 0:
                        LOG.info("Processed %d VCF records (%d ALT alleles)", record_count, variant_count)

        close_writer(tsv_writer)

        for annotator in self.annotators:
            annotator.close()

        LOG.info("Completed annotation of %d records.", record_count)

    def _annotate_record(self, record, alts: List[str]) -> List[AnnotationResult]:
        results: List[AnnotationResult] = []
        for alt_index, alt in enumerate(alts):
            context = VariantContext(record, alt, alt_index)
            aggregate = AnnotationResult()
            for annotator in self.annotators:
                res = annotator.annotate(context)
                aggregate.merge(res)
            results.append(aggregate)
        return results

    def _collect_per_field(self, alts: List[str], per_alt_results: List[AnnotationResult]) -> Dict[str, List[str]]:
        num_alts = len(alts)
        by_field: Dict[str, List[str]] = {}

        for alt_index, result in enumerate(per_alt_results):
            for key, value in result.info.items():
                by_field.setdefault(key, ["NA"] * num_alts)
                by_field[key][alt_index] = _as_string(value)

        return by_field

    def _build_tsv_rows(self, record, alts: List[str], per_field_values: Dict[str, List[str]]) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        extra_fields = set(per_field_values.keys())

        for alt_index, alt in enumerate(alts):
            row = {
                "CHROM": record.chrom,
                "POS": str(record.pos),
                "ID": record.id or ".",
                "REF": record.ref,
                "ALT": alt,
            }

            for field in extra_fields:
                values = per_field_values.get(field, ["NA"] * len(alts))
                row[field] = values[alt_index] if alt_index < len(values) else "NA"

            rows.append(row)

        return rows

    def _init_tsv_writer(self, path: Optional[Path], annotators: Iterable) -> Optional[csv.DictWriter]:
        if path is None:
            return None

        extra_fields: List[str] = []
        for annotator in annotators:
            if hasattr(annotator, "output_fields"):
                extra_fields.extend(list(annotator.output_fields()))

        columns = ["CHROM", "POS", "ID", "REF", "ALT"] + sorted(set(extra_fields))
        handle = Path(path).open("w", newline="")
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", restval="NA")
        writer.writeheader()
        writer._handle = handle  # type: ignore[attr-defined]
        writer._columns = columns  # type: ignore[attr-defined]
        return writer


def _as_string(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (list, tuple)):
        return "|".join(str(item) for item in value)
    return str(value)


def close_writer(writer: Optional[csv.DictWriter]) -> None:
    if writer is None:
        return
    handle = getattr(writer, "_handle", None)
    if handle:
        handle.close()
