from __future__ import annotations

import csv
import gzip
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
            self._ensure_record_contigs(header)
            for annotator in self.annotators:
                annotator.register_fields(header)

            LOG.info("Writing annotated VCF to %s", self.config.output_path)
            with pysam.VariantFile(str(self.config.output_path), "w", header=header) as writer:
                tsv_writer = self._init_tsv_writer(self.config.tsv_path, self.annotators)

                record_count = 0
                variant_count = 0
                for record in reader:
                    record_count += 1
                    alts = list(record.alts or [])

                    if not alts:
                        writer.write(record)
                        if tsv_writer is not None:
                            row = self._build_tsv_row(
                                record,
                                ".",
                                {},
                                getattr(tsv_writer, "_extra_fields", []),
                            )
                            tsv_writer.writerow(row)
                        continue

                    per_alt_results = self._annotate_record(record, alts)

                    for alt_index, alt in enumerate(alts):
                        aggregate = per_alt_results[alt_index]
                        row_infos = aggregate.rows or [dict()]

                        for info_map in row_infos:
                            new_record = self._create_single_alt_record(
                                writer,
                                record,
                                alt,
                                alt_index,
                                alts,
                            )
                            for key, value in info_map.items():
                                new_record.info[key] = value
                            writer.write(new_record)

                            if tsv_writer is not None:
                                row = self._build_tsv_row(
                                    record,
                                    alt,
                                    info_map,
                                    getattr(tsv_writer, "_extra_fields", []),
                                )
                                tsv_writer.writerow(row)

                        variant_count += len(row_infos)

                    if record_count % self.config.log_every == 0:
                        LOG.info("Processed %d VCF records (%d annotated rows)", record_count, variant_count)

        close_writer(tsv_writer)

        for annotator in self.annotators:
            annotator.close()

        LOG.info("Completed annotation of %d records.", record_count)

    def _ensure_record_contigs(self, header) -> None:
        if header.contigs:
            return

        contigs = _collect_record_contigs(self.config.input_path)
        for contig in contigs:
            header.contigs.add(contig)

        if contigs:
            LOG.info(
                "Input VCF header has no contig definitions; added %d contigs from records.",
                len(contigs),
            )

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

    def _build_tsv_row(
        self,
        record,
        alt: str,
        info_map,
        extra_fields: List[str],
    ) -> Dict[str, str]:
        row = {
            "CHROM": record.chrom,
            "POS": str(record.pos),
            "ID": record.id or ".",
            "REF": record.ref,
            "ALT": alt,
        }

        for field in extra_fields:
            value = info_map.get(field, "NA") if info_map else "NA"
            row[field] = str(value)

        return row

    def _create_single_alt_record(
        self,
        writer,
        base_record,
        alt: str,
        alt_index: int,
        alts: List[str],
    ):
        new_record = writer.new_record(
            contig=base_record.chrom,
            start=base_record.start,
            stop=base_record.stop,
            alleles=(base_record.ref, alt),
        )

        new_record.id = base_record.id
        new_record.qual = base_record.qual
        new_record.pos = base_record.pos

        new_record.filter.clear()
        for filter_id in base_record.filter.keys():
            new_record.filter.add(filter_id)

        original_header = base_record.header

        for key in base_record.info.keys():
            info_def = original_header.info.get(key)
            value = base_record.info[key]

            if info_def is not None and info_def.number == 0:
                if value:
                    new_record.info[key] = True
                continue

            selected = self._extract_info_value(value, info_def, alt_index)
            if selected is not None:
                new_record.info[key] = selected

        for sample in base_record.samples:
            new_record.samples[sample] = base_record.samples[sample]

        return new_record

    def _extract_info_value(self, value, info_def, alt_index: int):
        if value is None:
            return None

        if info_def is None:
            return value

        number = info_def.number

        if number == "A":
            if isinstance(value, tuple):
                if len(value) > alt_index:
                    return value[alt_index]
                return value[-1] if value else None
            return value

        if number == "R":
            if isinstance(value, tuple):
                idx = alt_index + 1
                if len(value) > idx:
                    return value[idx]
                return value[-1] if value else None
            return value

        if number == "G":
            return value

        if number == ".":
            return value

        if isinstance(number, int):
            if number == 0:
                return None
            if number == 1:
                if isinstance(value, tuple):
                    return value[0] if value else None
                return value
            return value

        return value

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
        writer._extra_fields = columns[5:]  # type: ignore[attr-defined]
        return writer


def _collect_record_contigs(path: Path) -> List[str]:
    seen = set()
    contigs: List[str] = []

    with Path(path).open("rb") as raw_handle:
        magic = raw_handle.read(2)

    opener = gzip.open if magic == b"\x1f\x8b" else open
    with opener(path, "rt", newline="") as handle:
        for line in handle:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#"):
                continue

            contig = line.split("\t", 1)[0]
            if contig and contig not in seen:
                seen.add(contig)
                contigs.append(contig)

    return contigs


def close_writer(writer: Optional[csv.DictWriter]) -> None:
    if writer is None:
        return
    handle = getattr(writer, "_handle", None)
    if handle:
        handle.close()
