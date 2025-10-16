from __future__ import annotations

import csv
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Region:
    """Represents an exon or intron within a transcript (1-based inclusive coordinates)."""

    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    number: int  # 1-based within the transcript for the region type
    region_type: str  # "exon" or "intron"

    def contains(self, pos1: int) -> bool:
        """Return True if the 1-based position falls within the region."""
        return self.start <= pos1 <= self.end


@dataclass
class Transcript:
    """Single transcript parsed from a genePred record."""

    name: str
    chrom: str
    strand: str
    gene: str
    tx_start: int
    tx_end: int
    exons: Sequence[Region]
    introns: Sequence[Region]
    mane: bool = False

    def overlaps(self, start0: int, end0: int) -> bool:
        """Return True if the transcript overlaps the interval [start0, end0)."""
        return not (end0 <= self.tx_start or start0 >= self.tx_end)

    def locate(self, pos1: int) -> Optional[Region]:
        """
        Identify the region (exon or intron) containing the 1-based position.
        """
        for region in self.exons:
            if region.contains(pos1):
                return region
        for region in self.introns:
            if region.contains(pos1):
                return region
        return None

    def donors(self) -> Sequence[int]:
        """Return donor site positions (1-based coordinate at exon boundary)."""
        donors: List[int] = []
        if self.strand == "+":
            for exon in self.exons[:-1]:
                donors.append(exon.end)
        else:
            for exon in self.exons[1:]:
                donors.append(exon.start)
        return donors

    def acceptors(self) -> Sequence[int]:
        """Return acceptor site positions (1-based coordinate at intron boundary)."""
        acceptors: List[int] = []
        if self.strand == "+":
            for exon in self.exons[1:]:
                acceptors.append(exon.start)
        else:
            for exon in self.exons[:-1]:
                acceptors.append(exon.end)
        return acceptors


def load_mane_transcripts(path: Optional[Path]) -> Tuple[set, Dict[str, str]]:
    """
    Load MANE transcript identifiers.

    Returns:
        Tuple of (set of transcript IDs, mapping from transcript ID to gene symbol).
    """
    if path is None:
        return set(), {}

    mane_ids: set = set()
    gene_map: Dict[str, str] = {}

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            refseq = row.get("RefSeq_nuc")
            ensembl = row.get("Ensembl_nuc")
            gene = row.get("symbol") or ""
            for tid in (refseq, ensembl):
                if tid:
                    mane_ids.add(tid.strip())
                    gene_map[tid.strip()] = gene
    return mane_ids, gene_map


def parse_gene_pred(path: Path, mane_ids: Optional[set] = None) -> Iterator[Transcript]:
    """
    Parse a genePred formatted file and yield Transcript objects.
    """
    mane_ids = mane_ids or set()
    required_cols = 15

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < required_cols:
                raise ValueError(f"Invalid genePred record (expected {required_cols} columns): {line}")

            name = parts[0]
            chrom = parts[1]
            strand = parts[2]
            tx_start = int(parts[3])
            tx_end = int(parts[4])
            exon_count = int(parts[7])
            exon_starts = _parse_int_list(parts[8])
            exon_ends = _parse_int_list(parts[9])
            gene = parts[11] if len(parts) > 11 else name

            if len(exon_starts) != exon_count or len(exon_ends) != exon_count:
                raise ValueError(f"Inconsistent exon count for transcript {name}")

            exons = list(zip(exon_starts, exon_ends))
            transcript_exons = _build_regions(exons, strand, region_type="exon")
            transcript_introns = _build_introns(exons, strand)

            transcript = Transcript(
                name=name,
                chrom=chrom,
                strand=strand,
                gene=gene,
                tx_start=tx_start,
                tx_end=tx_end,
                exons=transcript_exons,
                introns=transcript_introns,
                mane=name in mane_ids,
            )
            yield transcript


def _build_regions(exons: Iterable[Tuple[int, int]], strand: str, region_type: str) -> List[Region]:
    ordered = sorted(exons, key=lambda item: item[0])
    if strand == "-":
        ordered = list(reversed(ordered))

    regions: List[Region] = []
    for idx, (start, end) in enumerate(ordered, start=1):
        regions.append(Region(start=start + 1, end=end, number=idx, region_type=region_type))
    return regions


def _build_introns(exons: Iterable[Tuple[int, int]], strand: str) -> List[Region]:
    ordered = sorted(exons, key=lambda item: item[0])
    intron_pairs: List[Tuple[int, int]] = []
    for (start_a, end_a), (start_b, end_b) in zip(ordered, ordered[1:]):
        intron_pairs.append((end_a, start_b))

    if strand == "-":
        intron_pairs = list(reversed(intron_pairs))

    regions: List[Region] = []
    for idx, (start, end) in enumerate(intron_pairs, start=1):
        regions.append(Region(start=start + 1, end=end, number=idx, region_type="intron"))
    return regions


def _parse_int_list(field: str) -> List[int]:
    items = [item for item in field.strip().split(",") if item]
    return [int(item) for item in items]


class TranscriptIndex:
    """
    Provides fast lookup of transcripts overlapping a genomic interval.
    """

    def __init__(self, transcripts: Iterable[Transcript]):
        self.by_chrom: Dict[str, List[Transcript]] = {}
        for transcript in transcripts:
            self.by_chrom.setdefault(transcript.chrom, []).append(transcript)

        # Ensure deterministic ordering (MANE first, then by name).
        for chrom, items in self.by_chrom.items():
            self.by_chrom[chrom] = sorted(
                items,
                key=lambda tx: (
                    0 if tx.mane else 1,
                    tx.gene,
                    tx.name,
                ),
            )

    def fetch(self, chrom: str, start0: int, end0: int) -> List[Transcript]:
        transcripts = self.by_chrom.get(chrom, [])
        return [tx for tx in transcripts if tx.overlaps(start0, end0)]


def build_transcript_index(
    gene_pred_path: Path, chrom_mapper, mane_path: Optional[Path] = None
) -> TranscriptIndex:
    mane_ids, _ = load_mane_transcripts(mane_path)
    transcripts = []
    for transcript in parse_gene_pred(gene_pred_path, mane_ids):
        mapped = replace(transcript, chrom=chrom_mapper.to_vcf(transcript.chrom))
        transcripts.append(mapped)
    return TranscriptIndex(transcripts)
