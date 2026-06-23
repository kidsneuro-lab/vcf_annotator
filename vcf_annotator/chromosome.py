from __future__ import annotations

from typing import Iterable


UCSC = "ucsc"
PLAIN = "plain"


def to_ucsc(chrom: str) -> str:
    """
    Map a chromosome name into UCSC style for internal annotation lookups.
    """
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"MT", "M"}:
        return "chrM"
    if chrom.isdigit() or chrom in {"X", "Y"}:
        return f"chr{chrom}"
    return chrom


def from_ucsc(chrom: str, target_style: str) -> str:
    """
    Map a UCSC chromosome name into a target naming style.
    """
    if target_style == UCSC:
        return chrom
    if target_style == PLAIN and chrom.startswith("chr"):
        if chrom == "chrM":
            return "MT"
        return chrom[3:]
    return chrom


def detect_style(contigs: Iterable[str]) -> str:
    """
    Detect the dominant chromosome naming style from a set of contigs.
    """
    names = list(contigs)
    has_chr = any(name.startswith("chr") for name in names)
    return UCSC if has_chr else PLAIN


class ChromosomeMapper:
    """
    Handles mapping between UCSC ('chr1') and plain ('1') contig naming conventions.
    """

    def __init__(self, contigs: Iterable[str]):
        contigs = list(contigs)
        self.style = detect_style(contigs)
        self.vcf_has_chr = any(name.startswith("chr") for name in contigs)
        self.vcf_has_plain = any(not name.startswith("chr") for name in contigs)

    def to_ucsc(self, chrom: str) -> str:
        return to_ucsc(chrom)

    def from_ucsc(self, chrom: str) -> str:
        return from_ucsc(chrom, self.style)

    def to_vcf(self, chrom: str) -> str:
        """
        Map a chromosome name from input data to the VCF contig style.
        """
        return self.from_ucsc(to_ucsc(chrom))

    def to_external(self, chrom: str) -> str:
        """
        Map a chromosome name from internal UCSC style to this mapper's style.
        """
        return self.from_ucsc(to_ucsc(chrom))
