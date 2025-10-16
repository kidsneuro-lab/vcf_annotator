from __future__ import annotations

from typing import Iterable


class ChromosomeMapper:
    """
    Handles mapping between UCSC ('chr1') and Ensembl ('1') contig naming conventions.
    """

    def __init__(self, contigs: Iterable[str]):
        contigs = list(contigs)
        self.vcf_has_chr = any(name.startswith("chr") for name in contigs)
        self.vcf_has_plain = any(not name.startswith("chr") for name in contigs)

    def to_vcf(self, chrom: str) -> str:
        """
        Map a chromosome name from input data to the VCF contig style.
        """
        if self.vcf_has_chr and not chrom.startswith("chr"):
            if chrom in {"MT", "M"}:
                return "chrM"
            return f"chr{chrom}"
        if not self.vcf_has_chr and chrom.startswith("chr"):
            if chrom == "chrM":
                return "MT"
            return chrom[3:]
        return chrom

    def to_external(self, chrom: str) -> str:
        """
        Map a chromosome name from the VCF to external data style, assuming
        the external data uses the opposite convention.
        """
        if self.vcf_has_chr:
            if chrom.startswith("chr"):
                return chrom[3:] if chrom != "chrM" else "MT"
            return chrom
        return chrom if not chrom.startswith("chr") else chrom

