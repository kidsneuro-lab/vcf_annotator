from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class VariantCoordinates:
    """
    Represents variant genomic coordinates in 0-based half-open form.
    """

    start0: int
    end0: int


def compute_variant_bounds(pos: int, ref: str, alt: str) -> VariantCoordinates:
    """
    Compute 0-based half-open coordinates for the variant based on REF/ALT.

    For insertions, the start coordinate corresponds to the base preceding the
    insertion (consistent with VCF conventions) and end equals start + 1.
    """
    start0 = pos - 1
    ref_len = len(ref)
    alt_len = len(alt)

    if ref_len == alt_len == 1:
        return VariantCoordinates(start0=start0, end0=start0 + 1)

    if ref_len == 1 and alt_len > 1:
        # insertion occurs between positions start and start+1
        return VariantCoordinates(start0=start0, end0=start0 + 1)

    # deletions and complex alleles cover the reference span
    end0 = start0 + ref_len
    return VariantCoordinates(start0=start0, end0=end0)


def distance_to_site(position0: int, site0: int) -> int:
    """
    Compute absolute distance between a genomic position and a splice site.
    """
    return abs(position0 - site0)


def closest_site(position_bounds: Tuple[int, int], candidates: Tuple[int, ...]) -> Tuple[int, int]:
    """
    Find the closest site to either the start or end of a variant interval.

    Args:
        position_bounds: Tuple of (start0, end0) for the variant.
        candidates: Tuple of site coordinates.

    Returns:
        Tuple of (distance, site_position). If no candidates, returns (None, None).
    """
    if not candidates:
        return None, None  # type: ignore[return-value]

    start0, end0 = position_bounds
    # Evaluate distance from both the start and end coordinates.
    best_distance = None
    best_site = None

    for site in candidates:
        for anchor in (start0, end0 - 1):
            distance = abs(anchor - site)
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_site = site
    return best_distance, best_site  # type: ignore[return-value]

