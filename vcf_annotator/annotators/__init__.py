"""
Annotator implementations for vcf_annotator.
"""

from .base import Annotator, AnnotationResult  # noqa: F401
from .splice_distance import SpliceJunctionDistanceAnnotator  # noqa: F401
from .custom_vcf import CustomVcfAnnotator  # noqa: F401

__all__ = [
    "Annotator",
    "AnnotationResult",
    "SpliceJunctionDistanceAnnotator",
    "CustomVcfAnnotator",
]

