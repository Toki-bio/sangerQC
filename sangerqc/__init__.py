"""sangerQC: chromatogram quality from raw traces, not Phred scores alone."""

from .v0_1 import classify as classify_v01
from .v0_2 import ALPHA, BETA, classify as classify_v02

__all__ = ["classify_v01", "classify_v02", "ALPHA", "BETA"]
