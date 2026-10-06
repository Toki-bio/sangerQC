"""sangerQC: chromatogram quality from raw traces, not Phred scores alone."""

__version__ = "0.3.0"

from .v0_1 import classify as classify_v01
from .v0_2 import ALPHA, BETA, classify as classify_v02
from .v0_3 import classify as classify_v03
from .export import to_sequence

__all__ = ["classify_v01", "classify_v02", "classify_v03", "to_sequence", "ALPHA", "BETA"]
