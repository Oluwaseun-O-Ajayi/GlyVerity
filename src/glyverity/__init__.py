"""GlyVerity: artifact-aware, FDR-controlled glycan composition assignment.

The package is organised as one engine with pluggable evidence layers:

* :mod:`glyverity.chemistry` holds monosaccharide, label and adduct masses.
* :mod:`glyverity.evidence` defines the evidence-layer plugin interface.
"""

__version__ = "0.0.1"
