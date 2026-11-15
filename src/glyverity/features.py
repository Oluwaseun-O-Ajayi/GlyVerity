"""Feature-first MS1 processing with OpenMS.

Centroided MS1 scans are placed in an OpenMS ``MSExperiment`` and processed with the OpenMS
metabolomics feature finder: MassTraceDetection, ElutionPeakDetection and FeatureFindingMetabo
(Kenar et al., Mol. Cell. Proteomics 2014). Each feature is one isotope envelope with a charge,
monoisotopic m/z, apex retention time and intensity. pyOpenMS is the optional extra ``ms``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from glyverity import chemistry as ch


@dataclass(frozen=True)
class MSFeature:
    """One isotope-envelope feature."""

    feature_id: int
    mz: float
    charge: int
    rt: float  # minutes
    intensity: float
    n_isotopes: int

    @property
    def neutral_mass(self) -> float:
        """Neutral monoisotopic mass assuming protonation (H+ only)."""
        return (self.mz - ch.PROTON_MASS) * self.charge


def to_msexperiment(scans: Sequence[tuple[float, np.ndarray, np.ndarray]]):
    """Build an MS1 ``MSExperiment`` from (time in minutes, m/z array, intensity array) tuples."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    for t, mz, it in scans:
        s = oms.MSSpectrum()
        s.setMSLevel(1)
        s.setRT(float(t) * 60.0)
        order = np.argsort(mz)
        s.set_peaks((np.asarray(mz, float)[order], np.asarray(it, float)[order]))
        exp.addSpectrum(s)
    exp.updateRanges()
    return exp


def scans_from_centroid(centroid_file, keep: np.ndarray | None, t_lo: float | None = None,
                        t_hi: float | None = None):
    """Yield (time, m/z, intensity) for kept scans inside the time window."""
    times = np.asarray(centroid_file.xlabels, float)
    for i, t in enumerate(times):
        if keep is not None and not keep[i]:
            continue
        if (t_lo is not None and t < t_lo) or (t_hi is not None and t > t_hi):
            continue
        mz, it = centroid_file.scan(i)
        yield float(t), np.asarray(mz, float), np.asarray(it, float)


def find_features(exp, ppm: float = 10.0, noise: float = 100.0, chrom_fwhm_s: float = 10.0,
                  charges: tuple[int, int] = (1, 4), min_trace_length_s: float = 5.0,
                  min_sample_rate: float = 0.5,
                  remove_single_traces: bool = True) -> list[MSFeature]:
    """Run MassTraceDetection, ElutionPeakDetection and FeatureFindingMetabo.

    ``isotope_filtering_model`` is set to ``none`` because the default model is trained on
    small metabolites; 13C spacing scoring (``mz_scoring_13C``) is switched on.

    ``noise`` must stay below the background level: OpenMS rejects traces by sampling rate
    (``min_sample_rate``), and when every background peak is below the threshold no trace
    passes (found on synthetic data with pyOpenMS 3.6.0).
    """
    import pyopenms as oms

    mtd = oms.MassTraceDetection()
    p = mtd.getDefaults()
    p.setValue("mass_error_ppm", float(ppm))
    p.setValue("noise_threshold_int", float(noise))
    p.setValue("min_trace_length", float(min_trace_length_s))
    p.setValue("min_sample_rate", float(min_sample_rate))
    mtd.setParameters(p)
    try:  # pyOpenMS >= 3.5 returns the traces
        traces = list(mtd.run(exp, 0))
    except TypeError:  # older pyOpenMS fills an output list
        traces = []
        mtd.run(exp, traces, 0)

    epd = oms.ElutionPeakDetection()
    p = epd.getDefaults()
    p.setValue("chrom_fwhm", float(chrom_fwhm_s))
    p.setValue("width_filtering", "fixed")
    epd.setParameters(p)
    try:
        split = list(epd.detectPeaks(traces))
    except TypeError:
        split = []
        epd.detectPeaks(traces, split)

    ffm = oms.FeatureFindingMetabo()
    p = ffm.getDefaults()
    p.setValue("charge_lower_bound", int(charges[0]))
    p.setValue("charge_upper_bound", int(charges[1]))
    p.setValue("chrom_fwhm", float(chrom_fwhm_s))
    p.setValue("isotope_filtering_model", "none")
    p.setValue("mz_scoring_13C", "true")
    p.setValue("remove_single_traces", "true" if remove_single_traces else "false")
    p.setValue("report_convex_hulls", "true")
    ffm.setParameters(p)
    fmap = oms.FeatureMap()
    try:
        ffm.run(split, fmap)
    except TypeError:
        ffm.run(split, fmap, [])

    out = []
    for i, f in enumerate(fmap):
        n_iso = len(f.getConvexHulls())
        out.append(MSFeature(i, float(f.getMZ()), int(f.getCharge()), float(f.getRT()) / 60.0,
                             float(f.getIntensity()), n_iso))
    return out


def write_mzml(exp, path: str) -> None:
    """Write an ``MSExperiment`` to mzML."""
    import pyopenms as oms

    oms.MzMLFile().store(str(path), exp)
