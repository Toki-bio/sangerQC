"""v0.1 Sanger chromatogram classifier: periodicity + spike + valley-depth."""
import numpy as np
from Bio import SeqIO

PERIOD_HALF_WIN = 60
PERIOD_THRESH = 0.3
SPIKE_THRESH = 1.5
VALLEY_THRESH = 0.5


def load_ab1(ab1_path):
    rec = SeqIO.read(ab1_path, "abi")
    ann = rec.annotations["abif_raw"]
    seq = str(rec.seq)
    base_order = ann["FWO_1"].decode() if isinstance(ann["FWO_1"], bytes) else ann["FWO_1"]
    channels = {
        b: np.array(ann[ch], dtype=float)
        for ch, b in zip(["DATA9", "DATA10", "DATA11", "DATA12"], base_order)
    }
    ploc = np.array(ann.get("PLOC2") or ann.get("PLOC1"))
    top_h = np.array(
        [channels[max("ACGT", key=lambda b: channels[b][p])][p] for p in ploc]
    )
    phred = list(rec.letter_annotations.get("phred_quality", [None] * len(seq)))
    return {
        "seq": seq,
        "ploc": ploc,
        "channels": channels,
        "top_h": top_h,
        "base_order": base_order,
        "phred": phred,
    }


def periodicity(ploc, channels, half_win=PERIOD_HALF_WIN):
    n_scans = len(channels["A"])
    envelope = np.maximum.reduce([channels[b] for b in "ACGT"])
    avg_spacing = np.mean(np.diff(ploc))
    expected_freq = 1.0 / avg_spacing

    def strength(center_scan):
        lo = max(0, center_scan - half_win)
        hi = min(n_scans, center_scan + half_win)
        seg = envelope[lo:hi] - envelope[lo:hi].mean()
        if len(seg) < 10 or seg.std() == 0:
            return 0.0
        fft = np.fft.rfft(seg)
        freqs = np.fft.rfftfreq(len(seg))
        power = np.abs(fft) ** 2
        total = power.sum()
        if total == 0:
            return 0.0
        tol = expected_freq * 0.4
        mask = (freqs >= expected_freq - tol) & (freqs <= expected_freq + tol)
        return float(power[mask].sum() / total)

    return np.array([strength(p) for p in ploc])


def amplitude_spike_z(top_h, signed=False):
    log_h = np.log1p(top_h)
    sd = log_h.std()
    z = (log_h - log_h.mean()) / sd if sd > 0 else np.zeros_like(log_h)
    return z if signed else np.abs(z)


def valley_ratio(seq, ploc, channels):
    n = len(ploc)
    envelope = np.maximum.reduce([channels[b] for b in "ACGT"])

    def gap_ratio(i):
        if seq[i].upper() == seq[i + 1].upper() and seq[i].upper() in "ACGT":
            return 0.0
        lo, hi = ploc[i], ploc[i + 1]
        if hi <= lo:
            return 1.0
        valley_min = envelope[lo:hi + 1].min()
        peak_avg = (envelope[lo] + envelope[hi]) / 2.0
        return float(valley_min / peak_avg) if peak_avg > 0 else 1.0

    gap = [gap_ratio(i) for i in range(n - 1)]
    out = np.zeros(n)
    for i in range(n):
        left = gap[i - 1] if i > 0 else gap[0]
        right = gap[i] if i < n - 1 else gap[-1]
        out[i] = (left + right) / 2.0
    return out


def classify(ab1_path):
    d = load_ab1(ab1_path)
    per = periodicity(d["ploc"], d["channels"])
    spike = amplitude_spike_z(d["top_h"])
    valley = valley_ratio(d["seq"], d["ploc"], d["channels"])
    bad = (per < PERIOD_THRESH) | (spike > SPIKE_THRESH) | (valley > VALLEY_THRESH)
    seq = d["seq"]
    return {
        i + 1: dict(
            bad=bool(bad[i]),
            periodicity=float(per[i]),
            spike_z=float(spike[i]),
            valley_ratio=float(valley[i]),
        )
        for i in range(len(seq))
    }, d
