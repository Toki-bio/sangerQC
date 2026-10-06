"""v0.3: v0.2 flags plus secondary peaks, all-reasons reporting, 5' start-up rule, primer ends.

Bad/ok flags are identical to v0.2 (checked on the six Artemia reads and two
Thermocyclops reads). Everything is additive:
- a leading low-amplitude run (before the local level first reaches alpha) is
  relabelled "low_amp_5prime" (not bad): v0.2 called it low_amp_keep, a class
  designed for the slowly decaying 3' tail. It can be start-up junk or real but
  weak gene, so export writes it lowercase and only a reference check trims it.
- `reasons`: every test that fired, not only the first. The single `reason`
  label now prefers "spike" over "relative" when both fire (oversaturation is
  the more specific class).
- secondary-peak record per position (see secondary.py).
- file diagnostics: low_snr when HQ < 2 x noise floor (the 30 RFU floor then
  swallows the whole read and the amplitude gates stop meaning anything), and
  the read's own background under called peaks.
"""
import numpy as np

from .primers import primer_mask
from .secondary import secondary_peaks
from .v0_1 import PERIOD_THRESH, SPIKE_THRESH, VALLEY_THRESH
from .v0_2 import ALPHA, BETA, classify as classify_v02, rolling_median

LOW_SNR_X = 2.0


def file_diagnostics(raw, hq, floor, island):
    ch, ploc = raw["channels"], raw["ploc"]
    a, b = island
    second = [sorted(ch[x][p] for x in "ACGT")[-2] for p in ploc[a - 1:b]]
    background = float(np.median(second)) if second else float("nan")
    return dict(
        hq=hq,
        noise_floor=floor,
        snr=hq / floor if floor else float("inf"),
        low_snr=bool(hq < LOW_SNR_X * floor),
        background_under_peaks=background,
        snr_vs_background=hq / max(background, 1.0),
    )


def classify(ab1_path, alpha=ALPHA, beta=BETA, forward=None, reverse=None):
    pred, raw = classify_v02(ab1_path, alpha=alpha, beta=beta)
    n = len(raw["seq"])
    one = pred[1]
    hq, alpha_rfu, stop_rfu = one["hq_body"], one["alpha_rfu"], one["stop_rfu"]
    level = rolling_median(raw["top_h"])

    # leading low-amplitude run
    lead = next((i for i in range(n) if level[i] >= alpha_rfu), n)

    for i in range(1, n + 1):
        p = pred[i]
        reasons = []
        if p["local_level"] < stop_rfu:
            reasons.append("below_stop_level")
        elif p["local_level"] < alpha_rfu:
            reasons.append("low_amp")
        if p["periodicity"] < PERIOD_THRESH:
            reasons.append("periodicity")
        if p["valley_ratio"] > VALLEY_THRESH:
            reasons.append("valley")
        if p["spike_z"] > SPIKE_THRESH:
            reasons.append("spike")
        p["reasons"] = reasons
        p["v0_2_bad"], p["v0_2_reason"] = p["bad"], p["reason"]
        if p["bad"] and p["reason"] == "relative" and "spike" in reasons:
            p["reason"] = "spike"
        if p["reason"] == "low_amp_keep" and i - 1 < lead:
            # caution, not bad: Artemia A2 44-65 is real, BLAST-confirmed gene in exactly this
            # state; A3 49-144 and Thermocyclops F2 4-30 are not the locus. Only reference
            # evidence (export.reference_trim) decides; without it these stay lowercase.
            p["reason"] = "low_amp_5prime"

    sec = secondary_peaks(raw, hq)
    for r in sec:
        pred[r["pos"]]["secondary"] = {k: v for k, v in r.items() if k != "probe"}

    island = one["hq_island"]
    diag = file_diagnostics(raw, hq, one["noise_floor"], island)
    diag["lead_low_amp_end"] = lead
    if forward or reverse:
        s, e, notes = primer_mask(raw["seq"], forward, reverse)
        diag["primer_span"] = (s, e)
        diag["primer_notes"] = notes
        for i in range(1, n + 1):
            pred[i]["in_primer_or_beyond"] = not (s <= i <= e)
    raw["secondary_records"] = sec
    raw["diagnostics"] = diag
    return pred, raw
