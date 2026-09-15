"""v0.2: v0.1 plus HQ-body relative amplitude gate for long decaying reads."""
import numpy as np

from .gates import amplitude_gates, noise_floor
from .v0_1 import (
    PERIOD_THRESH,
    SPIKE_THRESH,
    VALLEY_THRESH,
    amplitude_spike_z,
    classify as classify_v01,
    load_ab1,
    periodicity,
    valley_ratio,
)

ROLL_WIN = 21
ALPHA = 0.35
BETA = 0.18
PERIOD_HALF_WIN_WIDE = 60
PERIOD_HALF_WIN_LOCAL = 8
VALLEY_CLEAN_AFTER_BLOB = 0.42
VALLEY_BLOB_BEFORE = 0.45
RECOVERY_MIN_RUN = 5


def rolling_median(x, win=ROLL_WIN):
    x = np.asarray(x, dtype=float)
    n = len(x)
    half = win // 2
    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(x[max(0, i - half):min(n, i + half + 1)])
    return out


def longest_clean_island(v01, n):
    best = (1, 0, 0)
    i = 1
    while i <= n:
        if not v01[i]["bad"]:
            j = i
            while j <= n and not v01[j]["bad"]:
                j += 1
            if j - i > best[2]:
                best = (i, j - 1, j - i)
            i = j
        else:
            i += 1
    return best


def _clear_period_smear_after_blob(out, seq, ploc, channels, top_h, stop_rfu):
    """Wide FFT periodicity can flag clean peaks still in the smear of a prior blob (A2 294+)."""
    n = len(seq)
    per_w = periodicity(ploc, channels, PERIOD_HALF_WIN_WIDE)
    per_l = periodicity(ploc, channels, PERIOD_HALF_WIN_LOCAL)
    valley = valley_ratio(seq, ploc, channels)
    level = rolling_median(top_h)
    i = 1
    while i < n:
        if valley[i - 1] <= VALLEY_BLOB_BEFORE:
            i += 1
            continue
        j = i
        while j < n and valley[j] <= VALLEY_CLEAN_AFTER_BLOB and float(level[j]) >= stop_rfu:
            j += 1
        if j - i >= RECOVERY_MIN_RUN:
            for k in range(i, j):
                p = out[k + 1]
                if (
                    p["bad"]
                    and p["reason"] == "relative"
                    and valley[k] <= VALLEY_THRESH
                    and per_w[k] < PERIOD_THRESH
                    and per_l[k] >= PERIOD_THRESH
                ):
                    p["bad"] = False
                    p["reason"] = "ok"
        i = j if j > i else i + 1


def classify(ab1_path, alpha=ALPHA, beta=BETA):
    v01, raw = classify_v01(ab1_path)
    seq, ploc, channels, top_h = raw["seq"], raw["ploc"], raw["channels"], raw["top_h"]
    n = len(seq)
    a, b, w = longest_clean_island(v01, n)
    hq = float(np.median(top_h[a - 1:b])) if w else float(np.median(top_h))
    alpha_rfu, stop_rfu, floor = amplitude_gates(top_h, hq, alpha, beta)
    level = rolling_median(top_h)
    spike_z = amplitude_spike_z(top_h, signed=True)
    per = periodicity(ploc, channels)
    valley = valley_ratio(seq, ploc, channels)

    out = {}
    for i in range(n):
        rel_bad = (per[i] < PERIOD_THRESH) or (valley[i] > VALLEY_THRESH)
        spike_bad = spike_z[i] > SPIKE_THRESH
        lv = float(level[i])
        frac = lv / hq if hq else 0.0
        if lv < stop_rfu:
            bad, reason = True, "below_stop_level"
        elif lv < alpha_rfu:
            bad = bool(spike_bad)
            reason = "spike" if bad else "low_amp_keep"
        else:
            bad = bool(rel_bad or spike_bad)
            reason = "relative" if rel_bad else ("spike" if spike_bad else "ok")
        out[i + 1] = dict(
            bad=bool(bad),
            reason=reason,
            periodicity=float(per[i]),
            spike_z=float(spike_z[i]),
            valley_ratio=float(valley[i]),
            top_h=float(top_h[i]),
            local_level=lv,
            hq_body=hq,
            level_frac=float(frac),
            noise_floor=floor,
            alpha_rfu=alpha_rfu,
            stop_rfu=stop_rfu,
            v0_1_bad=bool(v01[i + 1]["bad"]),
            hq_island=(a, b),
            alpha=alpha,
            beta=beta,
        )
    _clear_period_smear_after_blob(out, seq, ploc, channels, top_h, stop_rfu)
    return out, raw
