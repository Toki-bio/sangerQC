"""Secondary peaks: is the second channel a real co-located base, or a neighbour's tail?

For every called position the slot is the trace between the midpoints to the
neighbouring calls. For each non-primary channel we take its highest interior
local maximum inside the slot (its own apex) and measure how far that apex sits
from the called peak, in units of mean peak spacing.

Classes (per position, for the strongest secondary channel):
  none         secondary / primary < RATIO_MIN
  spill        no interior apex, or apex offset > COLOCATE_MAX: the tail of a
               neighbouring peak (in practice almost always the neighbour's base)
  colocated    own apex within COLOCATE_MAX of the called peak: a real second
               signal at this position
  near_sat     co-located, but within SAT_WINDOW positions of an oversaturated
               peak (top_h >= SAT_X_HQ * HQ): spectral pull-up shows up exactly
               here, at the same scan in the other channels

Calibrated on two 18S reads (Thermocyclops, 2026-10-06): ~80 % of positions with
secondary/primary >= 0.33 were spill; a random wrong base passed the co-location
test at 0.9 % / 1.3 % of ordinary positions. Not yet a validated rule.
"""
import numpy as np

RATIO_MIN = 0.20      # below this nothing is recorded as a secondary
RATIO_CALL = 0.33     # sangerseqR makeBaseCalls default; IUPAC from a single read
COLOCATE_MAX = 0.35   # |apex offset| in units of mean peak spacing
SAT_X_HQ = 3.0
SAT_WINDOW = 10

IUPAC = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S",
    frozenset("AT"): "W", frozenset("GT"): "K", frozenset("AC"): "M",
}


def _interior_apex(x):
    """Index of the highest interior local maximum, or None."""
    if len(x) < 3:
        return None
    d = np.diff(x)
    idx = []
    k = 1
    while k < len(x) - 1:
        if d[k - 1] > 0:
            j = k
            while j < len(d) and d[j] == 0:   # plateau: a peak only if it then falls
                j += 1
            if j < len(d) and d[j] < 0:
                idx.append((k + j) // 2)
            k = j + 1
        else:
            k += 1
    if not idx:
        return None
    return max(idx, key=lambda k: x[k])


def slot(ploc, i, n_scans):
    n = len(ploc)
    left = (ploc[i - 1] + ploc[i]) // 2 if i > 0 else max(0, ploc[i] - 6)
    right = (ploc[i] + ploc[i + 1]) // 2 if i < n - 1 else min(n_scans - 1, ploc[i] + 6)
    if i > 0 and i < n - 1:
        spacing = (ploc[i + 1] - ploc[i - 1]) / 2.0
    elif n > 1:
        spacing = float(np.mean(np.diff(ploc)))
    else:
        spacing = 12.0
    return int(left), int(right), max(spacing, 1.0)


def channel_probe(raw, i):
    """Per-channel numbers for position i (0-based)."""
    ch, ploc = raw["channels"], raw["ploc"]
    n_scans = len(ch["A"])
    L, R, sp = slot(ploc, i, n_scans)
    h = {b: float(ch[b][L:R + 1].max()) for b in "ACGT"}
    apex = {}
    for b in "ACGT":
        k = _interior_apex(ch[b][L:R + 1])
        if k is not None:
            apex[b] = (float(ch[b][L + k]), (L + k - ploc[i]) / sp)
    # Primary = the tallest channel with a real apex on the call. The tallest value anywhere in the
    # slot can be the tail of a taller neighbouring peak (e.g. a small G after a big C), which
    # made a true single base look like a double (found on the allergology plates, 2026-10-07).
    on_call = {b: v[0] for b, v in apex.items() if abs(v[1]) <= COLOCATE_MAX}
    primary = max(on_call, key=on_call.get) if on_call else max(h, key=h.get)
    top = on_call[primary] if on_call else h[primary]
    out = {"primary": primary, "top": top, "channels": {}}
    for b in "ACGT":
        if b not in apex:
            out["channels"][b] = dict(height=h[b], ratio=h[b] / top if top else 0.0, apex=False, offset=None)
        else:
            out["channels"][b] = dict(height=apex[b][0], ratio=apex[b][0] / top if top else 0.0,
                                      apex=True, offset=float(apex[b][1]))
    return out


def is_colocated(c, ratio_min=RATIO_MIN, colocate_max=COLOCATE_MAX):
    return c["apex"] and c["offset"] is not None and abs(c["offset"]) <= colocate_max and c["ratio"] >= ratio_min


def saturated_zone(top_h, hq, x_hq=SAT_X_HQ, window=SAT_WINDOW):
    top_h = np.asarray(top_h, dtype=float)
    hot = np.where(top_h >= x_hq * hq)[0]
    zone = np.zeros(len(top_h), dtype=bool)
    for k in hot:
        zone[max(0, k - window):k + window + 1] = True
    return zone


def secondary_peaks(raw, hq, ratio_min=RATIO_MIN, ratio_call=RATIO_CALL, colocate_max=COLOCATE_MAX):
    """List (0-based) of per-position secondary-peak records."""
    n = len(raw["seq"])
    zone = saturated_zone(raw["top_h"], hq)
    recs = []
    for i in range(n):
        p = channel_probe(raw, i)
        prim = p["primary"]
        others = [b for b in "ACGT" if b != prim]
        sec = max(others, key=lambda b: p["channels"][b]["height"])
        c = p["channels"][sec]
        ratio_slot = p["channels"][sec]["height"] / p["top"] if p["top"] else 0.0
        # strongest co-located candidate, if any (may differ from the tallest secondary)
        coloc = [b for b in others if is_colocated(p["channels"][b], ratio_min, colocate_max)]
        if coloc:
            sec = max(coloc, key=lambda b: p["channels"][b]["ratio"])
            c = p["channels"][sec]
            cls = "near_sat" if zone[i] else "colocated"
        elif ratio_slot < ratio_min:
            cls = "none"
        else:
            cls = "spill"
        iupac = None
        if cls == "colocated" and c["ratio"] >= ratio_call:
            iupac = IUPAC[frozenset((prim, sec))]
        recs.append(dict(
            pos=i + 1, primary=prim, secondary=sec, h_primary=p["top"], h_secondary=c["height"],
            ratio=round(float(c["ratio"]), 3), apex=c["apex"],
            offset=None if c["offset"] is None else round(c["offset"], 3),
            cls=cls, iupac=iupac, probe=p,
        ))
    return recs


def background_rate(raw, recs, positions, ratio_min=RATIO_MIN, colocate_max=COLOCATE_MAX):
    """Control: how often does a non-called base pass the co-location test at ordinary positions?"""
    tot = hit = 0
    for i in positions:
        r = recs[i]
        for b in "ACGT":
            if b == r["primary"]:
                continue
            tot += 1
            hit += bool(is_colocated(r["probe"]["channels"][b], ratio_min, colocate_max))
    return hit / tot if tot else float("nan")
