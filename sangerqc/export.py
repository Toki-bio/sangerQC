"""Turn v0.3 output into a sequence plus a per-position decision table.

Letter rules (inside the kept span):
  co-located secondary, ratio >= 0.33 (or reproduced across reads)  -> IUPAC
  basecaller IUPAC but the secondary is a neighbour's tail / absent -> primary base
  co-located secondary inside a saturated zone (pull-up)            -> caller letter, lowercase
  position flagged bad, or leading low-amplitude run (low_amp_5prime) -> lowercase
Ends: primer sites (if primers given) and leading/trailing stretches with no run
of MIN_CLEAN_RUN consecutive unflagged positions are trimmed.
Interior bad stretches are kept lowercase, never N-filled here: N needs reference
(BLAST) evidence of disagreement, which this package does not fetch.
"""
from .secondary import IUPAC

MIN_CLEAN_RUN = 10
IUPAC_LETTERS = set("RYSWKMBDHVN")

FIELDS = ["pos", "scan", "caller", "phred", "out", "decision", "bad", "reason", "reasons", "local_level",
          "primary", "secondary", "ratio", "apex", "offset", "sec_class", "lbf", "p_mixed", "gt", "gq", "pl",
          "concordance",
          "h_A", "h_C", "h_G", "h_T", "off_A", "off_C", "off_G", "off_T"]


def kept_span(pred, n, min_clean_run=MIN_CLEAN_RUN):
    ok = [not pred[i]["bad"] and not pred[i].get("in_primer_or_beyond", False) for i in range(1, n + 1)]
    start = end = None
    run = 0
    for i, v in enumerate(ok):
        run = run + 1 if v else 0
        if run == min_clean_run:
            start = i - min_clean_run + 2
            break
    run = 0
    for i in range(n - 1, -1, -1):
        run = run + 1 if ok[i] else 0
        if run == min_clean_run:
            end = i + min_clean_run
            break
    if start is None or end is None or end < start:
        return None
    return start, end


def reference_trim(span, ref_status, window=10, max_bad=0.15):
    """Shrink span from both ends while the end window disagrees with the reference too often.

    ref_status: {read_pos1: True (match/compatible) | False (mismatch, or flanks a gap)};
    positions outside the local alignment count as disagreement. 0.15 = the Artemia
    project's bar for calling a run really bad against BLAST.
    """
    a, b = span

    def bad_frac(lo, hi):
        vals = [ref_status.get(i, False) for i in range(lo, hi + 1)]
        return sum(1 for v in vals if not v) / len(vals)

    while b - a + 1 > window and bad_frac(a, a + window - 1) > max_bad:
        a += 1
    while b - a + 1 > window and bad_frac(b - window + 1, b) > max_bad:
        b -= 1
    while a < b and not ref_status.get(a, False):   # a lone disagreeing end base
        a += 1
    while b > a and not ref_status.get(b, False):
        b -= 1
    return a, b


def to_sequence(pred, raw, reproduced=None, min_clean_run=MIN_CLEAN_RUN, ref_status=None):
    """reproduced: {pos(1-based): frozenset pair} from concordance.reproduced_sites().
    ref_status: per-position agreement with a reference (concordance.reference_status), used to trim ends."""
    seq = raw["seq"]
    n = len(seq)
    reproduced = reproduced or {}
    span = kept_span(pred, n, min_clean_run)
    if span is not None and ref_status is not None:
        span = reference_trim(span, ref_status)
    rows, out = [], []
    if span is None:
        return "", rows, None
    for i in range(span[0], span[1] + 1):
        p = pred[i]
        s = p["secondary"]
        caller = seq[i - 1].upper()
        letter, decision = caller, "caller"
        if i in reproduced:
            letter, decision = IUPAC[reproduced[i]], "IUPAC: co-located, reproduced across reads"
        elif s["cls"] == "colocated" and s["iupac"]:
            letter, decision = s["iupac"], "IUPAC: co-located secondary >= 0.33"
        elif s["cls"] == "near_sat":
            # inside a blob the tallest channel is the saturated one bleeding over; keep the caller's letter
            letter, decision = caller.lower(), "caller: co-located but saturated zone (pull-up)"
        elif caller in IUPAC_LETTERS and s["cls"] in ("spill", "none"):
            letter, decision = s["primary"], f"primary: caller IUPAC but secondary is {s['cls']}"
        elif caller in IUPAC_LETTERS and s["cls"] == "colocated":
            letter, decision = s["primary"], "primary: co-located secondary < 0.33"
        if p["bad"] or p["reason"] == "low_amp_5prime":
            letter = letter.lower()
        out.append(letter)
        if decision != "caller" or s["cls"] != "none" or p["bad"] or p["reason"] == "low_amp_5prime":
            chans = raw["secondary_records"][i - 1]["probe"]["channels"]
            row = dict(
                pos=i, scan=int(raw["ploc"][i - 1]), caller=caller, phred=raw["phred"][i - 1], out=letter,
                decision=decision, bad=p["bad"], reason=p["reason"], reasons="+".join(p.get("reasons", [])),
                local_level=round(p["local_level"]), primary=s["primary"], secondary=s["secondary"],
                ratio=s["ratio"], apex=s["apex"], offset=s["offset"], sec_class=s["cls"],
                lbf=p.get("lbf", ""), p_mixed=p.get("p_mixed", ""), gt=p.get("gt", ""), gq=p.get("gq", ""),
                pl=p.get("pl", ""), concordance="reproduced" if i in reproduced else "",
            )
            for b in "ACGT":
                row[f"h_{b}"] = round(chans[b]["height"])
                row[f"off_{b}"] = "" if chans[b]["offset"] is None else round(chans[b]["offset"], 2)
            rows.append(row)
    return "".join(out), rows, span
