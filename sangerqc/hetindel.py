"""Heterozygous indel: split a mixed trace into its two allele sequences (prototype).

Model. Before the indel (call `pos`) both alleles agree. From `pos` on, one
allele lags the other by d calls:
    S_i = {P_i, Q_i},   Q_i = P_{i-d}   for i >= pos + d
and Q_pos .. Q_{pos+d-1} are the d bases present in only one allele.
S_i is the set of bases with a peak *on* call i: the primary channel plus every
channel whose own apex is co-located (secondary.py, ratio >= RATIO_MIN). Using
co-location rather than the caller's IUPAC letters matters: a neighbour's tail
(spill) puts a false second base into S_i and breaks the chain logic.

Positions split into d independent chains (i = pos + r, pos + r + d, ...).
Choosing P at the head of a chain fixes the whole chain, so each chain has two
phasings; we keep the one with fewer contradictions (Q_i not in S_i). The
indel size d and onset pos are chosen by the lowest contradiction rate.
This is the idea of Dmitriev & Rakitov 2008 (PLoS Comput Biol;
reference-free decoding from the peak string) with chain-wise rather than
full dynamic-programming search. Which allele is longer is not decided here;
a reference (or the wild-type read) names it afterwards (`name_variant`).
"""
import numpy as np

from .secondary import is_colocated

RATIO_MIN = 0.20
MAX_D = 40


def peak_sets(raw, ratio_min=RATIO_MIN):
    sets = []
    for r in raw["secondary_records"]:
        p = r["probe"]
        s = {p["primary"]}
        for b in "ACGT":
            if b != p["primary"] and is_colocated(p["channels"][b], ratio_min):
                s.add(b)
        sets.append(s)
    return sets


def onset(sets, start=0, end=None, before=30, min_after=0.4, max_before=0.25):
    """Start (0-based) of the persistent double-peak run: the step that maximises
    (double fraction from i to end) - (double fraction in the `before` calls before i).
    None if no step with >= min_after after and <= max_before before."""
    end = len(sets) if end is None else end
    dbl = np.array([len(s) >= 2 for s in sets[:end]], dtype=float)
    best, best_gain = None, 0.0
    for i in range(start + 5, end - 10):
        after = dbl[i:end].mean()
        b0 = max(start, i - before)
        prior = dbl[b0:i].mean() if i > b0 else 0.0
        if after >= min_after and prior <= max_before and after - prior > best_gain and dbl[i]:
            best, best_gain = i, after - prior
    return best


def decode(sets, pos, d, end=None):
    """Return (P, Q, contradictions, n_checked) for onset pos (0-based) and lag d >= 1."""
    end = len(sets) if end is None else end
    P = [None] * end
    Q = [None] * end
    for i in range(pos):
        P[i] = Q[i] = min(sets[i]) if len(sets[i]) == 1 else sorted(sets[i])[0]
    bad_total = 0
    checked = 0
    for r in range(d):
        head = pos + r
        if head >= end:
            continue
        best = None
        for first in sorted(sets[head]):
            p_chain, q_chain, bad = {}, {}, 0
            i = head
            p_chain[i] = first
            q_chain[i] = next(iter(sets[i] - {first}), first)
            i += d
            while i < end:
                q = p_chain[i - d]
                if q not in sets[i]:
                    bad += 1
                q_chain[i] = q
                rest = sets[i] - {q}
                p_chain[i] = next(iter(rest)) if len(rest) == 1 else (q if not rest else sorted(rest)[0])
                i += d
            if best is None or bad < best[2]:
                best = (p_chain, q_chain, bad)
        for i, v in best[0].items():
            P[i] = v
        for i, v in best[1].items():
            Q[i] = v
        bad_total += best[2]
        checked += len(best[0]) - 1
    return "".join(P), "".join(Q), bad_total, checked


def analysis_span(pred, raw, min_clean_run=10):
    """Calls usable for decoding: from the first run of min_clean_run unflagged calls to the last (0-based, end exclusive)."""
    from .export import kept_span
    span = kept_span(pred, len(raw["seq"]), min_clean_run)
    return (0, len(raw["seq"])) if span is None else (span[0] - 1, span[1])


def solve(sets, max_d=MAX_D, pos=None, start=0, end=None, search=6):
    """Best (pos, d) by contradiction rate. Returns dict or None if no double-peak run."""
    end = len(sets) if end is None else end
    p0 = onset(sets, start=start, end=end) if pos is None else pos
    if p0 is None:
        return None
    results = []
    for p in range(max(start, p0 - search), min(end - 1, p0 + search) + 1):
        for d in range(1, max_d + 1):
            if p + 2 * d >= end:
                break
            P, Q, bad, n = decode(sets, p, d, end)
            if n < 10:
                continue
            results.append((bad / n, d, p, P, Q, bad, n))
    if not results:
        return None
    results.sort(key=lambda x: (x[0], x[1], -x[2]))
    rate, d, p, P, Q, bad, n = results[0]
    second = next((x for x in results if x[1] != d), None)
    return dict(pos=p, d=d, P=P, Q=Q, contradiction_rate=rate, contradictions=bad, checked=n,
                runner_up=None if second is None else dict(d=second[1], pos=second[2], rate=second[0]))


def name_variant(P, Q, reference):
    """Which allele matches the reference better, and the other one's indel against it."""
    from Bio.Align import PairwiseAligner
    al = PairwiseAligner()
    al.mode = "local"
    al.match_score, al.mismatch_score, al.open_gap_score, al.extend_gap_score = 2, -3, -6, -1
    sp, sq = al.score(P, reference), al.score(Q, reference)
    ref_like, other = (P, Q) if sp >= sq else (Q, P)
    aln = al.align(other, reference)[0]
    gaps = []
    (qb, rb) = aln.aligned
    for j in range(1, len(qb)):
        dq = qb[j][0] - qb[j - 1][1]
        dr = rb[j][0] - rb[j - 1][1]
        if dq > 0 and dr == 0:
            gaps.append(("insertion", int(rb[j][0]) + 1, other[qb[j - 1][1]:qb[j][0]]))
        elif dr > 0 and dq == 0:
            gaps.append(("deletion", int(rb[j - 1][1]) + 1, reference[rb[j - 1][1]:rb[j][0]]))
    return dict(ref_like_allele="P" if sp >= sq else "Q", score_P=sp, score_Q=sq, indels_in_other=gaps)
