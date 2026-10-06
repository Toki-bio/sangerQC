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


def peak_sets(raw, ratio_min=RATIO_MIN, max_bases=2):
    """Bases with a peak on each call: primary + co-located secondaries >= ratio_min.
    max_bases=2: a two-allele mixture has at most two bases per call, so only the strongest
    co-located secondary is kept (extra ones are noise and make every shift fit a little)."""
    sets = []
    for r in raw["secondary_records"]:
        p = r["probe"]
        cand = [(p["channels"][b]["ratio"], b) for b in "ACGT"
                if b != p["primary"] and is_colocated(p["channels"][b], ratio_min)]
        cand.sort(reverse=True)
        s = {p["primary"]} | {b for _, b in cand[:max_bases - 1]}
        sets.append(s)
    return sets


def adaptive_floor(m, start, onset_pos, q=0.9, floor=RATIO_MIN, cap=0.45):
    """Per-read floor for counting a secondary: the q-quantile of the mixing signal in the
    clean stretch before the onset (that read's own spurious co-located seconds), never
    below `floor`, never above `cap`."""
    seg = m[start:onset_pos]
    if len(seg) < 10:
        return floor
    return float(min(cap, max(floor, np.quantile(seg, q) + 0.02)))


def mixing_signal(raw, ratio_min=0.0):
    """Per call: ratio of the strongest co-located secondary (0 if none). Continuous, so a
    weak minor allele still raises it (the idea of Tracy's best-minus-second signal, but
    with neighbour tails excluded)."""
    m = []
    for r in raw["secondary_records"]:
        p = r["probe"]
        best = 0.0
        for b in "ACGT":
            c = p["channels"][b]
            if b != p["primary"] and is_colocated(c, ratio_min):
                best = max(best, c["ratio"])
        m.append(best)
    return np.array(m)


def onset_continuous(m, start=0, end=None, win=25, min_step=0.08, min_persist=0.35, sets=None):
    """Breakpoint = call maximising mean(m, right window) - mean(m, left window) (Tracy's
    findBreakpoint in decompose.h uses 25-call windows). Accepted only if the step is
    >= min_step and double peaks persist (>= min_persist of calls from there to the end): a
    block of double peaks that ends again (e.g. intragenomic variants) is not an indel shift."""
    end = len(m) if end is None else end
    best, best_step = None, 0.0
    for i in range(start + win, end - win):
        step = m[i:i + win].mean() - m[i - win:i].mean()
        if step > best_step:
            best, best_step = i, step
    if best is None or best_step < min_step:
        return None, best_step
    if sets is not None:
        persist = np.mean([len(s) >= 2 for s in sets[best:end]])
        if persist < min_persist:
            return None, best_step
    return best, best_step


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
    start, end = (0, len(raw["seq"])) if span is None else (span[0] - 1, span[1])
    # a leading low-amplitude run (low_amp_5prime) is not a clean baseline either
    lead = raw.get("diagnostics", {}).get("lead_low_amp_end", 0) or 0
    return max(start, lead), end


def size_call(rates, tol=0.02, z_min=3.0, frac_max=0.5):
    """Pick the indel size from {d: best contradiction rate}.

    A real shift drops far below what other sizes reach by chance (Tracy: median - MAD
    cut-off; Indelligent: large shifts fit by chance, so prefer the smallest size). Returns
    (d, z, accepted): smallest d within tol of the minimum; accepted if its rate is
    <= frac_max x median and >= z_min robust SDs below the median."""
    ds = sorted(rates)
    vals = np.array([rates[d] for d in ds])
    med = float(np.median(vals))
    mad = float(np.median(np.abs(vals - med))) * 1.4826
    lo = float(vals.min())
    d = next(d for d in ds if rates[d] <= lo + tol)
    z = (med - rates[d]) / max(mad, 1e-3)
    return d, z, bool(rates[d] <= frac_max * med and z >= z_min)


def allele_fraction(raw, P, Q, start, end):
    """Mean share of the Q allele's peak height where the alleles differ (peak heights are
    only semi-quantitative)."""
    fr = []
    for i in range(start, min(end, len(P), len(Q))):
        if P[i] != Q[i] and P[i] in "ACGT" and Q[i] in "ACGT":
            ch = raw["secondary_records"][i]["probe"]["channels"]
            hp, hq = ch[P[i]]["height"], ch[Q[i]]["height"]
            if hp + hq > 0:
                fr.append(hq / (hp + hq))
    return float(np.mean(fr)) if fr else None


def solve(sets, max_d=MAX_D, pos=None, start=0, end=None, search=6, raw=None, _floor_tried=False):
    """Best (pos, d) by contradiction rate, with a significance call. None if no indel-like step.

    With `raw`, the onset comes from the continuous mixing signal plus a persistence check;
    without it, from the binary double-peak step (older behaviour)."""
    end = len(sets) if end is None else end
    step = None
    if pos is not None:
        p0 = pos
    elif raw is not None:
        m = mixing_signal(raw)
        p0, step = onset_continuous(m, start=start, end=end, sets=sets)
        if p0 is not None and not _floor_tried:
            # a per-read floor removes noise seconds but can also remove a weak minor allele
            # (20:80 mix -> ratio ~0.25): run both and keep the more significant answer
            fl = adaptive_floor(m, start, p0)
            if fl > RATIO_MIN + 1e-9:
                a = solve(sets, max_d, p0, start, end, search, raw, _floor_tried=True)
                b = solve(peak_sets(raw, ratio_min=fl), max_d, p0, start, end, search, raw, _floor_tried=True)
                for r_ in (a, b):
                    if r_ is not None:
                        r_["step"] = step
                if b is not None:
                    b["floor"] = fl
                cands = [r_ for r_ in (a, b) if r_ is not None]
                return max(cands, key=lambda r_: (r_["indel"], r_["z"])) if cands else None
    else:
        p0 = onset(sets, start=start, end=end)
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
    best_by_d = {}
    for x in results:
        if x[1] not in best_by_d or (x[0], -x[2]) < (best_by_d[x[1]][0], -best_by_d[x[1]][2]):
            best_by_d[x[1]] = x
    d, z, accepted = size_call({k: v[0] for k, v in best_by_d.items()})
    rate, d, p, P, Q, bad, n = best_by_d[d]
    second = sorted((v for k, v in best_by_d.items() if k != d), key=lambda x: x[0])
    out = dict(pos=p, d=d, P=P, Q=Q, contradiction_rate=rate, contradictions=bad, checked=n,
               z=round(z, 2), indel=accepted, step=step,
               runner_up=None if not second else dict(d=second[0][1], pos=second[0][2], rate=second[0][0]))
    if raw is not None:
        out["allele_fraction_Q"] = allele_fraction(raw, P, Q, p, end)
    return out


def solve_reference(sets, reference, call_seq, start=0, end=None, max_d=MAX_D, search=6, raw=None, _floor_tried=False):
    """Reference-guided decoding (the idea of Poly Peak Parser / Tracy decompose): one allele is
    the reference; scan signed shifts s of the other allele against it.

    Calls before the onset are placed on the reference by local alignment of the clean prefix;
    from there call i sits on reference index r(i) = i + offset (the reference allele has no
    indel by definition). s > 0: the sample allele skips s reference bases (deletion);
    s < 0: it carries |s| extra bases (insertion), read from the peak sets.
    """
    from Bio.Align import PairwiseAligner
    end = len(sets) if end is None else end
    m = mixing_signal(raw) if raw is not None else np.array([float(len(s) >= 2) for s in sets])
    p0, step = onset_continuous(m, start=start, end=end, sets=sets)
    if p0 is None:
        return None
    if raw is not None and not _floor_tried:
        fl = adaptive_floor(m, start, p0)
        if fl > RATIO_MIN + 1e-9:
            a = solve_reference(sets, reference, call_seq, start, end, max_d, search, raw, _floor_tried=True)
            b = solve_reference(peak_sets(raw, ratio_min=fl), reference, call_seq, start, end, max_d, search, raw,
                                _floor_tried=True)
            if b is not None:
                b["floor"] = fl
            cands = [r_ for r_ in (a, b) if r_ is not None]
            return max(cands, key=lambda r_: (r_["indel"], r_["z"])) if cands else None
    al = PairwiseAligner()
    al.mode = "local"
    al.match_score, al.mismatch_score, al.open_gap_score, al.extend_gap_score = 2, -3, -6, -1
    prefix = call_seq[start:p0]
    aln = al.align(prefix, reference)[0]
    qa, ra = int(aln.aligned[0][-1][1]), int(aln.aligned[1][-1][1])   # end of last aligned block
    offset = ra - (start + qa)
    ref = reference.upper()
    results = []
    for p in range(max(start, p0 - search), min(end - 1, p0 + search) + 1):
        for s in [x for d in range(1, max_d + 1) for x in (d, -d)]:
            bad = n = 0
            first = p + (abs(s) if s < 0 else 0)
            for i in range(first, end):
                ri, qi = i + offset, i + offset + s
                if not (0 <= ri < len(ref) and 0 <= qi < len(ref)):
                    continue
                n += 1
                if ref[ri] not in sets[i] or ref[qi] not in sets[i]:
                    bad += 1
            # calls between the aligned prefix and p must still look like the reference
            for i in range(start + qa, p):
                ri = i + offset
                if 0 <= ri < len(ref):
                    n += 1
                    bad += ref[ri] not in sets[i]
            if n >= 10:
                results.append((bad / n, s, p, bad, n))
    if not results:
        return None
    best_by_s = {}
    for x in results:
        if x[1] not in best_by_s or x[0] < best_by_s[x[1]][0]:
            best_by_s[x[1]] = x
    rates = {}
    for s, x in best_by_s.items():
        rates[abs(s)] = min(rates.get(abs(s), 1.0), x[0])
    d, z, accepted = size_call(rates)
    rate, s, p, bad, n = min((best_by_s[x] for x in (d, -d) if x in best_by_s), key=lambda x: x[0])
    P, Q = [], []
    for i in range(end):
        ri = i + offset
        pr = ref[ri] if 0 <= ri < len(ref) else min(sets[i])
        P.append(pr)
        if i < p:
            Q.append(pr)
        elif s < 0 and i < p - s:
            rest = sets[i] - {pr}
            Q.append(next(iter(rest)) if rest else pr)       # inserted base
        else:
            qi = ri + s
            Q.append(ref[qi] if 0 <= qi < len(ref) else "N")
    P, Q = "".join(P), "".join(Q)
    out = dict(pos=p, shift=s, d=d, kind="deletion" if s > 0 else "insertion", P=P, Q=Q,
               contradiction_rate=rate, contradictions=bad, checked=n, z=round(z, 2), indel=accepted,
               step=step, ref_offset=offset)
    if raw is not None:
        out["allele_fraction_Q"] = allele_fraction(raw, P, Q, p, end)
    return out


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


def decode_read(pred, raw, reference=None, max_d=MAX_D):
    """Entry point: decode a v0.3-classified read inside its unflagged span only.

    Decoding a dying 3' tail invents shifts (Artemia H4: a "1-bp indel" at call 807, inside the
    BLAST-confirmed degraded zone; none within the unflagged span). With a reference the
    reference-guided mode is used (stronger in every benchmark condition); otherwise the
    reference-free mode. Returns None when there is no indel-like step."""
    start, end = analysis_span(pred, raw)
    sets = peak_sets(raw)
    if reference:
        return solve_reference(sets, reference, raw["seq"], start=start, end=end, max_d=max_d, raw=raw)
    return solve(sets, max_d=max_d, start=start, end=end, raw=raw)
