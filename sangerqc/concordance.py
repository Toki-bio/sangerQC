"""Multi-read concordance for secondary peaks.

Reads of the same fragment (other specimens, other primers, re-reads) are put on
a shared reference coordinate. A site is *reproduced* when at least two reads
carry the same pair of bases at it, each with a co-located secondary apex at
ratio >= RATIO_MIN, outside bad and saturated positions. Which base dominates may
differ between reads.

Evidence strength (record it, see docs/SECONDARY_PEAKS.md):
  FR  forward + reverse of the same template: independent sequence context — strong
  DS  different specimens, same primer: independent templates, same sequence
      context — a context-driven artifact would reproduce too — weaker
"""
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from .secondary import is_colocated

RATIO_MIN = 0.20
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def _aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score, a.mismatch_score = 2, -3
    a.open_gap_score, a.extend_gap_score = -5, -2
    return a


def map_to_reference(seq, reference):
    """{ref_index0: (read_index0, strand)} using the better strand."""
    al = _aligner()
    best = None
    for strand, q in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        aln = al.align(q, reference)[0]
        if best is None or aln.score > best[0]:
            best = (aln.score, strand, aln)
    _, strand, aln = best
    n = len(seq)
    m = {}
    for (qa, qb), (ra, rb) in zip(*aln.aligned):
        for k in range(qb - qa):
            qi = qa + k
            m[ra + k] = (qi if strand == "+" else n - 1 - qi, strand)
    return m


IUPAC_SET = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
             "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def reference_status(seq, reference):
    """{read_pos1: bool}: base agrees (IUPAC-compatible) with the reference in the local alignment.
    A base flanking an alignment gap counts as disagreement."""
    al = _aligner()
    best = None
    for strand, q in (("+", seq), ("-", str(Seq(seq).reverse_complement()))):
        aln = al.align(q, reference)[0]
        if best is None or aln.score > best[0]:
            best = (aln.score, strand, q, aln)
    _, strand, q, aln = best
    n = len(seq)
    st = {}
    blocks = list(zip(*aln.aligned))
    for bi, ((qa, qb), (ra, rb)) in enumerate(blocks):
        for k in range(qb - qa):
            a, r = q[qa + k].upper(), reference[ra + k].upper()
            ok = bool(set(IUPAC_SET.get(a, "")) & set(IUPAC_SET.get(r, "")))
            if (k == 0 and bi > 0) or (k == qb - qa - 1 and bi < len(blocks) - 1):
                ok = False
            st[(qa + k + 1) if strand == "+" else (n - qa - k)] = ok
    return st


def reproduced_sites(items, reference, ratio_min=RATIO_MIN, min_reads=2):
    """items: [(name, pred, raw)] from v0_3.classify. Returns (sites, per_read_positions).

    sites: {ref_pos1: {"pair": frozenset, "reads": {name: (read_pos1, primary, secondary, ratio)}}}
    per_read_positions: {name: {read_pos1: frozenset pair}} ready for export.to_sequence(reproduced=...).
    """
    maps = {name: map_to_reference(raw["seq"], reference) for name, _, raw in items}
    sites = {}
    for j in range(len(reference)):
        calls = {}
        for name, pred, raw in items:
            if j not in maps[name]:
                continue
            i, strand = maps[name][j]
            p = pred[i + 1]
            if p["bad"] or p["secondary"]["cls"] == "near_sat":
                continue
            probe = raw["secondary_records"][i]["probe"]
            prim = probe["primary"]
            for b in "ACGT":
                if b == prim or not is_colocated(probe["channels"][b], ratio_min):
                    continue
                pr, sb = (prim, b) if strand == "+" else (COMP[prim], COMP[b])
                calls.setdefault(frozenset((pr, sb)), {})[name] = (i + 1, pr, sb, round(probe["channels"][b]["ratio"], 3))
        ok = [(pair, reads) for pair, reads in calls.items() if len(reads) >= min_reads]
        if ok:
            # several pairs can pass at one site; keep the best supported, strongest one
            ok.sort(key=lambda pr: (len(pr[1]), sum(v[3] for v in pr[1].values())), reverse=True)
            sites[j + 1] = {"pair": ok[0][0], "reads": ok[0][1], "also": [p for p, _ in ok[1:]]}
    per_read = {name: {} for name, _, _ in items}
    for j, s in sites.items():
        for name, (i, _, _, _) in s["reads"].items():
            per_read[name][i] = s["pair"] if maps[name][j - 1][1] == "+" else frozenset(COMP[x] for x in s["pair"])
    return sites, per_read
