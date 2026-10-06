"""Primer-aware ends: find primer sites in a read so primer-derived bases are not reported as sample sequence.

A forward read that runs off the end of the amplicon ends in the reverse
complement of the reverse primer, often followed by one non-templated A.
"""
from Bio.Seq import Seq

IUPAC_SET = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _compatible(a, b):
    return bool(set(IUPAC_SET.get(a.upper(), "")) & set(IUPAC_SET.get(b.upper(), "")))


def find_site(seq, primer, max_mismatch=2):
    """Best (start, end, mismatches) 0-based half-open match of primer in seq, or None."""
    best = None
    L = len(primer)
    for s in range(0, len(seq) - L + 1):
        mm = 0
        for k in range(L):
            if not _compatible(seq[s + k], primer[k]):
                mm += 1
                if mm > max_mismatch:
                    break
        if mm <= max_mismatch and (best is None or mm < best[2] or (mm == best[2] and s > best[0])):
            best = (s, s + L, mm)
    return best


def primer_mask(seq, forward=None, reverse=None, max_mismatch=2):
    """Return (keep_start, keep_end, notes): 1-based inclusive span of insert sequence.

    Looks for the forward/reverse primer at the 5' end and the reverse complement
    of the other primer at the 3' end; works for reads from either primer.
    Everything from a 3' primer site onward (including a non-templated +A) is dropped.
    """
    n = len(seq)
    start, end, notes = 1, n, []
    for name, p in (("forward", forward), ("reverse", reverse)):
        if not p:
            continue
        rc = str(Seq(p).reverse_complement())
        hit5 = find_site(seq, p, max_mismatch)
        if hit5 and hit5[0] < len(p) + 40:
            start = max(start, hit5[1] + 1)
            notes.append(f"{name} primer at 5' {hit5[0]+1}-{hit5[1]} ({hit5[2]} mm)")
        hit3 = find_site(seq, rc, max_mismatch)
        if hit3 and hit3[1] > n - (len(p) + 40):
            end = min(end, hit3[0])
            tail = n - hit3[1]
            notes.append(f"{name} primer RC at 3' {hit3[0]+1}-{hit3[1]} ({hit3[2]} mm); {tail} base(s) after it")
    return start, end, notes
