"""Synthetic heterozygous-indel traces with known truth, made from a real clean read.

Allele A is the read as recorded. Allele B carries an indel at call `pos`
(1-based): for every call i >= pos, B shows the peak of A's call i + k.
  k > 0: B lacks A's calls pos .. pos+k-1 (a k-bp deletion)
  k < 0: B repeats A's |k| calls before pos (a |k|-bp tandem duplication)
B's trace is A's trace warped peak-interval by peak-interval, so peak shape,
spacing and dye behaviour stay real. The mixed trace w*A + (1-w)*B is written
into a byte copy of the .ab1 (DATA9-12 only; calls and PLOC are left as they
were, as a basecaller run on a mixed template would leave them roughly).
"""
import struct

import numpy as np

from .v0_1 import load_ab1


def _abif_entries(buf):
    if buf[:4] != b"ABIF":
        raise ValueError("not an ABIF file")
    # root directory entry starts at byte 6
    _, _, _, _, n, _, off, _ = struct.unpack(">4sihhiiii", buf[6:34])
    entries = {}
    for j in range(n):
        e = buf[off + 28 * j: off + 28 * (j + 1)]
        name, num, etype, esize, nelem, dsize, doff, _ = struct.unpack(">4sihhiiii", e)
        entries[(name.decode("latin1"), num)] = (etype, esize, nelem, dsize, doff)
    return entries


def allele_b(channels, ploc, pos, k):
    """Warp A's channels so call i (>= pos) shows A's call i+k."""
    n_scans = len(channels["A"])
    ploc = np.asarray(ploc, dtype=float)
    n = len(ploc)
    out = {b: channels[b].astype(float).copy() for b in "ACGT"}
    start = int(round((ploc[pos - 2] + ploc[pos - 1]) / 2)) if pos >= 2 else 0
    scans = np.arange(start, n_scans)
    # fractional call index of every scan, then shift by k and map back to a scan of A
    call_idx = np.interp(scans, ploc, np.arange(n), left=0, right=n - 1)
    # beyond the last call, extrapolate with mean spacing
    sp = float(np.mean(np.diff(ploc)))
    tail = scans > ploc[-1]
    call_idx[tail] = (n - 1) + (scans[tail] - ploc[-1]) / sp
    src_idx = call_idx + k
    src_scan = np.interp(src_idx, np.arange(n), ploc)
    valid = (src_idx >= 0) & (src_idx <= n - 1)
    for b in "ACGT":
        v = np.interp(src_scan, np.arange(n_scans), channels[b].astype(float))
        v[~valid] = 0.0
        out[b][start:] = v
    return out


def make_het_indel(ab1_in, ab1_out, pos, k, weight_a=0.5):
    raw = load_ab1(ab1_in)
    ch, ploc, seq = raw["channels"], raw["ploc"], raw["seq"]
    b = allele_b(ch, ploc, pos, k)
    mixed = {x: weight_a * ch[x].astype(float) + (1 - weight_a) * b[x] for x in "ACGT"}
    buf = bytearray(open(ab1_in, "rb").read())
    ent = _abif_entries(bytes(buf))
    order = raw["base_order"]
    for tag_num, base in zip((9, 10, 11, 12), order):
        etype, esize, nelem, dsize, doff = ent[("DATA", tag_num)]
        assert esize == 2 and nelem == len(mixed[base])
        arr = np.clip(np.round(mixed[base]), -32768, 32767).astype(">i2")
        buf[doff:doff + dsize] = arr.tobytes()
    open(ab1_out, "wb").write(bytes(buf))
    # truth: allele sequences in call space
    a_seq = seq
    n = len(seq)
    b_seq = seq[:pos - 1] + "".join(seq[i + k] if 0 <= i + k < n else "" for i in range(pos - 1, n))
    return dict(pos=pos, k=k, weight_a=weight_a, allele_a=a_seq, allele_b=b_seq)
