"""Benchmark: synthetic heterozygous indels from clean reads.
sangerqc.hetindel (reference-free and reference-guided) vs tracy decompose.

    python validation/bench_hetindel.py OUTDIR TRACY_EXE READ.ab1:START-END [...] [--neg NEG.ab1 ...]

START-END = a clean (BLAST-confirmed) call range of the source read, 1-based; the first
START-END stretch is also the reference given to the reference-guided mode and to Tracy.
--neg: reads that contain no indel shift (they must come back "no indel").
Writes OUTDIR/bench.tsv and OUTDIR/negatives.tsv.
"""
import csv
import json
import os
import subprocess
import zlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sangerqc.hetindel import peak_sets, solve, solve_reference  # noqa: E402
from sangerqc.synth import make_het_indel  # noqa: E402
from sangerqc.v0_3 import classify  # noqa: E402

SIZES = [1, 2, 3, 5, 8, 13, 21]
# (weight of allele A, noise, jitter)
CONDITIONS = [(0.5, 0.0, 0.0), (0.3, 0.0, 0.0), (0.2, 0.0, 0.0), (0.5, 0.1, 0.3), (0.3, 0.1, 0.3)]
if os.environ.get("BENCH_QUICK"):   # end-to-end smoke test
    SIZES, CONDITIONS = [3], [(0.5, 0.0, 0.0), (0.3, 0.1, 0.3)]


def tracy_call(tracy, ab1, ref, prefix):
    subprocess.run([tracy, "decompose", "-r", ref, "-o", prefix, ab1], capture_output=True, check=False)
    j = Path(prefix + ".json")
    if not j.exists():
        return None, None
    d = json.loads(j.read_text())
    dec = d.get("decomposition", {})
    xs, ys = dec.get("x", []), dec.get("y", [])
    nonzero = [(y, x) for x, y in zip(xs, ys) if x != 0]
    size = min(nonzero)[1] if nonzero else None
    # tracy reports a het indel when some shift beats the unshifted fit
    return size, d.get("hetindel")


def acc_pair(P, Q, A, B, a0, a1):
    acc = lambda x, y: sum(x[i] == y[i] for i in range(a0, a1)) / max(1, a1 - a0)
    return max((acc(P, A) + acc(Q, B)) / 2, (acc(P, B) + acc(Q, A)) / 2)


def main(argv):
    out, tracy = Path(argv[0]), argv[1]
    rest = argv[2:]
    negs = []
    if "--neg" in rest:
        k = rest.index("--neg")
        rest, negs = rest[:k], rest[k + 1:]
    out.mkdir(parents=True, exist_ok=True)
    rows, nrows = [], []
    for spec in rest:
        path, rng = spec.rsplit(":", 1)
        lo, hi = (int(x) for x in rng.split("-"))
        stem = Path(path).stem[:2]
        src_seq = classify(path)[1]["seq"]
        refseq = src_seq[lo - 1:hi]
        ref = out / f"{stem}_ref.fa"
        ref.write_text(f">{stem}\n{refseq}\n")
        pos = lo + int(0.4 * (hi - lo))
        for (w, noise, jit) in CONDITIONS:
            for d in SIZES:
                for sign in (1, -1):
                    k = sign * d
                    name = f"{stem}_k{k:+d}_w{int(w * 100)}_n{int(noise * 100)}_j{int(jit * 100)}"
                    f = out / f"{name}.ab1"
                    t = make_het_indel(path, str(f), pos, k, w, noise, jit, seed=zlib.crc32(name.encode()))
                    pred, raw = classify(str(f))
                    S = peak_sets(raw)
                    A, B = t["allele_a"], t["allele_b"]
                    a0, a1 = pos - 1, hi - d - 2
                    rf = solve(S, start=lo - 1, end=hi, raw=raw)
                    rr = solve_reference(S, refseq, raw["seq"], start=lo - 1, end=hi, raw=raw)
                    ts, th = tracy_call(tracy, str(f), str(ref), str(out / f"{name}_tracy"))
                    row = dict(read=stem, k=k, w=w, noise=noise, jitter=jit,
                               free_indel=bool(rf and rf["indel"]), free_d=rf["d"] if rf else "",
                               free_acc=round(acc_pair(rf["P"], rf["Q"], A, B, a0, a1), 3) if rf else "",
                               ref_indel=bool(rr and rr["indel"]), ref_shift=rr["shift"] if rr else "",
                               ref_acc=round(acc_pair(rr["P"], rr["Q"], A, B, a0, a1), 3) if rr else "",
                               ref_frac=round(rr["allele_fraction_Q"], 2) if rr and rr.get("allele_fraction_Q") else "",
                               tracy_size=ts, tracy_het=th)
                    row["free_ok"] = row["free_indel"] and row["free_d"] == d
                    row["ref_ok"] = row["ref_indel"] and row["ref_shift"] == k   # k>0 deletion, k<0 insertion
                    row["tracy_ok"] = ts is not None and abs(ts) == d and th == 1
                    rows.append(row)
                    print("\t".join(str(v) for v in row.values()), flush=True)
    for f in negs:
        pred, raw = classify(f)
        S = peak_sets(raw)
        rf = solve(S, start=50, end=len(S) - 5, raw=raw)
        nrows.append(dict(read=Path(f).name, free_indel=bool(rf and rf["indel"]), free_d=rf["d"] if rf else ""))
        print("NEG", nrows[-1], flush=True)
    with open(out / "bench.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    if nrows:
        with open(out / "negatives.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(nrows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(nrows)
    print("\ncondition (w, noise, jitter): n | free size+call | ref shift+call | tracy size | free acc | ref acc")
    for c in CONDITIONS:
        rs = [r for r in rows if (r["w"], r["noise"], r["jitter"]) == c]
        mean = lambda key: sum(r[key] for r in rs if r[key] != "") / max(1, sum(1 for r in rs if r[key] != ""))
        print(c, len(rs), "|", sum(r["free_ok"] for r in rs), "|", sum(r["ref_ok"] for r in rs), "|",
              sum(r["tracy_ok"] for r in rs), "| %.3f | %.3f" % (mean("free_acc"), mean("ref_acc")))
    if nrows:
        print("negatives called as indel (free):", sum(r["free_indel"] for r in nrows), "/", len(nrows))


if __name__ == "__main__":
    main(sys.argv[1:])
