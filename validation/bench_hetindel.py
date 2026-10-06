"""Benchmark: synthetic heterozygous indels from clean reads; sangerqc.hetindel vs tracy decompose.

    python validation/bench_hetindel.py OUTDIR TRACY_EXE READ.ab1:START-END [READ.ab1:START-END ...]

START-END = a clean (BLAST-confirmed) call range of the source read, 1-based.
Writes OUTDIR/bench.tsv. Synthetic files and tracy outputs stay in OUTDIR.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sangerqc.hetindel import peak_sets, solve  # noqa: E402
from sangerqc.synth import make_het_indel  # noqa: E402
from sangerqc.v0_3 import classify  # noqa: E402

SIZES = [1, 2, 3, 5, 8, 13, 21]
WEIGHTS = [0.5, 0.3]


def tracy_size(tracy, ab1, ref, prefix):
    subprocess.run([tracy, "decompose", "-r", ref, "-o", prefix, ab1], capture_output=True, check=False)
    j = Path(prefix + ".json")
    if not j.exists():
        return None, None
    d = json.loads(j.read_text())
    dec = d.get("decomposition", {})
    xs, ys = dec.get("x", []), dec.get("y", [])
    nonzero = [(y, x) for x, y in zip(xs, ys) if x != 0]
    size = min(nonzero)[1] if nonzero else None
    return size, d.get("hetindel")


def main(out, tracy, specs):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in specs:
        path, rng = spec.rsplit(":", 1)
        lo, hi = (int(x) for x in rng.split("-"))
        stem = Path(path).stem[:6]
        src_seq = classify(path)[1]["seq"]
        ref = out / f"{stem}_ref.fa"
        ref.write_text(f">{stem}\n{src_seq[lo - 1:hi]}\n")
        pos = lo + int(0.4 * (hi - lo))
        for d in SIZES:
            for sign in (1, -1):
                k = sign * d
                for w in WEIGHTS:
                    name = f"{stem}_k{k:+d}_w{int(w * 100)}"
                    f = out / f"{name}.ab1"
                    t = make_het_indel(path, str(f), pos, k, w)
                    pred, raw = classify(str(f))
                    S = peak_sets(raw)
                    r = solve(S, start=lo - 1, end=hi)
                    A, B = t["allele_a"], t["allele_b"]
                    a0, a1 = pos - 1, hi - d - 2
                    acc = lambda x, y: sum(x[i] == y[i] for i in range(a0, a1)) / max(1, a1 - a0)
                    if r:
                        m1 = (acc(r["P"], A) + acc(r["Q"], B)) / 2
                        m2 = (acc(r["P"], B) + acc(r["Q"], A)) / 2
                        allele_acc = max(m1, m2)
                        ours = dict(found_pos=r["pos"] + 1, found_d=r["d"], rate=round(r["contradiction_rate"], 3),
                                    allele_acc=round(allele_acc, 3))
                    else:
                        ours = dict(found_pos="", found_d="", rate="", allele_acc="")
                    ts, th = tracy_size(tracy, str(f), str(ref), str(out / f"{name}_tracy"))
                    row = dict(read=stem, k=k, weight_a=w, pos=pos, **ours,
                               ours_size_ok=(ours["found_d"] == d), ours_pos_ok=(ours["found_pos"] == pos),
                               tracy_size=ts, tracy_size_ok=(ts is not None and abs(ts) == d), tracy_hetindel=th)
                    rows.append(row)
                    print("\t".join(str(v) for v in row.values()), flush=True)
    with open(out / "bench.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    print(f"ours size correct {sum(r['ours_size_ok'] for r in rows)}/{n}, onset exact {sum(r['ours_pos_ok'] for r in rows)}/{n}, "
          f"mean allele acc {sum(r['allele_acc'] or 0 for r in rows) / n:.3f}; "
          f"tracy size correct {sum(r['tracy_size_ok'] for r in rows)}/{n}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
