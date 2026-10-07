"""Bayesian second-peak model on two real reads (Thermocyclops F2, G2; different specimens, same primer).

    python validation/eval_bayes_cyclops.py BENCH_DIR F2.ab1 G2.ab1 REF.fasta

H1 comes from the synthetic mixtures (all three source reads, pooled); H0 is fitted per read from
its own kept span, robustly. Per reference coordinate the log10 Bayes factors of the two reads are
added (independent specimens, same sequence context, so this slightly over-counts: the DS caveat in
SECONDARY_PEAKS.md). Reference for "truth": sites where BOTH reads show a co-located second peak
(the 16 manual/automatic reproduced sites). That is evidence, not proof, so the comparison reports
how the two ways of using the same two reads rank those sites against every other position.
"""
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import eval_bayes as E  # noqa: E402
from sangerqc import bayes  # noqa: E402
from sangerqc.concordance import map_to_reference  # noqa: E402
from sangerqc.secondary import secondary_peaks  # noqa: E402
from sangerqc.v0_1 import load_ab1  # noqa: E402

SPANS = {"F2": (31, 284), "G2": (52, 248)}   # kept spans from sangerqc.cli (primers + reference trim)


def main(bench, f2, g2, ref_fa):
    data = {k: E.collect(bench, k) for k in E.SRC}
    alt = bayes.fit_alt([r["obs"][i] for k in E.SRC for r in data[k] for i, t in r["truth"].items() if t])
    print("H1 (pooled synthetic):", {k: round(v, 3) for k, v in alt.to_dict().items()})
    ref = str(next(SeqIO.parse(ref_fa, "fasta")).seq).upper()
    res = {}
    for name, path in (("F2", f2), ("G2", g2)):
        raw = load_ab1(path)
        lo, hi = SPANS[name]
        hq = float(np.median(raw["top_h"][lo - 1:hi]))
        raw["secondary_records"] = secondary_peaks(raw, hq)
        obs = bayes.observations(raw)
        m = bayes.fit_null(obs[lo - 1:hi], alt)
        print(f"{name}: null q0={m.q0:.2f} median x={np.exp(m.mu0):.2f} s0={m.s0:.2f}")
        mp = map_to_reference(raw["seq"], ref)             # ref index0 -> (read index0, strand)
        lbf = {j: bayes.log10_bf(obs[i], m) for j, (i, _) in mp.items() if lo - 1 <= i <= hi - 1}
        rule = {j: int(E.rule_call(obs[i])) for j, (i, _) in mp.items() if lo - 1 <= i <= hi - 1}
        res[name] = (lbf, rule, obs, mp)
    common = sorted(set(res["F2"][0]) & set(res["G2"][0]))
    both_rule = {j for j in common if res["F2"][1][j] and res["G2"][1][j]}
    both_coloc = set()
    for j in common:
        o1 = res["F2"][2][res["F2"][3][j][0]]
        o2 = res["G2"][2][res["G2"][3][j][0]]
        if o1 and o2 and o1[0] >= 0.20 and o2[0] >= 0.20:
            both_coloc.add(j)
    summed = {j: res["F2"][0][j] + res["G2"][0][j] for j in common}
    print(f"\ncommon reference positions: {len(common)}; both reads co-located >=0.20: {len(both_coloc)}; "
          f"both pass the 0.33 rule: {len(both_rule)}")
    ref_set = both_coloc
    others = [j for j in common if j not in ref_set]
    for label, score in (("sum of LBF (F2+G2)", summed),
                         ("F2 LBF alone", res["F2"][0]), ("G2 LBF alone", res["G2"][0])):
        s_in = np.array([score[j] for j in ref_set])
        s_out = np.array([score[j] for j in others])
        thr = [1, 2, 3]
        print(f"{label:20s} reference sites: median {np.median(s_in):.2f} min {s_in.min():.2f} | "
              f"others: median {np.median(s_out):.2f} max {s_out.max():.2f} | "
              + " ".join(f"LBF>{t}: sites {int((s_in > t).sum())}/{len(s_in)}, others {int((s_out > t).sum())}/{len(s_out)}"
                         for t in thr))
    print("\nreference positions whose summed LBF > 2 but are not in the reproduced set:",
          sorted(j + 1 for j in others if summed[j] > 2))
    print("reproduced sites with summed LBF <= 2:", sorted((j + 1, round(summed[j], 2)) for j in ref_set if summed[j] <= 2))


if __name__ == "__main__":
    main(*sys.argv[1:5])
