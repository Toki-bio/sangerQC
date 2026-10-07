"""Evaluate sangerqc.bayes against the old ratio rule, with truth from synthetic mixtures.

    python validation/eval_bayes.py BENCH_DIR OUT_DIR

BENCH_DIR holds the .ab1 files written by validation/bench_hetindel.py (names
READ_k+3_w50_n0_j0.ab1) and the source reads' paths are given below. Truth at a call:
two bases if the two alleles differ there (both A,C,G,T), one base otherwise.

Protocol: leave one source read out; H1 is fitted on the other two reads' mixtures;
H0 is fitted per file from the clean stretch before the indel onset (what the clean prefix of
a real read offers). Reported: AUC, sensitivity/specificity of the old rule and of LBF cut-offs,
calibration of the posterior (prior = the training mixed fraction), and a control on unmodified
reads (every call single).
"""
import csv
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sangerqc import bayes  # noqa: E402
from sangerqc.secondary import secondary_peaks  # noqa: E402
from sangerqc.v0_1 import load_ab1  # noqa: E402

A = "C:/work/Raks_COI/Artemia_parthenogenetica/"
SRC = {"C1": (A + "C1_C1_20241118_143915.ab1", 67, 271), "A2": (A + "A2_A2_20241119_111258.ab1", 66, 285),
       "A3": (A + "A3_A3_20241119_125733.ab1", 200, 410)}


def raw_of(path, lo, hi):
    raw = load_ab1(path)
    hq = float(np.median(raw["top_h"][lo - 1:hi]))
    raw["secondary_records"] = secondary_peaks(raw, hq)
    return raw


def rule_call(obs, ratio_call=0.33):
    return obs is not None and obs[0] >= ratio_call and obs[0] >= 0.20


def auc(scores, labels):
    s = np.asarray(scores, float)
    y = np.asarray(labels, bool)
    order = np.argsort(s)
    ranks = np.empty(len(s))
    ranks[order] = np.arange(1, len(s) + 1)
    # average ranks over ties
    for v in np.unique(s):
        m = s == v
        if m.sum() > 1:
            ranks[m] = ranks[m].mean()
    n1, n0 = y.sum(), (~y).sum()
    return (ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def collect(bench, readkey):
    """Per file: observations in the zone, truth, null-fit positions."""
    path, lo, hi = SRC[readkey]
    src = load_ab1(path)
    seq = src["seq"]
    pos = lo + int(0.4 * (hi - lo))
    rows = []
    for f in sorted(Path(bench).glob(f"{readkey}_k*.ab1")):
        k, w, n, j = (int(x) for x in re.match(rf"{readkey}_k([+-]\d+)_w(\d+)_n(\d+)_j(\d+)", f.name).groups())
        raw = raw_of(str(f), lo, hi)
        obs = bayes.observations(raw)
        d = abs(k)
        zone = range(lo - 1, hi - d - 1)
        truth = {}
        for i in zone:
            a = seq[i]
            b = seq[i] if i < pos - 1 else (seq[i + k] if 0 <= i + k < len(seq) else None)
            if a in "ACGT" and b is not None and b in "ACGT":
                truth[i] = bool(i >= pos - 1 and a != b)
        rows.append(dict(file=f.name, k=k, w=w / 100, noise=n / 100, jitter=j / 100, obs=obs, truth=truth,
                         null_idx=list(range(lo - 1, pos - 1)), pos=pos))
    return rows


def reliability(post, y, bins=(0, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98, 1.0001)):
    post, y = np.asarray(post), np.asarray(y, float)
    out = []
    for a, b in zip(bins[:-1], bins[1:]):
        m = (post >= a) & (post < b)
        if m.sum():
            out.append((a, b, int(m.sum()), float(post[m].mean()), float(y[m].mean())))
    return out


def main(bench, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = {k: collect(bench, k) for k in SRC}
    allrows, summary = [], []
    for held in SRC:
        train = [r for k in SRC if k != held for r in data[k]]
        alt_obs = [r["obs"][i] for r in train for i, t in r["truth"].items() if t]
        alt = bayes.fit_alt(alt_obs)
        n_mixed = sum(t for r in train for t in r["truth"].values())
        n_all = sum(len(r["truth"]) for r in train)
        prior = n_mixed / n_all
        print(f"held out {held}: H1 fitted on {len(alt_obs)} mixed calls: q1={alt.q1:.2f} "
              f"mu1={math.exp(alt.mu1):.2f} s1={alt.s1:.2f} so={alt.so:.3f} eo={alt.eo:.2f}; training mixed fraction {prior:.2f}")
        for r in data[held]:
            lbfs = []
            m = bayes.fit_null([r["obs"][i] for i in r["null_idx"]], alt)
            pi_r = prior
            if os.environ.get("EVAL_MODE") == "em":
                m, pi_r = bayes.adapt_alt([r["obs"][i] for i in r["truth"]], m)
            elif os.environ.get("EVAL_MODE") == "em2":
                # H0 from the global clone baseline, then both hypotheses adapt (no clean prefix assumed)
                base = bayes.Model(**{**alt.to_dict(), **{k: getattr(bayes.DEFAULT_MODEL, k) for k in ("q0", "mu0", "s0")}})
                m, pi_r = bayes.adapt_model([r["obs"][i] for i in r["truth"]], base, n0_h0=float(os.environ.get("N0H0", 30)))
            for i, t in r["truth"].items():
                o = r["obs"][i]
                lbf = bayes.log10_bf(o, m)
                allrows.append(dict(held=held, file=r["file"], k=r["k"], w=r["w"], noise=r["noise"], i=i + 1,
                                    truth=int(t), x=o[0] if o else 0.0, off=o[1] if o else "", lbf=lbf,
                                    post=bayes.posterior(lbf, pi_r if os.environ.get("EVAL_MODE") in ("em", "em2") else prior),
                                    rule=int(rule_call(o))))
    with open(outdir / "eval_positions.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(allrows)

    def report(rows, label):
        y = [r["truth"] for r in rows]
        if sum(y) == 0 or sum(y) == len(y):
            return
        lbf = [r["lbf"] for r in rows]
        xr = [r["x"] for r in rows]
        rule = np.array([r["rule"] for r in rows], bool)
        yb = np.array(y, bool)
        sens = lambda c: (c & yb).sum() / yb.sum()
        spec = lambda c: (~c & ~yb).sum() / (~yb).sum()
        lb = np.array(lbf)
        print(f"{label:34s} n={len(rows):6d} mixed={yb.mean():.2f} | AUC LBF {auc(lbf, y):.3f} vs ratio-only {auc(xr, y):.3f}"
              f" | rule (>=0.33) sens {sens(rule):.2f} spec {spec(rule):.3f}"
              f" | LBF>1 sens {sens(lb > 1):.2f} spec {spec(lb > 1):.3f}"
              f" | LBF>2 sens {sens(lb > 2):.2f} spec {spec(lb > 2):.3f}")

    print()
    report(allrows, "all conditions")
    for w_ in (0.5, 0.3, 0.2):
        report([r for r in allrows if r["w"] == w_ and r["noise"] == 0], f"allele A {int(w_*100)}%, no noise")
    report([r for r in allrows if r["noise"] > 0], "noise 0.1 + jitter 0.3")
    for held in SRC:
        report([r for r in allrows if r["held"] == held], f"held-out read {held}")
    post = [r["post"] for r in allrows]
    y = [r["truth"] for r in allrows]
    p = np.clip(post, 1e-6, 1 - 1e-6)
    ll = -np.mean([math.log(pp) if t else math.log(1 - pp) for pp, t in zip(p, y)])
    base = np.mean(y)
    ll0 = -np.mean([math.log(base) if t else math.log(1 - base) for t in y])
    print(f"\nlog-loss of posterior {ll:.3f} vs prior-only {ll0:.3f}")
    print("reliability (bin, n, mean posterior, observed mixed fraction):")
    for a, b, n, mp, of in reliability(post, y):
        print(f"  [{a:.2f},{min(b,1):.2f}) n={n:6d} predicted {mp:.3f} observed {of:.3f}")

    # control: unmodified reads, every call single; H0 fitted from the whole read (robust)
    print("\ncontrol, unmodified reads (all calls single), null fitted robustly from the whole zone:")
    alt = bayes.fit_alt([r["obs"][i] for k in SRC for r in data[k] for i, t in r["truth"].items() if t])
    for key, (path, lo, hi) in SRC.items():
        raw = raw_of(path, lo, hi)
        obs = bayes.observations(raw)
        m = bayes.fit_null(obs[lo - 1:hi], alt)
        if os.environ.get("EVAL_MODE") == "em2":
            m, pi_c = bayes.adapt_model(obs[lo - 1:hi], bayes.DEFAULT_MODEL)
            print(f"  {key}: estimated fraction of two-base calls pi = {pi_c:.3f}")
        elif os.environ.get("EVAL_MODE") == "em":
            m, pi_c = bayes.adapt_alt(obs[lo - 1:hi], m)
            print(f"  {key}: estimated fraction of two-base calls pi = {pi_c:.3f}")
        lbf = [bayes.log10_bf(obs[i], m) for i in range(lo - 1, hi)]
        rule = [rule_call(obs[i]) for i in range(lo - 1, hi)]
        print(f"  {key}: {sum(l > 1 for l in lbf)}/{len(lbf)} calls LBF>1, {sum(l > 2 for l in lbf)} LBF>2; old rule flags {sum(rule)}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
