"""Bayesian second-peak model: P(position carries two bases | the four channels), instead of a ratio cut-off.

For every call the evidence is the strongest *co-located* second channel (secondary.py geometry):
its height ratio x to the primary and its apex offset |o| (in peak spacings), or "none".

  H0  one base:   none with prob q0; else log x ~ Normal(mu0, s0) and |o| ~ Uniform(0, COLOCATE_MAX)
                  (a second channel that happens to peak near the called peak)
  H1  two bases:  none with prob q1 (minor peak missed); else log x ~ Normal(mu1, s1) and
                  |o| ~ (1-eo) half-Normal(so) + eo Uniform(0, COLOCATE_MAX): real peaks sit on the
                  called peak, but in a noisy trace the apex wanders (heavy tail, weight fitted by EM)

H0 is estimated from the file's own single-base positions (robustly, so a few real double
positions do not matter); H1 is a global default fitted on labelled mixtures. The output that
does not depend on any prior is log10 Bayes factor (LBF); the posterior additionally needs a prior
pi that the caller states. Evidence over several reads adds in log space (correlated reads, e.g.
same primer, should be down-weighted by the caller).

Estimated from 2026-10 synthetic mixtures only; see validation/eval_bayes.py before relying on it.
"""
import math
from dataclasses import dataclass, asdict

import numpy as np

from .secondary import COLOCATE_MAX

X_MIN = 0.05   # ignore co-located seconds weaker than this (below the trace's own wobble)
# Probability that a real second peak is missed (allele dropout / strong imbalance in one read).
# Synthetic mixtures never miss it (fit gives ~0), real heterozygote reads do: on the allergology
# plates 1 of 20 reads of known heterozygotes showed no second peak. Floor at 2%, default 5%.
Q1_MIN = 0.02
S_MIN = 0.5    # floor for the null's spread in log x: keeps the null from claiming certainty


@dataclass
class Model:
    q0: float = 0.9
    mu0: float = math.log(0.12)
    s0: float = 0.7
    q1: float = 0.05
    mu1: float = math.log(0.5)
    s1: float = 0.6
    so: float = 0.1
    eo: float = 0.0

    def to_dict(self):
        return asdict(self)


# H1 fitted on 2026-10 synthetic mixtures (3 Artemia source reads, allele A share 0.2-0.5, with and
# without noise); real heterozygotes are not in it. Treat the numbers as a starting point.
DEFAULT_ALT = Model(q0=0.9, mu0=math.log(0.12), s0=0.7, q1=0.05, mu1=-0.68, s1=0.49, so=0.047, eo=0.16)


def observe(probe, x_min=X_MIN, colocate_max=COLOCATE_MAX):
    """(x, |offset|) of the strongest co-located non-primary channel, or None."""
    best = None
    for b in "ACGT":
        if b == probe["primary"]:
            continue
        c = probe["channels"][b]
        if c["apex"] and c["offset"] is not None and abs(c["offset"]) <= colocate_max and c["ratio"] >= x_min:
            if best is None or c["ratio"] > best[0]:
                best = (float(c["ratio"]), abs(float(c["offset"])))
    return best


def observations(raw):
    return [observe(r["probe"]) for r in raw["secondary_records"]]


def _lognorm(lx, mu, s):
    return -0.5 * ((lx - mu) / s) ** 2 - math.log(s) - 0.5 * math.log(2 * math.pi)


def _halfnorm_trunc(o, so, cmax=COLOCATE_MAX):
    norm = math.erf(cmax / (so * math.sqrt(2)))
    return math.log(2.0 / (so * math.sqrt(2 * math.pi)) * math.exp(-0.5 * (o / so) ** 2) / norm + 1e-300)


def loglik(obs, m: Model, hyp: int):
    q = m.q0 if hyp == 0 else m.q1
    if obs is None:
        return math.log(max(q, 1e-9))
    x, o = obs
    lx = math.log(x)
    if hyp == 0:
        return math.log(max(1 - q, 1e-9)) + _lognorm(lx, m.mu0, m.s0) - math.log(COLOCATE_MAX)
    off = math.log((1 - m.eo) * math.exp(_halfnorm_trunc(o, m.so)) + m.eo / COLOCATE_MAX)
    return math.log(max(1 - q, 1e-9)) + _lognorm(lx, m.mu1, m.s1) + off


def log10_bf(obs, m: Model):
    """log10 [ P(obs | two bases) / P(obs | one base) ]; > 0 favours two bases."""
    return (loglik(obs, m, 1) - loglik(obs, m, 0)) / math.log(10)


def posterior(lbf, prior=0.05):
    z = lbf * math.log(10) + math.log(prior / (1 - prior))
    return 1.0 / (1.0 + math.exp(-z))


def fit_alt(obs_list, base=None):
    """Fit H1 from observations at positions known to carry two bases."""
    m = Model(**(base.to_dict() if base else {}))
    seen = [o for o in obs_list if o is not None]
    m.q1 = max(1 - len(seen) / max(1, len(obs_list)), Q1_MIN)
    if len(seen) >= 10:
        lx = np.log([o[0] for o in seen])
        m.mu1, m.s1 = float(np.mean(lx)), float(max(np.std(lx), 0.2))
        offs = np.array([o[1] for o in seen])
        so, eo = 0.08, 0.3
        for _ in range(50):   # EM for the two-component offset mixture
            h = np.exp([_halfnorm_trunc(o, so) for o in offs])
            r = (1 - eo) * h / ((1 - eo) * h + eo / COLOCATE_MAX)
            so = float(max(math.sqrt(np.sum(r * offs ** 2) / max(np.sum(r), 1e-9)), 0.02))
            eo = float(min(max(1 - np.mean(r), 0.0), 0.9))
        m.so, m.eo = so, eo
    return m


def fit_null(obs_list, base: Model):
    """Fit H0 robustly from a stretch of (mostly) single-base positions."""
    m = Model(**base.to_dict())
    seen = [o for o in obs_list if o is not None]
    m.q0 = min(max(1 - len(seen) / max(1, len(obs_list)), 0.02), 0.98)
    if len(seen) >= 8:
        lx = np.log([o[0] for o in seen])
        med = float(np.median(lx))
        mad = float(np.median(np.abs(lx - med))) * 1.4826
        m.mu0, m.s0 = med, max(mad, S_MIN)
    return m


def score_read(raw, null_positions=None, alt=None, prior=0.05):
    """Per-position LBF and posterior for a read.

    null_positions: indices (0-based) of calls used to fit the file's own H0; default: all calls
    (the fit is robust to a minority of real double positions).
    """
    alt = alt or Model()
    obs = observations(raw)
    idx = list(range(len(obs))) if null_positions is None else list(null_positions)
    m = fit_null([obs[i] for i in idx], alt)
    out = []
    for o in obs:
        lbf = log10_bf(o, m)
        out.append((lbf, posterior(lbf, prior)))
    return out, m, obs


def adapt_alt(obs_list, m: Model, n0=30.0, iters=30, prior0=0.05):
    """Per-read H1 by self-training (EM) with shrinkage toward the global H1 held in `m`.

    H0 in `m` stays fixed (fitted from the read's own single-base stretch). Each position gets a
    responsibility r_i = P(two bases | data); H1's parameters are re-estimated from the
    r-weighted observations, with `n0` pseudo-observations pulling them back to the global
    values, so a read with no real double positions keeps the global H1 instead of drifting
    onto the tail of H0. Returns (adapted Model, pi), pi = fraction of calls with two bases.
    """
    base = m
    cur = Model(**m.to_dict())
    pi = prior0
    for _ in range(iters):
        r = []
        for o in obs_list:
            lbf = (loglik(o, cur, 1) - loglik(o, cur, 0))
            z = lbf + math.log(pi / (1 - pi))
            r.append(1.0 / (1.0 + math.exp(-max(min(z, 50), -50))))
        r = np.array(r)
        seen = [(i, o) for i, o in enumerate(obs_list) if o is not None]
        w_all = r.sum()
        w_seen = np.array([r[i] for i, _ in seen])
        lx = np.array([math.log(o[0]) for _, o in seen])
        offs = np.array([o[1] for _, o in seen])
        W = w_seen.sum()
        cur.q1 = float(min(max((n0 * base.q1 + (w_all - W)) / (n0 + w_all), Q1_MIN), 0.9))
        cur.mu1 = float((n0 * base.mu1 + np.sum(w_seen * lx)) / (n0 + W))
        var = (n0 * (base.s1 ** 2) + np.sum(w_seen * (lx - cur.mu1) ** 2)) / (n0 + W)
        cur.s1 = float(max(math.sqrt(var), 0.2))
        h = np.array([(1 - cur.eo) * math.exp(_halfnorm_trunc(o, cur.so)) /
                      ((1 - cur.eo) * math.exp(_halfnorm_trunc(o, cur.so)) + cur.eo / COLOCATE_MAX) for o in offs])
        cur.so = float(max(math.sqrt((n0 * base.so ** 2 + np.sum(w_seen * h * offs ** 2)) /
                                      (n0 + np.sum(w_seen * h))), 0.02))
        cur.eo = float(min(max((n0 * base.eo + np.sum(w_seen * (1 - h))) / (n0 + W), 0.0), 0.9))
        pi = float(min(max(r.mean(), 1e-3), 0.9))
    return cur, pi


# ---- genotype probabilities and Phred likelihoods (the Clair3 output convention) --------------
GENOTYPES = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"]   # Clair3's GT21 SNP order


def _gt(a, b):
    return "".join(sorted(a + b))


def genotype_probs(probe, p_mixed, x_min=X_MIN):
    """P over the 10 unordered base pairs at one call, as a zygosity x identity product:
    homozygous X = (1 - p_mixed) for the primary base; heterozygous pairs share p_mixed in
    proportion to each other channel's co-located ratio (the strongest co-located second
    channel dominates)."""
    prim = probe["primary"]
    w = {}
    for b in "ACGT":
        if b == prim:
            continue
        c = probe["channels"][b]
        r = c["ratio"] if (c["apex"] and c["offset"] is not None and abs(c["offset"]) <= COLOCATE_MAX) else c["ratio"] * 0.1
        w[b] = max(r, x_min * 0.1)
    tot = sum(w.values())
    probs = {g: 1e-6 for g in GENOTYPES}
    probs[_gt(prim, prim)] = max(1 - p_mixed, 1e-6)
    for b, v in w.items():
        probs[_gt(prim, b)] = max(p_mixed * v / tot, 1e-6)
    s = sum(probs.values())
    return {g: v / s for g, v in probs.items()}


def phred_likelihoods(probs):
    """PL as in VCF / Clair3 compute_PL: -10 log10 of normalised likelihood, shifted so the best is 0."""
    pl = {g: -10.0 * math.log10(p + 1e-8) for g, p in probs.items()}
    best = min(pl.values())
    return {g: int(math.ceil(v - best)) for g, v in pl.items()}


def call_genotype(probs):
    """(genotype, GQ): GQ = Phred of the probability that the call is wrong, capped at 99."""
    g = max(probs, key=probs.get)
    return g, int(min(99, round(-10.0 * math.log10(max(1 - probs[g], 1e-10)))))
