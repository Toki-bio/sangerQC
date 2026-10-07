"""Joint genotype at one reference position from several reads of the same template.

Reads (forward + reverse, re-reads, different primers) are put on a shared reference. Each read
contributes a likelihood over the 10 unordered base-pair genotypes at that position:

    homozygous X     1 if X is the read's primary base, else EPS
    heterozygous XY  10**LBF(read)  x identity term (pair = primary + strongest co-located second: 1;
                     pair contains the primary only: 1/3 of the weight if no second was seen, 0.1 if a
                     different second was seen; otherwise EPS)

The LBF is the read's own prior-free Bayes factor for "two bases" (bayes.py), so a weak, noisy read
says little and a clean read says a lot. Reads of one template are combined as a weighted product,
log L(g) = sum_r w_r log l_r(g), then multiplied by a prior over genotypes:

  - flat prior: P(homozygous X) = (1 - pi) / 4, P(het XY) = pi / 6, with pi the site's heterozygous
    fraction (default: the mean of the reads' own estimates);
  - known biallelic SNP: only AA, AB, BB, with Hardy-Weinberg weights from the allele frequency.

Weights encode how independent the reads are (RELATION_WEIGHT): opposite strands see different
sequence context and count fully; the same primer re-run shares template and context artefacts and
counts less. These are stated assumptions, not fitted values; the allergology plates
(validation/eval_multiread.py) can check the forward+reverse weight only.
"""
import math
from dataclasses import dataclass

from . import bayes
from .bayes import GENOTYPES, adapt_model, log10_bf, observations
from .export import kept_span
from .secondary import COLOCATE_MAX, secondary_peaks
from .v0_3 import classify

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
EPS = 1e-3
RELATION_WEIGHT = {"FR": 1.0, "DP": 1.0, "RR": 0.5, "DM": 0.7, "SAME": 0.5}


def gt_key(a, b):
    return "".join(sorted(a + b))


def comp_gt(g):
    return gt_key(COMP[g[0]], COMP[g[1]])


@dataclass
class ReadAnalysis:
    name: str
    path: str
    raw: dict
    pred: dict
    model: object
    pi: float
    span: tuple


def analyze_read(path, name=None, forward=None, reverse=None):
    """classify (v0.3) + per-read Bayesian model (H0 from the unflagged span, H1 adapted)."""
    pred, raw = classify(path, forward=forward, reverse=reverse)
    obs = observations(raw)
    sp = kept_span(pred, len(raw["seq"]))
    span_obs = obs[sp[0] - 1:sp[1]] if sp else obs
    model, pi = adapt_model(span_obs, bayes.DEFAULT_MODEL)
    raw["bayes"] = dict(pi=pi, model=model)
    return ReadAnalysis(name or path, path, raw, pred, model, pi, sp or (1, len(raw["seq"])))


def read_likelihood(ra: ReadAnalysis, i):
    """{genotype: likelihood} at read index i (0-based), on the READ's strand; also returns diagnostics."""
    probe = ra.raw["secondary_records"][i]["probe"]
    obs = bayes.observe(probe)
    prim = probe["primary"]
    sec = None
    if obs is not None:
        best = -1.0
        for b in "ACGT":
            c = probe["channels"][b]
            if b != prim and c["apex"] and c["offset"] is not None and abs(c["offset"]) <= COLOCATE_MAX \
                    and c["ratio"] >= bayes.X_MIN and c["ratio"] > best:
                best, sec = c["ratio"], b
    lbf = log10_bf(obs, ra.model)
    bf = 10.0 ** max(min(lbf, 12.0), -12.0)
    lik = {}
    for g in GENOTYPES:
        if g[0] == g[1]:
            lik[g] = 1.0 if g[0] == prim else EPS
        else:
            pair = set(g)
            if sec is not None and pair == {prim, sec}:
                ident = 1.0
            elif prim in pair and sec is None:
                ident = 1.0 / 3.0
            elif prim in pair:
                ident = 0.1
            else:
                ident = EPS
            lik[g] = max(bf * ident, 1e-12)
    return lik, dict(prim=prim, sec=sec, lbf=lbf)


def oriented(lik, strand):
    """Move a read-strand likelihood onto the reference strand."""
    return lik if strand == "+" else {comp_gt(g): v for g, v in lik.items()}


def prior_over_genotypes(pi, alleles=None, freq=0.5):
    """Flat (pi het) or known biallelic SNP (Hardy-Weinberg with allele A at `freq`)."""
    if alleles:
        a, b = alleles
        p = {g: 1e-6 for g in GENOTYPES}
        p[gt_key(a, a)] = freq ** 2
        p[gt_key(b, b)] = (1 - freq) ** 2
        p[gt_key(a, b)] = 2 * freq * (1 - freq)
    else:
        p = {g: ((1 - pi) / 4 if g[0] == g[1] else pi / 6) for g in GENOTYPES}
    s = sum(p.values())
    return {g: v / s for g, v in p.items()}


def joint_genotype(items, pi=0.05, alleles=None, freq=0.5, weights=None):
    """items: [(ReadAnalysis, read_index0, strand, relation)] -> dict(posterior, call, gq, per-read)

    Relation names the read's role for weighting ("FR", "DP", "RR", "DM", "SAME"); `weights` overrides
    RELATION_WEIGHT. Reads flagged bad at that position are skipped (their evidence is not usable)."""
    wts = dict(RELATION_WEIGHT)
    wts.update(weights or {})
    logL = {g: 0.0 for g in GENOTYPES}
    used, per_read = 0, []
    for ra, i, strand, rel in items:
        if ra.pred[i + 1]["bad"]:
            per_read.append((ra.name, i + 1, "skipped (flagged bad)"))
            continue
        lik, diag = read_likelihood(ra, i)
        lik = oriented(lik, strand)
        w = wts.get(rel, 1.0)
        for g in GENOTYPES:
            logL[g] += w * math.log(lik[g])
        used += 1
        per_read.append((ra.name, i + 1, diag))
    prior = prior_over_genotypes(pi, alleles, freq)
    lp = {g: math.log(prior[g]) + logL[g] for g in GENOTYPES}
    m = max(lp.values())
    z = sum(math.exp(v - m) for v in lp.values())
    post = {g: math.exp(v - m) / z for g, v in lp.items()}
    call, gq = bayes.call_genotype(post)
    return dict(posterior=post, call=call, gq=gq, pl=bayes.phred_likelihoods(post), reads_used=used,
                p_het=sum(v for g, v in post.items() if g[0] != g[1]), per_read=per_read)
