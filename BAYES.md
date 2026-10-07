# Bayesian second-peak model (`sangerqc/bayes.py`)

Replaces "ratio ≥ 0.33" with the evidence that a position carries two bases, as NGS genotype
callers do (likelihood of the data under each genotype, then a prior). The ratio rule is kept in
v0.3's IUPAC output; the Bayes factor is reported next to it (`lbf` column of `*_decisions.tsv`
and `reproduced_sites.tsv`) and does not change the FASTA yet.

## Model

Per call the data are the strongest *co-located* second channel (secondary.py geometry): its
height ratio x to the primary and its apex offset |o| in peak spacings, or "none" (x ≥ 0.05).

| | one base (H0) | two bases (H1) |
|---|---|---|
| no co-located second | q0 | q1 (minor peak missed) |
| ratio | log x ~ Normal(μ0, σ0) | log x ~ Normal(μ1, σ1) |
| offset | Uniform(0, 0.35) | (1−e) half-Normal(s) + e Uniform(0, 0.35) |

- **H0 is fitted from the read itself**, robustly (median/MAD), over its unflagged span: a read's
  own spill and noise define "nothing here" for that read. A few real double positions do not
  move a median.
- **H1 is a global default** (`DEFAULT_ALT`: second peak at x ≈ 0.5, offset s = 0.047, tail weight
  e = 0.16, fitted by EM) learned from synthetic mixtures only.
- Output: **LBF = log10 [P(data | two bases) / P(data | one base)]**; independent of any prior.
  LBF 1 = 10:1, 2 = 100:1. Posterior = LBF combined with a prior the *user* states (e.g. 0.05).
- Several reads: add LBFs. Specimens sharing a primer share sequence context, so the sum
  over-counts (the `DS` caveat in SECONDARY_PEAKS.md); forward + reverse would not.

## Results (2026-10-07)

`validation/eval_bayes.py`: 210 synthetic mixtures (three Artemia source reads, allele A share
0.5 / 0.3 / 0.2, with and without noise 0.1 + allele jitter 0.3), 42 010 calls, truth = the two
alleles differ at that call. Leave-one-source-read-out: H1 fitted on two reads, tested on the third;
H0 fitted per file from the clean stretch before the indel.

| | calls | AUC, LBF | AUC, ratio alone | old rule (≥ 0.33): sens / spec | LBF > 1: sens / spec |
|---|---|---|---|---|---|
| all | 42 010 | **0.986** | 0.960 | 0.79 / 0.942 | 0.72 / 0.995 |
| allele A 50 %, clean | 8 402 | 1.000 | 1.000 | 1.00 / 0.988 | 1.00 / 0.997 |
| allele A 30 %, clean | 8 402 | 0.999 | 0.994 | 0.80 / 0.988 | 0.93 / 0.997 |
| allele A 20 %, clean | 8 402 | 0.993 | 0.972 | 0.34 / 0.988 | 0.68 / 0.997 |
| noise + jitter | 16 804 | **0.951** | **0.964** | 0.92 / 0.873 | 0.50 / 0.992 |
| held out C1 / A2 / A3 | ~14 000 each | 0.987 / 0.989 / 0.981 | 0.965 / 0.960 / 0.956 | | |

- **Calibration** (posterior at the training mixed fraction 0.44, 8 bins): predicted vs observed
  mixed fraction 0.002/0.002, 0.050/0.055, 0.185/0.194, 0.396/0.419, 0.602/0.671, 0.817/0.878,
  0.950/0.976, 0.996/0.999: slightly under-confident in the middle, close to the diagonal. Log-loss
  0.154 vs 0.687 for the prior alone.
- **Clean unmodified reads** (every call single): C1 1/205 calls LBF > 1, A2 0/220, A3 9/211
  (7 with LBF > 2; the old rule flags 11 there).
- **Real data, two Thermocyclops specimens** (`validation/eval_bayes_cyclops.py`): of the 16 reference
  positions where both reads have a co-located second peak ≥ 0.20, summed LBF > 1 in 13 and > 2 in
  10; none of the other 181 shared positions exceeds 0.5. Weak or inconclusive: 204 (sum 1.5),
  212, 214, 220, 245, 251 (1.4). 204 and 251 are the two sites that separate Japanese
  *T. taihokuensis* from Taiwanese *Thermocyclops* sp. 1: **evidence for double peaks there,
  about 30:1 each, not conclusive.** (Caveat: the 16 are defined by the same co-location geometry the
  model uses, so this checks specificity and consistency, not independent truth.)

## What it did not fix, and what changed along the way

- **Under noise plus jitter the plain ratio ranks slightly better** (0.951 vs 0.964), because H1
  is a single global shape and noisy mixtures have broader ratios. Per-file H1 (EM on the read's
  own high-ratio calls) is the obvious next step; not done.
- First version (narrow offset under H1) gave AUC 0.943 in that condition. Adding the heavy-tailed
  offset was a structural choice made **after seeing aggregate results**; held-out reads only protect
  the parameter fits, not that choice. The real-read check above was run after.
- Truth comes from synthetic mixtures whose alleles are cut from one real trace: the same dye
  context and peak shapes. Real heterozygotes will break the H1 numbers first.
- Priors are not estimated. The posterior is only as good as the prior you state; quote LBF.

## Next

1. Per-read H1 (adapt μ1, σ1, s to the read) and re-test under noise.
2. Real labelled mixtures: colleagues' heterozygote traces; known-template mixes at known ratios.
3. Use the LBF in `export.to_sequence` to decide IUPAC (instead of 0.33) once step 2 supports it.
4. Replace the heterozygous-indel z-threshold with a posterior over (onset, size, no indel).
