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
- **H1 starts as a global default** (per-read adaptation: see "Per-read H1" below) (`DEFAULT_ALT`: second peak at x ≈ 0.5, offset s = 0.047, tail weight
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

1. ~~Per-read H1~~ done (below).
2. Real labelled mixtures: colleagues' heterozygote traces; known-template mixes at known ratios.
3. Use the LBF in `export.to_sequence` to decide IUPAC (instead of 0.33) once step 2 supports it.
4. Replace the heterozygous-indel z-threshold with a posterior over (onset, size, no indel).

## Per-read H1 (2026-10-07; the CLI uses this)

`adapt_alt`: self-training (EM) on the read's own calls. Each call gets a responsibility
r = P(two bases | data); H1's parameters are re-estimated from the r-weighted observations, with 30
pseudo-observations pulling them back to the global values, so a read with no real double positions
stays near the global H1. The fraction of two-base calls π is estimated as well, so the posterior
no longer needs a stated prior (LBF still does not depend on it).

Same 210 synthetic mixtures, same leave-one-source-read-out protocol (`EVAL_MODE=em`):

| | AUC global H1 | AUC per-read H1 | ratio alone |
|---|---|---|---|
| all | 0.986 | **0.992** | 0.960 |
| 20 % minor allele, clean | 0.993 | 0.994 | 0.972 |
| noise + jitter | 0.951 | **0.971** | 0.964 |
| LBF > 1 at 20 % allele: sens / spec | 0.68 / 0.997 | 0.83 / 0.989 | rule 0.34 / 0.988 |

Log-loss 0.154 → 0.121 (prior-only 0.687); calibration still close to the diagonal (predicted vs
observed 0.002/0.002, 0.049/0.052, 0.185/0.194, 0.395/0.450, 0.605/0.643, 0.818/0.882,
0.950/0.976, 0.996/0.999). Unmodified reads: estimated π = 0.001 (C1), 0.001 (A2), 0.042 (A3);
LBF > 1 at 1/205, 0/220, 9/211 calls, i.e. no double peaks invented. The 20 %-allele specificity
dropped slightly (0.997 → 0.989); that is the price of the extra sensitivity.

Real reads (Thermocyclops F2, G2): per-read H1 converged to the synthetic shape (second peak at 0.51
of the main, spread 0.48 in log, offset 0.049) without being told; π = 0.068 and 0.058, against
16-17 shared double positions in ~250 calls. Results at the 16 reproduced sites are unchanged
(10 with summed LBF > 2; 204 = 1.43, 251 = 1.41), so the weak sites are weak in the data, not in
the model.

## Genotype output (Clair3's convention; `gt`, `gq`, `pl`, `p_mixed` columns)

Per call, probability over the 10 unordered base pairs (AA, AC, AG, AT, CC, CG, CT, GG, GT, TT, the
order Clair3 uses for its SNP genotypes) as zygosity × identity: homozygous for the primary base with
1 − p_mixed; the heterozygous pairs share p_mixed in proportion to each other channel's co-located
ratio. PL = −10 log10 of the normalised likelihood shifted to 0 for the best genotype (VCF
convention, Clair3's `compute_PL`); GQ = Phred of the probability the call is wrong, capped at 99.
Example, Thermocyclops F2 call 211 (reference 204): read alone AA, GQ 17; the same site is IUPAC `R`
in the FASTA only because the second specimen reproduces it. The FASTA still follows the 0.33 / reproduction
rules; switching to the genotype call is the next step.

## What was taken from Clair3 (source read 2026-10-07, commit eb625c7, BSD-3; no code copied)

| in Clair3 | what it does | used here |
|---|---|---|
| `task/main.py`, `gt21.py`, `genotype.py` | outputs: 21 unordered genotypes + 3-way zygosity + two sorted per-allele length heads (−16…+16) | 10 unordered base-pair genotypes × zygosity; sorted-allele idea noted for the indel decoder's P/Q order |
| `CallVariants.py` `compute_PL`, `quality_score_from` | normalise the product of heads over candidate genotypes → Phred likelihoods; QUAL from P(not reference) | `phred_likelihoods`, `call_genotype` (GQ) |
| `shared/param_*.py` `min_af_dict` | candidate gate depends on platform noise (0.08 HiFi/Illumina, 0.15 ONT) | same principle as the per-read floor and the per-read null |
| pileup model → full-alignment model on low-QUAL candidates | cheap network on everything, expensive one where unsure | v0.3 flags first, indel decoder only on the unflagged span |
| `Train.py` FocalLoss with class-balanced weights (beta 0.999), RAdam | rare classes (het, indel) not drowned by the reference class | not used (no network); relevant if one is ever trained here |
| BiLSTM (pileup) / residual CNN (full alignment) over a 33-position window, 16 flanking bases | learn context effects from millions of labelled sites | not used: our labelled corpus is ~400 calls plus synthetic mixtures |

Not taken, and why: the networks themselves (data), per-strand count channels (a single Sanger
read has no strand replicate; a forward + reverse pair would), haplotype phasing channels (no
Sanger analogue except the two-allele decoder).

## Update 2026-10-07 (after the allergology plates; supersedes the numbers above where they differ)

Two corrections, both found on real reads (details in MULTIREAD.md):

1. **Primary channel.** It was the channel with the highest value anywhere in the slot; a small G after a big C
   then lost to the C's tail and the G looked like a second base. Primary is now the tallest channel with a
   real apex on the call. Synthetic results are unchanged (they reuse a trace's own clean peaks).
2. **Miss rate of a real second peak, q1.** Fitted at 0.001 on synthetic mixtures (the second peak is never
   absent there); on the plates 1 of 20 reads of known heterozygotes showed none. Floor 0.02, default 0.05
   (`Q1_MIN`, `DEFAULT_ALT`). The plates motivated this, so they are not independent evidence for it.

Re-run, same protocols: synthetic, per-read H1: AUC 0.991 (was 0.992), noise + jitter 0.970, log-loss 0.128
(was 0.121; the floor costs a little where nothing is ever missed); unmodified reads unchanged (LBF > 1 at
1/205, 0/220, 9/211). Thermocyclops F2/G2 (per-read H1, π = 0.038 and 0.048): 12 of the 16 reproduced sites
have summed LBF > 1 and 10 have > 2 (none of the other 181 shared positions exceeds 0.1); the weak ones are
204 (1.5), 212 (0.9), 214 (−0.9), 220 (0.9), 251 (1.6) and **269 (−0.9)**, which the primary-channel fix
moved from strong (+2.0 in F2) to against: it was a small G next to a larger A. Sites 204 and 251, the two
that separate Japanese *T. taihokuensis* from Taiwanese *Thermocyclops* sp. 1, stay at about 30:1 each.

## Global-start model (2026-10-07, `adapt_model`; supersedes the per-read fit above)

Problem found on real indel carriers (SINE-flank reads, see HET_INDEL.md): `fit_null` takes the single-base baseline
from the read by median/MAD. That works up to a minority of mixed positions; in a heterozygous-indel read about half
of the positions are mixed, so the baseline moves onto the mixture (median second-peak ratio 0.70 in 60 confident
carriers against 0.06 in 300 cloned, single-allele reads) and the model sees almost nothing (1.3 % of positions above
100:1). The synthetic test did not show it because its baseline came from the clean stretch before the onset.

Fix: H0 starts from values measured on real single-allele reads (300 cloned reads, 191,548 positions: 17 % of
positions carry a weak co-located second apex, median ratio 0.067, offsets uniform) and **both** hypotheses adapt by EM,
each shrunk to the global values by pseudo-observations (30 for H1, 10 for H0). Mixed positions are assigned to H1 by
their responsibilities, so they no longer train the baseline.

| check | read-fitted baseline | global-start model |
|---|---|---|
| 50 decoder-confident indel carriers: share of positions above 100:1 before / after the decoded onset (median) | 0.00 / 0.00 (step in 1/50) | 0.02 / 0.69 (step in 37/50) |
| cloned reads (single allele): positions above 100:1 | 0.14 % | 0.45 % |
| synthetic mixtures, AUC / log-loss, no clean stretch assumed | (needs the clean stretch) | 0.991 / 0.132 |
| synthetic mixtures, overall log-loss | 0.128 (baseline from the clean stretch) | 0.132 (noisy conditions: 0.141 and 0.326; with n0 = 30: 0.79 and 2.59) |
| allergology plates, zygosity vs report | 46/46 joint (biallelic) | 46/46 biallelic; flat prior 45/46 |
| Thermocyclops, 16 reproduced sites with summed LBF > 2 | 10 | 10 |

The pseudo-observation count for H0 was chosen on the synthetic set from three values (30, 10, 3: log-loss 0.716,
0.132, 0.131; with 30 the baseline cannot follow the synthetic white noise); 10 was taken as the middle of the flat
region. The step check on real carriers and the clone control were run afterwards and did not enter that choice.
Cost: more weak calls on clones (0.45 % against 0.14 % of positions above 100:1).
