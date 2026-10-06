# Heterozygous indels: re-phasing a mixed trace (prototype)

When a template carries two alleles that differ by an indel, the trace is clean up to the
indel and then shows two superimposed sequences, one shifted by d bases against the other,
to the end of the read. The caller turns that into a run of IUPAC letters and wrong bases.
The task is to recover both allele sequences and the indel.

## Prior methods (verified)

| method | input | how |
|---|---|---|
| Indelligent — Dmitriev & Rakitov 2008, *PLoS Comput Biol* | one mixed trace, **no reference** | dynamic programming over the string of peak letters |
| Poly Peak Parser — Hill et al. 2014, *Dev Dyn* | mixed trace + reference | separates wild-type and mutant from ambiguous calls |
| Tracy `decompose` — Rausch et al. 2020, *BMC Genomics* | mixed trace + reference fasta or **wild-type .ab1** | decomposition score per shift, allele fractions, allele-specific alignments; v0.9.1 is in the Artemia folder |

All three start from "which bases have a peak at this call". A neighbour's tail (`spill`) puts a
false second base into that set; in the Thermocyclops reads ~80 % of ratio ≥ 0.33 secondaries
were spill. sangerQC's co-location test (secondary.py) removes them, so it is a better front end
for any of these decoders, including ours.

## What was taken from each, after reading the sources (2026-10-06)

No code was copied (Tracy is BSD-3, sangerseqR GPL-2); the ideas below were re-implemented.

| source | read | idea taken | where |
|---|---|---|---|
| Tracy `src/decompose.h` (`findBreakpoint`, `decomposeAlleles`) | full source | breakpoint from a **continuous** mixing signal, step between two 25-call windows; reference-guided scan of deletion *and* insertion lengths counting calls the reference cannot explain; size chosen by a **median − MAD** cut-off, smallest size first; allele fraction from peak heights | `onset_continuous`, `solve_reference`, `size_call`, `allele_fraction` |
| sangerseqR `makeBaseCalls`, `getpeaks`, `setAllelePhase` (Poly Peak Parser's engine) | full source | confirms the slot (midpoint) window, true local maxima with plateau handling, the 0.33 ratio; phasing = reference base to one allele, remainder to the other | same conventions in secondary.py; reference mode |
| Indelligent paper (Dmitriev & Rakitov 2008; full text via Europe PMC — the web site and its source no longer respond) | methods section | **large shifts fit by chance** (few overlapping calls), so prefer the smallest adequate shift; SNPs between alleles and several indels handled by DP over per-site phase shifts | smallest-size rule in `size_call`; multi-indel DP **not** implemented |

Two further fixes came from the harder benchmark, not from the sources: at most **two bases per
call** (a two-allele mixture cannot have more), and a **per-read floor** for counting a secondary,
taken from the clean stretch before the onset (noise there otherwise makes 59 % of calls "double"
and lets a wrong large shift win).

## Model used here (`sangerqc/hetindel.py`)

`S_i` = bases with a peak on call i (primary + co-located secondaries ≥ 0.20).
Before the onset both alleles agree; from the onset `pos`, one allele lags the other by d:

    S_i = {P_i, Q_i},   Q_i = P_{i-d}   (i ≥ pos + d)

Calls split into d independent chains (`pos + r + j·d`). Fixing P at a chain's head fixes the
whole chain, so each chain has two phasings; keep the one with fewer contradictions
(`Q_i ∉ S_i`). Scan d = 1..40 and the onset ±6 calls around the double-peak step; keep the
lowest contradiction rate. Linear in read length × d, no reference needed. Which allele is the
longer one (deletion vs insertion) is then named against a reference (`name_variant`).

Onset = the call that maximises (double-peak fraction from there to the end) − (fraction in the
30 calls before), restricted to the read's unflagged span, so a noisy 5′ start is not mistaken for
the indel. In repeats the position of an indel is not unique (left/right-alignment); the decoder
reports where the two alleles first differ in the trace.

## Synthetic test bed (`sangerqc/synth.py`)

Takes a clean real read as allele A, warps its own trace peak-interval by peak-interval to make
allele B with a known deletion (k > 0) or tandem duplication (k < 0), and writes
`w·A + (1−w)·B` into a byte copy of the .ab1. Peak shapes, spacing and dye behaviour stay real,
and Tracy reads the files as ordinary traces (it recovered a synthetic 3-bp deletion exactly).

`validation/bench_hetindel.py` runs a grid (d = 1, 2, 3, 5, 8, 13, 21; deletion and duplication;
w = 0.5 and 0.3) on clean zones of several reads and compares with `tracy decompose`.
Results: see "Benchmark" below.

## Benchmark

`validation/bench_hetindel.py`, 2026-10-06: synthetic indels at 40 % of the clean zone of three
Artemia reads (C1 calls 67–271, A2 66–285, A3 200–410); d = 1, 2, 3, 5, 8, 13, 21 as deletion and
duplication; 42 cases per condition. "Correct" = indel called **and** size right (reference mode:
size and deletion/insertion right; Tracy: size right and `hetindel` flag set).

| condition (allele A share, noise, jitter) | reference-free, first version | reference-free, final | reference-guided, final | Tracy decompose |
|---|---|---|---|---|
| 50 %, none | 42 | 42 | 42 | 42 |
| 30 %, none | 42 | 42 | 42 | 42 |
| 20 %, none | 32 | 41 | 42 | 22 |
| 50 %, noise 0.1, jitter 0.3 | **0** | 42 | 42 | 42 |
| 30 %, noise 0.1, jitter 0.3 | **1** | 42 | 42 | 42 |
| **total / 210** | 117 | 209 | 210 | 190 |

Allele-sequence accuracy (final): reference-guided 0.992–1.000 per condition; reference-free
0.83–0.97 (lowest with noise and a 30 % allele). Through the automatic entry point
`decode_read()` (span chosen by the read's own flags, no hand-given range): 209/210 and 210/210.

Negatives (must return "no indel"): Thermocyclops F2 and G2 (a block of double peaks that ends
again), and Artemia C1, A2, A3, B3, G4, H4: **0/8 called**. The first final run called H4 a 1-bp
indel at call 807 — inside H4's BLAST-confirmed degraded tail; `decode_read` now stays inside the
unflagged span and the call disappears.

What changed between "first version" and "final": two bases per call at most; per-read floor for
secondaries tried alongside the plain 0.20 floor (keep the more significant answer — the floor
alone cost weak 20 % alleles); continuous breakpoint signal with a persistence check; size chosen
by the median − MAD rule with the smallest adequate size; decoding only inside the unflagged span.

Limits of this test: one indel per read, no SNPs between the alleles, both alleles cut from the
same real trace, white noise. It shows the logic works and where it breaks; it is not evidence on
real heterozygotes.

## What real data is needed next

Synthetic mixtures are too clean in one way that matters: both alleles have identical peak
shapes and the same dye context, and there is no allele-specific PCR bias. Real data to fit and
test against, in order of value:

1. **Mixed traces with a known answer**: heterozygous indel carriers whose alleles were
   established independently (cloning, NGS, or a parent/offspring genotype), with the indel
   size and position.
2. **Forward and reverse reads** of the same template: the reverse read puts the indel at the
   other end, so the clean stretch and the mixed stretch swap; together they solve what one
   read cannot (and are the `FR` evidence in SECONDARY_PEAKS.md).
3. A **homozygous wild-type trace** of the same amplicon (Tracy takes it as the reference, and it
   gives the per-position peak-height profile, i.e. what one allele alone looks like).
4. The amplicon reference sequence and the primers.
5. A few **heterozygous SNP-only** and **clean homozygous** traces from the same run as negatives,
   so the onset detector is tested on reads that must return "no indel".
6. Two or more indels in one read, and indels inside homopolymers or short repeats, if they exist
   in the material: these are the known failure cases of every method above.

## Known gaps

- Allele ratio far from 50 % (w = 0.3 already lowers accuracy on duplications): the minor allele's
  peaks fall towards the 0.20 floor and drop out of S_i.
- No second indel, no SNPs between the alleles after the onset beyond what chains tolerate as
  contradictions.
- The onset is searched once; a read where the clean prefix is short gives the step detector
  little to work with — the reverse read is then the fix.
