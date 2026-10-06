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

2026-10-06, 84 synthetic cases from three clean Artemia reads (C1 calls 67–271, A2 66–285,
A3 200–410), indel at 40 % of the zone:

| read | type | w (allele A) | size correct (ours) | allele accuracy mean / min (ours) | size correct (Tracy) |
|---|---|---|---|---|---|
| C1 | deletion | 0.5 / 0.3 | 7/7, 7/7 | 0.981 / 0.966, 0.985 / 0.958 | 7/7, 7/7 |
| C1 | duplication | 0.5 / 0.3 | 7/7, 7/7 | 0.982 / 0.966, 0.948 / 0.912 | 7/7, 7/7 |
| A2 | deletion | 0.5 / 0.3 | 7/7, 7/7 | 0.979 / 0.938, 0.972 / 0.930 | 7/7, 7/7 |
| A2 | duplication | 0.5 / 0.3 | 7/7, 7/7 | 0.981 / 0.938, 0.930 / 0.829 | 7/7, 7/7 |
| A3 | deletion | 0.5 / 0.3 | 7/7, 7/7 | 0.924 / 0.885, 0.942 / 0.858 | 7/7, 7/7 |
| A3 | duplication | 0.5 / 0.3 | 7/7, 7/7 | 0.926 / 0.868, 0.918 / 0.888 | 7/7, 7/7 |

Both find the indel size in 84/84: the synthetic set is too easy to separate the methods on size.
Allele accuracy (ours) is per call from the onset to the end of the clean zone, against the
known A and B; truth A is the source read's own calls, so its residual miscalls cap the score
(A3 is the noisiest source). Tracy's alleles come out in its own re-basecalled coordinates, so
only its size call is scored here. The onset is "exact" in 32/84; the rest sit a few calls away
inside repeats where the indel position is not unique, which does not affect the alleles.

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
