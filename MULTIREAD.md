# Several reads of the same template (`sangerqc/multiread.py`, `multicli.py`)

When one template is read more than once (forward + reverse, a re-run, a second primer), the reads are
combined into one genotype per reference position instead of being compared by eye.

```
python -m sangerqc.multicli --reference REF.fasta --out DIR \
    --template S01=S01_fwd.ab1,S01_rev.ab1 --template S02=S02_a.ab1,S02_b.ab1 [--snp 191:C/G] [--distinct-primers]
python -m sangerqc.multicli --reference REF.fasta --out DIR --sheet reads.tsv     # columns: template, path, [relation]
```

Writes per template a consensus FASTA, a table of positions with a heterozygous call, low genotype
quality or a conflict, and `summary.json`.

## Model

Each read gives a likelihood over the 10 unordered base-pair genotypes at a position (reverse reads are
complemented onto the reference strand):

- homozygous X: 1 if X is the read's primary base, else 0.001;
- heterozygous XY: 10^LBF (the read's own prior-free Bayes factor for two bases, BAYES.md) × an identity
  term (the pair is the primary plus the strongest co-located second: 1; the primary alone when no second
  was seen: 1/3; otherwise small).

Reads of one template combine as a weighted product, log L(g) = Σ w · log l(g), times a prior: flat
(P(het) = the reads' own estimated fraction of two-base calls) or, for a known biallelic SNP, Hardy–Weinberg
over the two alleles. Reads flagged bad at that position are skipped. A **conflict** is two reads each
individually confident (> 0.9) in different genotypes; it is reported, not settled by majority.

Weights (how independent the reads are) are assumptions:

| relation | weight | meaning |
|---|---|---|
| FR | 1.0 | opposite strands: different sequence context |
| DP | 1.0 | same direction, different primers |
| SAME / RR | 0.5 | same primer again (default for a same-strand read): shares template and context artefacts |
| DM | 0.7 | same primer, other machine/run |

Only **FR is tested** (below). SAME, RR, DP and DM have no real data here and are placeholders.

## Test: 23 samples × 2 SNPs × forward + reverse (human DAO/AOC1, nested PCR, SeqStudio)

Files: C:\work\Abdushukur\Plate_20200318_Allergology and Plate_20200402_172007_Allergology (96 wells).
Truth: genotypes in the study's report, called by eye in Lasergene/SnapGene (manual grade M1, not an
independent method). `validation/eval_multiread.py` reads only sample numbers and genotypes from the report.

**Layout was not given and was inferred from the reads.** Each plate holds one SNP: three columns of
forward reads, three of reverse; 23 samples fill the first 8 × 3 wells column by column. The order of the
three reverse columns was chosen by forward/reverse agreement alone (no report): plate 1 columns 4, 6, 5
(20/20 pairs agree, against 17/20 for the next order), plate 2 columns 7, 8, 9 (22/23 vs 17/22). I first
noticed the plate-1 swap by comparing double-peak wells with the report, then confirmed it with the
truth-free search; the 23/23 match of the forward block to the report's sample numbers supports the rest.
Plate 1 = rs10156191 (C/T), plate 2 = rs1049793 (G/C), identified from the number of double-peak
samples (6 and 10, as in the report).

Zygosity (heterozygous vs homozygous), samples:

| | rs10156191 (6 het of 23) | rs1049793 (10 het of 23) |
|---|---|---|
| basecaller's own IUPAC code, per read | 11/11 het reads, 32/32 hom | 20/20, 26/26 |
| sangerQC rule (second ≥ 0.33), per read | 11/11, 32/32 | 19/20, 26/26 |
| Bayes, single read p > 0.5 | 11/11, 32/32 | 19/20, 26/26 |
| **joint forward + reverse**, flat prior | 6/6, 17/17 | 10/10, 13/13 |
| joint, biallelic prior | 6/6, 17/17 | 10/10, 13/13 |

Log-loss of the joint p(het) against the report: rs1049793 flat 0.026, biallelic 0.001, weights 0.5 → 0.092,
forward only 0.029, reverse only 0.343; rs10156191 0.000 / 0.003 / 0.008. Full weight for the opposite-strand
pair fits better than half weight, as expected for independent contexts.

**What this does and does not show.** The plates are easy: the instrument's own IUPAC codes already get every
read right, so the Bayesian layer matches the caller; it does not beat it. What the joint call adds is a
calibrated probability and tolerance for one read missing the minor allele. The one disagreeing pair
(rs1049793 sample 2, a true heterozygote): forward p = 0.99, reverse p = 0.00. With a 0.1 % miss rate
(fitted on synthetic data, where the second peak is never missing) the joint call said homozygous;
the miss rate floor is now 2 % (default 5 %), motivated by this very observation, so this sample is no
longer independent evidence.

**False-positive rate on real data** (all 23 plate-2 templates, `multicli`, SNP 191 declared): 5 719
sample-positions covered; heterozygous calls outside the two SNPs (columns 148 and 191): 26 (0.45 %), 13 of
them inside reference positions 40–280, the rest near read ends; 25 conflict positions, clustered around
reference 242–254 and 282. The call at the declared SNP matches the report in 23/23 samples. Reference
column 148 (A/G, 10 samples heterozygous, forward and reverse agree 22/22) is a second real polymorphism
not in the report; probably the study's third SNP rs1049742, not checked.

## Bug found and fixed on this data (all earlier results re-run)

The primary channel at a call was the channel with the highest value anywhere in the slot. A small G after a
large C had the C's tail at the slot edge taller than its own apex, so the primary became C and the true G
looked like a second base. Primary is now the tallest channel with a real apex on the call (Python and JS;
parity re-checked on 9 reads). Synthetic benchmarks are unchanged to the digit (they reuse a trace's own
clean peaks); on the cyclops reads it removes a spurious strong site (reference 269: LBF +2.0 → −0.5).

## Not tested, stated plainly

- Same-primer repeats, different primers, other machines: no real data; weights untested.
- Only heterozygous SNPs; no indel carriers among these plates.
- The report's genotypes are manual calls from the same chromatograms; where the human and the model agree
  they may share the same blind spots.
- `--snp` takes alleles on the reference strand; on plate 2 the data fixes it (C/G), the report's "G/C"
  is strand-agnostic.
