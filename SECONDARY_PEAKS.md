# Recording secondary peaks

A double peak is three different things that are easy to collapse into one IUPAC letter:

1. **a measurement** — what the four channels do at one called position of one read;
2. **a site** — the same reference coordinate seen in several reads, which either agree or not;
3. **an interpretation** — what the second signal is biologically (or technically).

Record them separately, in that order, and never let a later layer overwrite an earlier one.
The IUPAC letter in a FASTA is only a summary of layer 2; it carries none of the evidence.

## Layer 1 — per-position measurement (machine-generated)

One row for every position where any non-called channel reaches `ratio >= 0.20`.
`sangerqc.cli` writes this as `<read>_decisions.tsv`; `secondary.py` holds the full numbers.

| field | meaning |
|---|---|
| `read`, `file_md5` | which trace; md5 so a renamed file is still the same evidence |
| `instrument`, `run_date`, `primer`, `strand` | from the ab1 header / your sample sheet |
| `pos`, `scan` | 1-based call index and its PLOC scan |
| `ref_id`, `ref_pos` | shared coordinate (filled once the read is placed on a reference) |
| `caller`, `phred` | the machine's letter and score, kept verbatim |
| `h_A h_C h_G h_T` | each channel's maximum inside the call's slot (RFU) |
| `apex_off_A … apex_off_T` | offset of each channel's own interior apex from the called peak, in peak spacings; empty if no apex |
| `primary`, `secondary`, `ratio` | strongest channel, strongest second channel, their height ratio |
| `sec_class` | `none` / `spill` (neighbour's tail: no apex or offset > 0.35) / `colocated` / `near_sat` (co-located inside ±10 calls of a peak ≥ 3 × HQ: pull-up) |
| `local_level`, `hq`, `low_snr` | amplitude context: a ratio at 15 RFU is not the same evidence as at 400 RFU |
| `sangerqc_version`, `params` | `ratio_min`, `ratio_call`, `colocate_max`; a class means nothing without them |

Keep the **ratio**, not only the class. Peak-height ratios in Sanger are only
semi-quantitative (dye- and context-dependent incorporation), so a ratio is not an
allele or copy proportion unless it has been calibrated with known template mixes.
Record it anyway; it is what a later calibration needs.

## Layer 2 — site record (one per reference coordinate across reads)

Schema: [`schema/secondary_site.schema.json`](schema/secondary_site.schema.json).
`sangerqc.cli --reference` writes `reproduced_sites.tsv` as the starting point.

| field | meaning |
|---|---|
| `ref_id`, `ref_pos` | coordinate |
| `states` | the two (or more) bases seen, e.g. `["A","G"]` |
| `reads[]` | per read: `read`, `pos`, `dominant`, `ratio`, `sec_class` — including reads where the site is **clean** (absence is evidence too) |
| `n_reads_with`, `n_reads_covering` | e.g. 2 of 2 |
| `evidence[]` | how the agreeing reads are related (codes below) |
| `control_pass_rate` | per read: how often a non-called base passes the same co-location test at ordinary positions (the false-positive background; 1–2 % on the 2026-10 Thermocyclops reads) |
| `iupac` | the summary letter written into the FASTA |

Evidence codes (extends the ground-truth codes in ROADMAP §3b):

| code | what agrees | strength |
|---|---|---|
| `FR` | forward and reverse read of the same template | strong: different sequence context, different chemistry direction |
| `RR` | re-reads, same primer, same template | weak: shares every context artifact |
| `DP` | different primers, same direction, same template | medium |
| `DM` | same primer, different machine/run | weak–medium |
| `DS` | **different specimens, same primer** | medium for "real signal", but a sequence-context artifact (compression, hairpin) reproduces too |
| `CLONE` | cloned copies show both states | strong for "two templates exist in the extract" |
| `NGS` | amplicon / genome reads show both states | strong, and gives proportions |
| `REF` | a database sequence carries the second state | supporting only |

## Layer 3 — interpretation

One label per site, with the evidence that justifies it, and **`unresolved` is a valid answer**.

| class | what would support it |
|---|---|
| `heterozygous_snp` | single-copy nuclear locus; FR; ratio near a stable value across specimens |
| `intragenomic_variant` | multi-copy locus (rDNA, ITS); reproduced in several specimens (DS) with varying ratios; clean single-copy markers (e.g. COI) in the same extracts |
| `heteroplasmy` | mitochondrial locus; FR; ratio differs between tissues/individuals |
| `mixed_template` | second states match another taxon at most of its diagnostic sites, including sites outside the double-peak block; single-copy markers also double |
| `het_indel_phase_shift` | double peaks begin at one point and continue to the end of the read; the second track is the first shifted by k bases (see HET_INDEL.md) |
| `artifact_pullup`, `artifact_dye_blob`, `artifact_compression` | saturated neighbour; broad shared hump; irregular peak spacing at a GC hairpin |
| `unresolved` | evidence does not separate the classes above |

## What goes into the deliverables

- **FASTA / GenBank**: IUPAC letters at layer-2 sites (GenBank accepts IUPAC codes).
  Describe them in a `/note` on the feature, e.g. *"IUPAC codes at 15 positions reflect
  double peaks reproduced in two specimens; likely intragenomic rDNA variation"*.
- **Paper**: a table of sites (`ref_pos`, states, ratio per specimen, evidence) and one
  chromatogram figure of the block. Keep the `.ab1` files with the data so the claim can be re-checked.
- **Never** overwrite the caller's letters or the per-position numbers; the decisions table is the audit trail.

## Worked example (Thermocyclops taihokuensis 18S, Sabirsay, 2026-10-06)

Two specimens, primer d5, reads F2 and G2. 16 sites in LC851021 positions 196–269 have
co-located second peaks in both reads (`DS`, 2/2), control pass rate 1.0 % / 2.3 %.
Second states at 10 sites equal the state most other cyclopids carry; 5 match no reference.
Interpretation: `unresolved` between `intragenomic_variant` and `mixed_template`;
an r6 (reverse) read would add `FR`, and COI from the same extracts would separate the two.
