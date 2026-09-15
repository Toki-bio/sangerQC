# Roadmap

sangerQC is not a fitted ML model on a few files. It is a **calibration instrument**: point at runs, name the class, keep only rules that survive a new file.

## 1. Graphical output

**Status: hosted viewer has PLOC grid, HQ/α/β lines, horizontal scroll, auto-Y on the visible window, Y and px/scan sliders, α/β sliders, multiple open files.**

- Vertical bars on the **PLOC grid** (scan of each called peak), not a uniform base/time axis. That is the “corrected position grid”.
- Horizontal lines for **significance thresholds across the four dye signals**, currently this read’s own HQ-body amplitude: HQ, α·HQ (relative tests off), β·HQ (stop).
- Still to do:
  - Overlay classifier reason codes (periodicity / valley / spike / low-amp / stop) on the same grid.
  - Per-channel secondary/primary ratio threshold (0.33) as an optional horizontal, in local RFU, not a global line that mis-scales.
  - Click a bar → numbers for that position (heights, ratio, periodicity, valley).
  - Export the visible window as SVG/PNG for the calibration loop.

## 2. Two layers in the learning process

Do not dump every mismatch into one `bad` flag. Separate:

### Layer A — stochastic

Irreducible, not a function of scan index.

- Baseline wobble, dye noise, polymerase hitchhiking, isolated single-base ratio flips.
- Characterized as **distributions** (noise floor vs peak height on this machine/dye set), not as a sinusoid.
- Isolated flags inside an otherwise clean stretch are this layer. They are not a “region”.

### Layer B — structured, writable as a process

Measurable, predictable, describable mathematically or as a short program. These get **named classes** and, where possible, a generating model (slow envelope, Fourier peak at 1/mean spacing, shared-apex / no-valley, …).

Known members (add, do not collapse):

| Class | What it is | Handle |
|---|---|---|
| `bump_start` | Injection / dye-blob / primer-dimer at 5′; amplitude often huge, timing not gene | span/start, not a peak-ratio disaster |
| `slow_decay_readable` | Long HQ 3′ tail: each letter still a peak, vertical resolution falls, **relative** channels wobble because the level is low | α·HQ gate (v0.2); not periodicity |
| `stop_level` | Same decay, below β·HQ: no longer data | trim |
| `shared_hump` / `g_raise` | One broad apex sliced into several fake calls; G-channel elevation | valley-depth; not FWHM |
| `merged_blob` | Unresolved disaster, often BLAST-wrong even when the letter is unambiguous | valley + BLAST density |
| `oversaturation` | 10–40× HQ, spectral pull-up | positive amplitude z |
| `homopolymer` | Real AA/GG/TTT does not dip | valley exemption |
| `periodicity_loss` | Timing gone (true mistimed peaks), distinct from AM leakage of a slow envelope | FFT power near 1/spacing, **after** amplitude is still in the HQ regime |
| `primer_not_locus` | Sequence before the gene; BLAST/span, not quality | span stage |

Layer B is where sinusoid / Fourier / slow-envelope models belong. Layer A is what remains after Layer B is subtracted. Mixing them produced the H4 failure: periodicity (Layer B, timing) firing on slow decay (Layer B, amplitude) plus stochastic relative wobble (Layer A).

A new labeled zone must be tagged with a class from this table (or a new named class), never only `0/1`.

## 3. Ground truth

A label is not automatically truth. Record **how** it was established.

### 3a. Manual (chromatogram)

The human looks at the trace — not a screenshot verdict from the model — and names the class and the span. This is Grade **M1**. A second pass on the same span, or a second person, is **M2**.

Current Artemia C1 / G4 / B3 / H4 zones are M1 unless noted.

### 3b. Concordance (independent molecules / primers / machines)

If two observations of the **same biological fragment** agree, that agreement is ground truth for the letter, independent of trace aesthetics.

| Code | Evidence |
|---|---|
| `FR` | Forward and reverse reads overlap and agree |
| `RR` | Two or more re-reads from the **same** primer agree |
| `DP` | Two reads same direction, **different** primers, agree |
| `DM` | Same primer, **different** machines / runs, agree |
| `REF` | BLAST / close reference agrees (supporting, not sufficient alone) |

Grade **C1** = one concordance type. **C2** = two independent types (e.g. FR + DM). Disagreement between replicates is also truth: that position is unresolved, not a silent majority vote.

### 3c. How a position enters the training set

```
position → {file, pos, called, class, grade, evidence[], annotator, date}
```

- `class` from the Layer B table, or `stochastic`, or `clean`.
- `grade` = `M1` / `M2` / `C1` / `C2`.
- Grey / “not totally bad” stays unlabeled. Do not force a binary.
- Concordance upgrades a manual call; it does not replace looking at a class (a blob that happens to match reverse is still a blob for the **quality** model; it may still be a true letter for the **FASTA**).

Keep those two jobs separate: **call reliability** (trace class) vs **letter identity** (concordance / reference). The converter uses both; the classifier is only the first.

Schema: [`schema/ground_truth.schema.json`](schema/ground_truth.schema.json).

## 4. Near-term calibration (this project)

- [x] H4 572–705 named as `slow_decay_readable` → v0.2 α gate
- [ ] H4 431–571: same class or not? (still flagged at α=0.35)
- [ ] H4 β stop vs old FASTA cut at 705
- [ ] **B3** with v0.2 as the pointer (open 489–511 was `slow_decay_readable` at ~51% of B3 HQ — milder than H4; α=0.35 does not clear it)
- [ ] Label A2, A3; do not add features until those exist
- [ ] Rebuild FASTA from one converter (v0.2 + span + letter policy), not two pipelines
