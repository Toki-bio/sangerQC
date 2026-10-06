# sangerQC

Quality-aware conversion of Sanger `.ab1` chromatograms to sequence, from the **raw four-channel traces**, not from the basecaller Phred score.

The GitHub repo page is source code only — there is no upload control there.

Interactive viewer:

**https://toki-bio.github.io/sangerQC/**

(Old link `/docs/index.html` redirects to the same viewer.)

Horizontal scroll along the full read, auto vertical zoom on the visible window, Y / px-per-scan / α / β sliders, several `.ab1` files at once. The file does not leave the machine.

Vertical bars are the called-peak scan positions (`PLOC`) — a corrected position grid, not evenly spaced time. Horizontal dashed lines are this read’s own HQ-body amplitude thresholds (α = relative tests off, β = stop).

## Python classifier (local)

```
pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

v0.1: periodicity (Phred-style local FFT) + oversaturation spike + valley-depth (shared-hump / G-raise), with a homopolymer exemption.

v0.2: same, plus a **level gate** for the common long high-quality read whose 3′ end slowly loses amplitude. Gates are **RFU lines** per file: noise floor `max(30, 8% of p90 peak height)`, stop `max(floor, β×HQ)`, α from HQ and (on high-SNR reads) headroom above the floor. Sliders still set β and α scale (defaults 0.18 / 0.35). Relative metrics are skipped below α; spike uses positive z only.

v0.3 (current; the hosted viewer runs it too): **same bad/ok flags as v0.2**, plus

- **secondary peaks**: for every call, does the second channel have its *own apex on the called peak* (`colocated`, a real second signal) or is it the tail of a neighbouring peak (`spill`)? Only co-located seconds ≥ 0.33 become IUPAC; a machine IUPAC whose second channel is a tail becomes the primary base. Co-located seconds within ±10 calls of a peak ≥ 3× HQ are `near_sat` (pull-up) and keep the machine letter. See [SECONDARY_PEAKS.md](SECONDARY_PEAKS.md).
- **all reasons** that fired, and `spike` instead of `relative` when both fire.
- **`low_amp_5prime`**: a leading low-amplitude run, shown lowercase. It can be start-up junk *or* real weak gene (both seen), so only a reference check trims it.
- **`low_snr`** file flag when HQ < 2 × noise floor (the gates stop meaning anything there).
- **primers**: give them and the primer site at the 3′ end (and the non-templated +A after it) is cut.
- **concordance**: several reads on a shared reference; a second peak reproduced in ≥ 2 reads (co-located, ≥ 0.20) is coded IUPAC.

```
python -m sangerqc.cli READ1.ab1 READ2.ab1 --out out/ --forward AAACTTAAAGGAATTGACG --reverse CTAAGGGCATCACAGACC --reference ref.fasta
```

writes `sangerqc_v0_3.fasta`, `<read>_decisions.tsv` (every non-trivial position with all four channel heights and apex offsets), `reproduced_sites.tsv` and `summary.json`.
With `--reference`, the ends are also trimmed where a 10-base window disagrees with the reference > 15 % (the Artemia BLAST bar).

Checks behind v0.3 (`validation/`): labelled C1/G4/B3 zones unchanged (accuracy 0.933 / 0.898 / 0.786 as v0.2); bad/ok flags identical to v0.2 on all six Artemia and two Thermocyclops reads; Python and browser JS agree on every position of those eight reads (`validation/parity.py`).
The co-location test was calibrated on two reads only; the false-positive control (a non-called base passing at ordinary positions) was 1.0 % and 2.3 % there. Treat it as a candidate rule until it survives more labelled files.

Heterozygous indels (prototype, not in the CLI yet): `sangerqc/hetindel.py` (`decode_read`) splits a mixed trace into its two alleles using the co-located peak sets, with or without a reference; `sangerqc/synth.py` makes synthetic het-indel .ab1 files from clean reads; benchmark against Tracy in [HET_INDEL.md](HET_INDEL.md).

Never trim from Phred alone. Never issue a verdict by eyeballing a rendered plot.

## What this repo is for

The conversion algorithm is being fit by a calibration loop: the software **points** at flagged runs with numbers; a human looks at the chromatogram and names the class; a rule is kept only if it survives a file it was not designed on. See [ROADMAP.md](ROADMAP.md) for the two-layer model (stochastic vs structured artifacts) and the ground-truth protocol (manual call vs forward/reverse / re-read concordance).

Unpublished chromatograms are not stored here. Bring your own `.ab1`.
