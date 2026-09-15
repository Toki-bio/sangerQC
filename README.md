# sangerQC

Quality-aware conversion of Sanger `.ab1` chromatograms to sequence, from the **raw four-channel traces**, not from the basecaller Phred score.

Interactive viewer (upload an `.ab1` in the browser; the file does not leave the machine):

**https://toki-bio.github.io/sangerQC/**

Vertical bars are the called-peak scan positions (`PLOC`) — a corrected position grid, not evenly spaced time. Horizontal dashed lines are this read’s own HQ-body amplitude thresholds (α = relative tests off, β = stop).

## Python classifier (local)

```
pip install -r requirements.txt
streamlit run app/streamlit_app.py
```

v0.1: periodicity (Phred-style local FFT) + oversaturation spike + valley-depth (shared-hump / G-raise), with a homopolymer exemption.

v0.2: same, plus a **level gate** for the common long high-quality read whose 3′ end slowly loses amplitude. Relative metrics are skipped below α × HQ-body; the read is treated as ended below β × HQ-body. Spike uses positive z only (decay is not oversaturation). Working values after H4 calibration: α = 0.35, β = 0.18.

Never trim from Phred alone. Never issue a verdict by eyeballing a rendered plot.

## What this repo is for

The conversion algorithm is being fit by a calibration loop: the software **points** at flagged runs with numbers; a human looks at the chromatogram and names the class; a rule is kept only if it survives a file it was not designed on. See [ROADMAP.md](ROADMAP.md) for the two-layer model (stochastic vs structured artifacts) and the ground-truth protocol (manual call vs forward/reverse / re-read concordance).

Unpublished chromatograms are not stored here. Bring your own `.ab1`.
