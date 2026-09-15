"""Chromatogram figure: PLOC vertical grid + HQ-relative horizontal thresholds."""
import numpy as np
import matplotlib.pyplot as plt

from .v0_2 import ALPHA, BETA

COLORS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#222222", "T": "#d62728"}


def chromatogram_axes(
    ax,
    raw,
    pred=None,
    pos0=1,
    pos1=None,
    show_ploc_grid=True,
    show_thresholds=True,
    alpha=ALPHA,
    beta=BETA,
):
    """Draw four dye channels with a *corrected* base grid.

    Vertical bars sit at PLOC (the actual scan of each called peak), not at
    evenly spaced base indices. Horizontal lines are HQ-body fractions used
    as amplitude significance thresholds (alpha = relative-tests-off,
    beta = stop).
    """
    seq = raw["seq"]
    ploc = raw["ploc"]
    channels = raw["channels"]
    n = len(seq)
    if pos1 is None:
        pos1 = n
    pos0 = max(1, pos0)
    pos1 = min(n, pos1)
    i0, i1 = pos0 - 1, pos1
    scan0 = max(0, int(ploc[i0]) - 8)
    scan1 = min(len(channels["A"]), int(ploc[i1 - 1]) + 8)
    x = np.arange(scan0, scan1)

    for b in "GATC":
        ax.plot(x, channels[b][x], color=COLORS[b], lw=0.9, label=b)

    ymax = max(channels[b][x].max() for b in "ACGT")
    ymin = min(0.0, min(channels[b][x].min() for b in "ACGT"))

    if show_ploc_grid:
        for i in range(i0, i1):
            sc = ploc[i]
            ax.axvline(sc, color="#888888", lw=0.4, alpha=0.55)
            letter = seq[i]
            color = "#222222"
            if pred is not None:
                d = pred[i + 1]
                if d.get("reason") == "low_amp_keep":
                    color = "#1a7a72"
                elif d.get("bad"):
                    color = "#aa0000"
                    ax.axvspan(sc - 3, sc + 3, color="#cc4444", alpha=0.12, lw=0)
            ax.text(sc, ymin - 0.06 * (ymax - ymin + 1), letter,
                    ha="center", va="top", fontsize=7, color=color)

    if show_thresholds and pred is not None:
        hq = next(iter(pred.values()))["hq_body"]
        ax.axhline(hq, color="#555555", lw=0.8, ls=":", label=f"HQ body ({hq:.0f} RFU)")
        ax.axhline(alpha * hq, color="#2a9d8f", lw=1.0, ls="--",
                   label=f"alpha x HQ ({alpha:.2f}, relative tests off)")
        ax.axhline(beta * hq, color="#c44", lw=1.0, ls="--",
                   label=f"beta x HQ ({beta:.2f}, stop)")

    ax.set_xlim(scan0, scan1)
    ax.set_ylim(ymin - 0.12 * (ymax - ymin + 1), ymax * 1.08)
    ax.set_xlabel("scan (PLOC grid)")
    ax.set_ylabel("RFU")
    ax.legend(loc="upper right", ncol=4, fontsize=8)
    ax.set_title(f"bases {pos0}–{pos1}  (bars = called-peak scans, not uniform time)")
    return ax


def figure_strips(raw, pred=None, bases_per_strip=80, **kwargs):
    n = len(raw["seq"])
    n_strips = int(np.ceil(n / bases_per_strip))
    fig, axes = plt.subplots(n_strips, 1, figsize=(18, 2.6 * n_strips), squeeze=False)
    for s in range(n_strips):
        p0 = s * bases_per_strip + 1
        p1 = min(n, (s + 1) * bases_per_strip)
        chromatogram_axes(axes[s][0], raw, pred, pos0=p0, pos1=p1, **kwargs)
        if s > 0:
            axes[s][0].legend().remove()
    fig.tight_layout()
    return fig
