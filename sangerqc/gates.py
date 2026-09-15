"""Per-file amplitude gates: HQ body + adaptive noise floor."""
import numpy as np

NOISE_MIN = 30.0
NOISE_PCT = 0.08
# Above this HQ/noise ratio, alpha uses headroom above floor (high-yield reads).
SNR_SPAN_ALPHA = 4.0


def noise_floor(top_h):
    """Same rule as Artemia analyze_ab1: max(30 RFU, 8% of p90 peak height)."""
    th = np.asarray(top_h, dtype=float)
    if th.size == 0:
        return NOISE_MIN
    return float(max(NOISE_MIN, NOISE_PCT * np.percentile(th, 90)))


def amplitude_gates(top_h, hq, alpha_scale=0.35, beta_scale=0.18):
    """
    Return (alpha_rfu, stop_rfu, noise_floor).

    Stop = beta_scale * HQ (unchanged from v0.2). Noise floor is computed per
    file and drives alpha on high-SNR reads.

    Alpha scale: below stop, relative tests are off. Between stop and alpha,
    low_amp_keep. Uses HQ fraction on low-SNR reads; blends toward
    floor + alpha*(HQ-floor) when HQ >> floor so the low-amp band tracks
    intensity headroom, not a fixed fraction of a 800+ RFU HQ body alone.
    """
    floor = noise_floor(top_h)
    hq = float(hq)
    # Stop line: relative to HQ (H4-calibrated). Floor is shown and used for alpha headroom only.
    stop = beta_scale * hq
    alpha_hq = alpha_scale * hq
    headroom = max(hq - floor, 0.0)
    alpha_span = floor + alpha_scale * headroom
    snr = hq / max(floor, 1.0)
    if snr > SNR_SPAN_ALPHA:
        w = min(1.0, (snr - SNR_SPAN_ALPHA) / 6.0)
        alpha = (1.0 - w) * alpha_hq + w * alpha_span
    else:
        alpha = alpha_hq
    alpha = max(stop, alpha)
    return float(alpha), float(stop), floor
