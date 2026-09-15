/**
 * v0.2 classifier in the browser — same rules as sangerqc/v0_2.py.
 * Periodicity (local FFT), valley-depth, positive spike z, then α/β HQ-body gate.
 */
(function (global) {
  const PERIOD_HALF_WIN = 60;
  const PERIOD_THRESH = 0.3;
  const SPIKE_THRESH = 1.5;
  const VALLEY_THRESH = 0.5;
  const ROLL_WIN = 21;
  const NOISE_MIN = 30;
  const NOISE_PCT = 0.08;
  const SNR_SPAN_ALPHA = 4;

  function noiseFloor(top) {
    if (!top.length) return NOISE_MIN;
    const sorted = top.slice().sort((a, b) => a - b);
    const idx = Math.floor(0.9 * (sorted.length - 1));
    const p90 = sorted[idx];
    return Math.max(NOISE_MIN, NOISE_PCT * p90);
  }

  function amplitudeGates(top, hq, alphaScale, betaScale) {
    const floor = noiseFloor(top);
    const stop = betaScale * hq;
    const alphaHq = alphaScale * hq;
    const headroom = Math.max(hq - floor, 0);
    const alphaSpan = floor + alphaScale * headroom;
    const snr = hq / Math.max(floor, 1);
    let alpha;
    if (snr > SNR_SPAN_ALPHA) {
      const w = Math.min(1, (snr - SNR_SPAN_ALPHA) / 6);
      alpha = (1 - w) * alphaHq + w * alphaSpan;
    } else {
      alpha = alphaHq;
    }
    alpha = Math.max(stop, alpha);
    return { alphaRfu: alpha, stopRfu: stop, noiseFloor: floor };
  }

  function envelopeOf(rec) {
    const n = rec.nScans;
    const env = new Float64Array(n);
    for (let s = 0; s < n; s++) {
      env[s] = Math.max(
        rec.channels.A?.[s] ?? 0,
        rec.channels.C?.[s] ?? 0,
        rec.channels.G?.[s] ?? 0,
        rec.channels.T?.[s] ?? 0
      );
    }
    return env;
  }

  function topH(rec) {
    return rec.ploc.map((sc) => {
      const i = Math.max(0, Math.min(sc, rec.nScans - 1));
      return Math.max(
        rec.channels.A?.[i] ?? 0,
        rec.channels.C?.[i] ?? 0,
        rec.channels.G?.[i] ?? 0,
        rec.channels.T?.[i] ?? 0
      );
    });
  }

  function mean(arr) {
    let s = 0;
    for (let i = 0; i < arr.length; i++) s += arr[i];
    return arr.length ? s / arr.length : 0;
  }

  function median(arr) {
    if (!arr.length) return 0;
    const a = arr.slice().sort((x, y) => x - y);
    const m = Math.floor(a.length / 2);
    return a.length % 2 ? a[m] : 0.5 * (a[m - 1] + a[m]);
  }

  function std(arr) {
    if (arr.length < 2) return 0;
    const mu = mean(arr);
    let v = 0;
    for (let i = 0; i < arr.length; i++) {
      const d = arr[i] - mu;
      v += d * d;
    }
    return Math.sqrt(v / arr.length);
  }

  function rollingMedian(x, win) {
    const n = x.length;
    const half = Math.floor(win / 2);
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const lo = Math.max(0, i - half);
      const hi = Math.min(n, i + half + 1);
      out[i] = median(Array.prototype.slice.call(x, lo, hi));
    }
    return out;
  }

  /** Fraction of spectral power near expectedFreq. Matches v0_1.periodicity. */
  function periodStrength(seg, expectedFreq) {
    const n = seg.length;
    if (n < 10) return 0;
    const mu = mean(seg);
    let varSum = 0;
    const x = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = seg[i] - mu;
      varSum += x[i] * x[i];
    }
    if (varSum === 0) return 0;
    const nFreq = Math.floor(n / 2) + 1;
    let total = 0;
    let band = 0;
    const tol = expectedFreq * 0.4;
    for (let k = 0; k < nFreq; k++) {
      const freq = k / n;
      let re = 0, im = 0;
      const ang = (-2 * Math.PI * k) / n;
      for (let t = 0; t < n; t++) {
        const a = ang * t;
        re += x[t] * Math.cos(a);
        im += x[t] * Math.sin(a);
      }
      const p = re * re + im * im;
      total += p;
      if (freq >= expectedFreq - tol && freq <= expectedFreq + tol) band += p;
    }
    return total > 0 ? band / total : 0;
  }

  function periodicity(rec, env) {
    const ploc = rec.ploc;
    const diffs = [];
    for (let i = 1; i < ploc.length; i++) diffs.push(ploc[i] - ploc[i - 1]);
    const avg = mean(diffs) || 12;
    const expectedFreq = 1 / avg;
    const nScans = rec.nScans;
    return ploc.map((center) => {
      const lo = Math.max(0, center - PERIOD_HALF_WIN);
      const hi = Math.min(nScans, center + PERIOD_HALF_WIN);
      return periodStrength(env.subarray(lo, hi), expectedFreq);
    });
  }

  function valleyRatio(rec, env) {
    const seq = rec.seq;
    const ploc = rec.ploc;
    const n = ploc.length;
    const gap = new Float64Array(Math.max(0, n - 1));
    for (let i = 0; i < n - 1; i++) {
      const a = (seq[i] || "").toUpperCase();
      const b = (seq[i + 1] || "").toUpperCase();
      if (a === b && "ACGT".includes(a)) {
        gap[i] = 0;
        continue;
      }
      const lo = ploc[i], hi = ploc[i + 1];
      if (hi <= lo) {
        gap[i] = 1;
        continue;
      }
      let vmin = Infinity;
      for (let s = lo; s <= hi; s++) if (env[s] < vmin) vmin = env[s];
      const peakAvg = (env[lo] + env[hi]) / 2;
      gap[i] = peakAvg > 0 ? vmin / peakAvg : 1;
    }
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const left = i > 0 ? gap[i - 1] : gap[0];
      const right = i < n - 1 ? gap[i] : gap[n - 2] || 0;
      out[i] = (left + right) / 2;
    }
    return out;
  }

  function spikeZ(top) {
    const logh = top.map((h) => Math.log1p(h));
    const sd = std(logh);
    const mu = mean(logh);
    return logh.map((v) => (sd > 0 ? (v - mu) / sd : 0));
  }

  function longestClean(v01bad) {
    let best = [0, -1, 0];
    let i = 0;
    const n = v01bad.length;
    while (i < n) {
      if (!v01bad[i]) {
        let j = i;
        while (j < n && !v01bad[j]) j++;
        if (j - i > best[2]) best = [i, j - 1, j - i];
        i = j;
      } else i++;
    }
    return best;
  }

  function classifyV02(rec, alpha, beta) {
    const env = envelopeOf(rec);
    const top = topH(rec);
    const per = periodicity(rec, env);
    const valley = valleyRatio(rec, env);
    const spike = spikeZ(top);
    const n = rec.seq.length;
    const v01bad = new Array(n);
    for (let i = 0; i < n; i++) {
      v01bad[i] = per[i] < PERIOD_THRESH || Math.abs(spike[i]) > SPIKE_THRESH || valley[i] > VALLEY_THRESH;
    }
    const [a, b, w] = longestClean(v01bad);
    const hq = w ? median(top.slice(a, b + 1)) : median(top);
    const gates = amplitudeGates(top, hq, alpha, beta);
    const level = rollingMedian(top, ROLL_WIN);
    const out = [];
    for (let i = 0; i < n; i++) {
      const relBad = per[i] < PERIOD_THRESH || valley[i] > VALLEY_THRESH;
      const spikeBad = spike[i] > SPIKE_THRESH;
      const lv = level[i];
      let bad, reason;
      if (lv < gates.stopRfu) {
        bad = true;
        reason = "below_stop_level";
      } else if (lv < gates.alphaRfu) {
        bad = spikeBad;
        reason = spikeBad ? "spike" : "low_amp_keep";
      } else {
        bad = relBad || spikeBad;
        reason = relBad ? "relative" : spikeBad ? "spike" : "ok";
      }
      out.push({
        pos: i + 1,
        bad,
        reason,
        periodicity: per[i],
        spike_z: spike[i],
        valley_ratio: valley[i],
        top_h: top[i],
        local_level: lv,
        hq_body: hq,
      });
    }
    return {
      pred: out,
      hq,
      island: w ? [a + 1, b + 1] : [1, n],
      noiseFloor: gates.noiseFloor,
      alphaRfu: gates.alphaRfu,
      stopRfu: gates.stopRfu,
    };
  }

  global.classifyV02 = classifyV02;
  global.amplitudeGates = amplitudeGates;
  global.noiseFloor = noiseFloor;
  global.PERIOD_THRESH = PERIOD_THRESH;
  global.SPIKE_THRESH = SPIKE_THRESH;
  global.VALLEY_THRESH = VALLEY_THRESH;
})(typeof window !== "undefined" ? window : globalThis);
