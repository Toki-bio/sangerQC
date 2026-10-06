/**
 * Classifier in the browser — same rules as sangerqc/v0_2.py and sangerqc/v0_3.py.
 * v0.2: periodicity (local FFT), valley-depth, positive spike z, then α/β HQ-body gate.
 * v0.3 (classifyV03): same bad/ok flags, plus all-reasons, spike-over-relative label,
 * low_amp_5prime caution label, and the secondary-peak co-location test (secondary.py).
 */
(function (global) {
  const PERIOD_HALF_WIN = 60;
  const PERIOD_HALF_WIN_LOCAL = 8;
  const VALLEY_CLEAN_AFTER_BLOB = 0.42;
  const VALLEY_BLOB_BEFORE = 0.45;
  const RECOVERY_MIN_RUN = 5;
  const PERIOD_THRESH = 0.3;
  const SPIKE_THRESH = 2.0; // synced with sangerqc/v0_1.py (commit e427690)
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

  function periodicity(rec, env, halfWin) {
    const hw = halfWin || PERIOD_HALF_WIN;
    const ploc = rec.ploc;
    const diffs = [];
    for (let i = 1; i < ploc.length; i++) diffs.push(ploc[i] - ploc[i - 1]);
    const avg = mean(diffs) || 12;
    const expectedFreq = 1 / avg;
    const nScans = rec.nScans;
    return ploc.map((center) => {
      const lo = Math.max(0, center - hw);
      const hi = Math.min(nScans, center + hw);
      return periodStrength(env.subarray(lo, hi), expectedFreq);
    });
  }

  function clearPeriodSmearAfterBlob(pred, rec, env, valley, level, stopRfu) {
    const perW = periodicity(rec, env, PERIOD_HALF_WIN);
    const perL = periodicity(rec, env, PERIOD_HALF_WIN_LOCAL);
    const n = pred.length;
    let i = 1;
    while (i < n) {
      if (valley[i - 1] <= VALLEY_BLOB_BEFORE) {
        i++;
        continue;
      }
      let j = i;
      while (j < n && valley[j] <= VALLEY_CLEAN_AFTER_BLOB && level[j] >= stopRfu) j++;
      if (j - i >= RECOVERY_MIN_RUN) {
        for (let k = i; k < j; k++) {
          const p = pred[k];
          if (
            p.bad &&
            p.reason === "relative" &&
            valley[k] <= VALLEY_THRESH &&
            perW[k] < PERIOD_THRESH &&
            perL[k] >= PERIOD_THRESH
          ) {
            p.bad = false;
            p.reason = "ok";
          }
        }
      }
      i = j > i ? j : i + 1;
    }
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
    clearPeriodSmearAfterBlob(out, rec, env, valley, level, gates.stopRfu);
    return {
      pred: out,
      hq,
      island: w ? [a + 1, b + 1] : [1, n],
      noiseFloor: gates.noiseFloor,
      alphaRfu: gates.alphaRfu,
      stopRfu: gates.stopRfu,
    };
  }

  // ---- v0.3: secondary peaks (port of sangerqc/secondary.py) ----
  const SEC_RATIO_MIN = 0.2;
  const SEC_RATIO_CALL = 0.33;
  const SEC_COLOC_MAX = 0.35;
  const SAT_X_HQ = 3.0;
  const SAT_WINDOW = 10;
  const IUPAC2 = { AG: "R", CT: "Y", CG: "S", AT: "W", GT: "K", AC: "M" };
  const iupacOf = (a, b) => IUPAC2[[a, b].sort().join("")];

  /** Highest interior local maximum of ch within [L, R] (plateaus count only if they then fall); -1 if none. */
  function interiorApex(ch, L, R) {
    let best = -1, bv = -Infinity;
    let k = L + 1;
    while (k < R) {
      if (ch[k] - ch[k - 1] > 0) {
        let j = k;
        while (j < R && ch[j + 1] - ch[j] === 0) j++;
        if (j < R && ch[j + 1] - ch[j] < 0) {
          const mid = (k + j) >> 1;
          if (ch[mid] > bv) { bv = ch[mid]; best = mid; }
        }
        k = j + 1;
      } else k++;
    }
    return best;
  }

  function channelProbe(rec, i) {
    const ploc = rec.ploc, n = ploc.length, nScans = rec.nScans;
    const L = i > 0 ? Math.floor((ploc[i - 1] + ploc[i]) / 2) : Math.max(0, ploc[i] - 6);
    const R = i < n - 1 ? Math.floor((ploc[i] + ploc[i + 1]) / 2) : Math.min(nScans - 1, ploc[i] + 6);
    let sp;
    if (i > 0 && i < n - 1) sp = (ploc[i + 1] - ploc[i - 1]) / 2;
    else if (n > 1) sp = (ploc[n - 1] - ploc[0]) / (n - 1);
    else sp = 12;
    sp = Math.max(sp, 1);
    const h = {};
    for (const b of "ACGT") {
      const ch = rec.channels[b];
      let m = -Infinity;
      for (let s = L; s <= R; s++) if (ch[s] > m) m = ch[s];
      h[b] = m;
    }
    let primary = "A";
    for (const b of "CGT") if (h[b] > h[primary]) primary = b;
    const top = h[primary];
    const channels = {};
    for (const b of "ACGT") {
      const ch = rec.channels[b];
      const k = interiorApex(ch, L, R);
      if (k < 0) channels[b] = { height: h[b], ratio: top ? h[b] / top : 0, apex: false, offset: null };
      else channels[b] = { height: ch[k], ratio: top ? ch[k] / top : 0, apex: true, offset: (k - ploc[i]) / sp };
    }
    return { primary, top, channels };
  }

  function isColocated(c) {
    return c.apex && c.offset !== null && Math.abs(c.offset) <= SEC_COLOC_MAX && c.ratio >= SEC_RATIO_MIN;
  }

  function secondaryPeaks(rec, top, hq) {
    const n = rec.seq.length;
    const zone = new Uint8Array(n);
    for (let k = 0; k < n; k++) {
      if (top[k] >= SAT_X_HQ * hq) {
        for (let j = Math.max(0, k - SAT_WINDOW); j <= Math.min(n - 1, k + SAT_WINDOW); j++) zone[j] = 1;
      }
    }
    const out = [];
    for (let i = 0; i < n; i++) {
      const p = channelProbe(rec, i);
      const others = "ACGT".split("").filter((b) => b !== p.primary);
      let sec = others.reduce((a, b) => (p.channels[b].height > p.channels[a].height ? b : a));
      const ratioSlot = p.top ? p.channels[sec].height / p.top : 0;
      const coloc = others.filter((b) => isColocated(p.channels[b]));
      let cls;
      if (coloc.length) {
        sec = coloc.reduce((a, b) => (p.channels[b].ratio > p.channels[a].ratio ? b : a));
        cls = zone[i] ? "near_sat" : "colocated";
      } else cls = ratioSlot < SEC_RATIO_MIN ? "none" : "spill";
      const c = p.channels[sec];
      out.push({
        primary: p.primary, secondary: sec, ratio: c.ratio, offset: c.offset, apex: c.apex, cls,
        iupac: cls === "colocated" && c.ratio >= SEC_RATIO_CALL ? iupacOf(p.primary, sec) : null,
      });
    }
    return out;
  }

  function classifyV03(rec, alpha, beta) {
    const res = classifyV02(rec, alpha, beta);
    const level = rollingMedian(topH(rec), ROLL_WIN);
    let lead = level.length;
    for (let i = 0; i < level.length; i++) if (level[i] >= res.alphaRfu) { lead = i; break; }
    const sec = secondaryPeaks(rec, topH(rec), res.hq);
    res.pred.forEach((p, i) => {
      const reasons = [];
      if (p.local_level < res.stopRfu) reasons.push("below_stop_level");
      else if (p.local_level < res.alphaRfu) reasons.push("low_amp");
      if (p.periodicity < PERIOD_THRESH) reasons.push("periodicity");
      if (p.valley_ratio > VALLEY_THRESH) reasons.push("valley");
      if (p.spike_z > SPIKE_THRESH) reasons.push("spike");
      p.reasons = reasons;
      p.v02_reason = p.reason;
      if (p.bad && p.reason === "relative" && reasons.includes("spike")) p.reason = "spike";
      if (p.reason === "low_amp_keep" && i < lead) p.reason = "low_amp_5prime";
      p.sec = sec[i];
    });
    const hqTop = res.hq;
    res.lowSnr = hqTop < 2 * res.noiseFloor;
    res.leadLowAmpEnd = lead;
    return res;
  }

  global.classifyV02 = classifyV02;
  global.classifyV03 = classifyV03;
  global.secondaryPeaks = secondaryPeaks;
  global.amplitudeGates = amplitudeGates;
  global.noiseFloor = noiseFloor;
  global.PERIOD_THRESH = PERIOD_THRESH;
  global.SPIKE_THRESH = SPIKE_THRESH;
  global.VALLEY_THRESH = VALLEY_THRESH;
})(typeof window !== "undefined" ? window : globalThis);
