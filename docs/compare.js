/**
 * ABI PBAS vs v0.2 kept letters vs a Chromas/edited FASTA mapped onto the same PLOC.
 * v0.2 is a trim/mask on the machine calls, not a new basecaller.
 */
(function (global) {
  const RC = {
    A: "T", T: "A", G: "C", C: "G",
    M: "K", K: "M", R: "Y", Y: "R",
    W: "W", S: "S", B: "V", V: "B",
    D: "H", H: "D", N: "N",
  };

  function parseFasta(text) {
    const recs = [];
    let name = "fasta";
    let seq = "";
    const lines = String(text).split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      if (line.startsWith(">")) {
        if (seq) recs.push({ name, seq });
        name = line.slice(1).trim().split(/\s+/)[0] || "fasta";
        seq = "";
      } else seq += line.trim();
    }
    if (seq) recs.push({ name, seq });
    return recs;
  }

  function revcomp(seq) {
    let o = "";
    for (let i = seq.length - 1; i >= 0; i--) {
      const c = seq[i];
      const u = c.toUpperCase();
      const r = RC[u] || "N";
      o += c === u ? r : r.toLowerCase();
    }
    return o;
  }

  function pairScore(a, b) {
    if (a === b) return 2;
    if (a === "N" || b === "N") return 1;
    return -1;
  }

  /** Smith–Waterman; map[i] is the FASTA char placed on PBAS index i, or '-'. */
  function smithWaterman(pbas, fasta) {
    const A = pbas.toUpperCase();
    const B = fasta.toUpperCase();
    const n = A.length;
    const m = B.length;
    const GAP = -2;
    const cols = m + 1;
    const H = new Int16Array((n + 1) * cols);
    const P = new Uint8Array((n + 1) * cols);
    let best = 0, bi = 0, bj = 0;
    for (let i = 1; i <= n; i++) {
      const row = i * cols;
      const prev = (i - 1) * cols;
      for (let j = 1; j <= m; j++) {
        const d = H[prev + j - 1] + pairScore(A[i - 1], B[j - 1]);
        const u = H[prev + j] + GAP;
        const l = H[row + j - 1] + GAP;
        let v = 0, p = 0;
        if (d >= v) { v = d; p = 1; }
        if (u > v) { v = u; p = 2; }
        if (l > v) { v = l; p = 3; }
        H[row + j] = v;
        P[row + j] = p;
        if (v > best) { best = v; bi = i; bj = j; }
      }
    }
    const map = new Array(n).fill("-");
    let i = bi, j = bj;
    while (i > 0 && j > 0) {
      const p = P[i * cols + j];
      const h = H[i * cols + j];
      if (!p || h <= 0) break;
      if (p === 1) {
        map[i - 1] = fasta[j - 1];
        i--;
        j--;
      } else if (p === 2) {
        map[i - 1] = "-";
        i--;
      } else {
        j--;
      }
    }
    return { map, score: best, pbasSpan: [i, bi], fastaSpan: [j, bj] };
  }

  function alignFastaToPbas(pbas, fasta) {
    const fwd = smithWaterman(pbas, fasta);
    const rcSeq = revcomp(fasta);
    const rc = smithWaterman(pbas, rcSeq);
    if (rc.score > fwd.score * 1.05) {
      rc.rc = true;
      rc.fasta = rcSeq;
      return rc;
    }
    fwd.rc = false;
    fwd.fasta = fasta;
    return fwd;
  }

  function v02Letters(seq, pred) {
    const out = new Array(seq.length);
    for (let i = 0; i < seq.length; i++) {
      const p = pred[i];
      if (!p) {
        out[i] = "N";
        continue;
      }
      const c = seq[i] || "N";
      if (!p || p.bad) out[i] = "N";
      else if (p.reason === "low_amp_keep") out[i] = c.toLowerCase();
      else out[i] = c;
    }
    return out;
  }

  function isCall(c) {
    return c && c !== "-" && c !== "." && c.toUpperCase() !== "N";
  }

  function compareStats(seq, pred, chromasMap) {
    const v = v02Letters(seq, pred);
    let vKeep = 0, chKeep = 0, recovered = 0, lost = 0, disagree = 0, bothN = 0;
    for (let i = 0; i < seq.length; i++) {
      const vl = v[i];
      const cl = chromasMap ? chromasMap[i] : null;
      const vOk = isCall(vl);
      const cOk = isCall(cl);
      if (vOk) vKeep++;
      if (cOk) chKeep++;
      if (vOk && !cOk && chromasMap) recovered++;
      if (!vOk && cOk) lost++;
      if (vOk && cOk && vl.toUpperCase() !== cl.toUpperCase()) disagree++;
      if (!vOk && chromasMap && !cOk) bothN++;
    }
    return {
      n: seq.length,
      vKeep,
      chKeep: chromasMap ? chKeep : null,
      recovered: chromasMap ? recovered : null,
      lost: chromasMap ? lost : null,
      disagree: chromasMap ? disagree : null,
      bothN: chromasMap ? bothN : null,
    };
  }

  function nameHint(s) {
    const m = String(s).match(/(?:^|[^A-Za-z0-9])([A-H]\d)(?:[^A-Za-z0-9]|$)/);
    return m ? m[1].toUpperCase() : "";
  }

  global.parseFasta = parseFasta;
  global.revcomp = revcomp;
  global.alignFastaToPbas = alignFastaToPbas;
  global.v02Letters = v02Letters;
  global.compareStats = compareStats;
  global.nameHint = nameHint;
})(typeof window !== "undefined" ? window : globalThis);
