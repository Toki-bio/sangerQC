const DYE = { A: "#2ca02c", C: "#1f77b4", G: "#111111", T: "#d62728" };

const els = {
  file: document.getElementById("file"),
  meta: document.getElementById("meta"),
  canvas: document.getElementById("trace"),
  start: document.getElementById("start"),
  width: document.getElementById("width"),
  alpha: document.getElementById("alpha"),
  beta: document.getElementById("beta"),
  err: document.getElementById("err"),
};

let rec = null;

function peakHeights(rec) {
  const { ploc, channels } = rec;
  return ploc.map((sc) => {
    const i = Math.max(0, Math.min(sc, rec.nScans - 1));
    return Math.max(
      channels.A?.[i] ?? 0,
      channels.C?.[i] ?? 0,
      channels.G?.[i] ?? 0,
      channels.T?.[i] ?? 0
    );
  });
}

function hqBody(heights) {
  if (!heights.length) return 1;
  const mid = heights.slice(
    Math.floor(heights.length * 0.25),
    Math.ceil(heights.length * 0.75)
  );
  const s = (mid.length ? mid : heights).slice().sort((a, b) => a - b);
  return s[Math.floor(s.length / 2)] || 1;
}

function draw() {
  if (!rec) return;
  const canvas = els.canvas;
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 1100;
  const cssH = 280;
  canvas.width = Math.floor(cssW * dpr);
  canvas.height = Math.floor(cssH * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  const n = rec.seq.length;
  const start = Math.max(1, Math.min(n, Number(els.start.value) || 1));
  const width = Math.max(20, Number(els.width.value) || 80);
  const end = Math.min(n, start + width - 1);
  const i0 = start - 1;
  const i1 = end - 1;
  const scan0 = Math.max(0, rec.ploc[i0] - 8);
  const scan1 = Math.min(rec.nScans - 1, rec.ploc[i1] + 8);
  const padL = 44;
  const padR = 12;
  const padT = 16;
  const padB = 36;
  const plotW = cssW - padL - padR;
  const plotH = cssH - padT - padB;
  const sx = (scan) => padL + ((scan - scan0) / Math.max(1, scan1 - scan0)) * plotW;

  let ymax = 1;
  for (const b of "GATC") {
    const ch = rec.channels[b];
    if (!ch) continue;
    for (let s = scan0; s <= scan1; s++) ymax = Math.max(ymax, ch[s] || 0);
  }
  const sy = (rfu) => padT + (1 - rfu / ymax) * plotH;

  const heights = peakHeights(rec);
  const hq = hqBody(heights);
  const alpha = Number(els.alpha.value) || 0.35;
  const beta = Number(els.beta.value) || 0.18;

  // corrected PLOC grid
  ctx.strokeStyle = "rgba(80,80,80,0.35)";
  ctx.lineWidth = 0.7;
  for (let i = i0; i <= i1; i++) {
    const x = sx(rec.ploc[i]);
    ctx.beginPath();
    ctx.moveTo(x, padT);
    ctx.lineTo(x, padT + plotH);
    ctx.stroke();
  }

  // horizontal significance thresholds
  const lines = [
    [hq, "#555555", "HQ"],
    [alpha * hq, "#2a9d8f", "α·HQ"],
    [beta * hq, "#c44", "β·HQ"],
  ];
  ctx.setLineDash([5, 4]);
  ctx.lineWidth = 1;
  for (const [y, color] of lines) {
    if (y <= 0 || y > ymax * 1.2) continue;
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(padL, sy(y));
    ctx.lineTo(padL + plotW, sy(y));
    ctx.stroke();
  }
  ctx.setLineDash([]);

  // dye traces
  for (const b of "GATC") {
    const ch = rec.channels[b];
    if (!ch) continue;
    ctx.beginPath();
    ctx.strokeStyle = DYE[b];
    ctx.lineWidth = 1.1;
    for (let s = scan0; s <= scan1; s++) {
      const x = sx(s);
      const y = sy(ch[s] || 0);
      if (s === scan0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  // base letters under the grid
  ctx.font = "11px ui-monospace, Consolas, monospace";
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  ctx.fillStyle = "#222";
  for (let i = i0; i <= i1; i++) {
    ctx.fillText(rec.seq[i] || "", sx(rec.ploc[i]), padT + plotH + 4);
  }

  ctx.fillStyle = "#666";
  ctx.font = "11px system-ui, sans-serif";
  ctx.textAlign = "left";
  ctx.fillText(`scan ${scan0}–${scan1}`, padL, 2);
  ctx.textAlign = "right";
  ctx.fillText(
    `HQ≈${hq.toFixed(0)} RFU  α=${alpha} (${(alpha * hq).toFixed(0)})  β=${beta} (${(beta * hq).toFixed(0)})`,
    cssW - padR,
    2
  );
}

function showError(msg) {
  els.err.hidden = !msg;
  els.err.textContent = msg || "";
}

els.file.addEventListener("change", async (ev) => {
  showError("");
  rec = null;
  const f = ev.target.files && ev.target.files[0];
  if (!f) return;
  try {
    const buf = await f.arrayBuffer();
    rec = parseABIF(buf);
    if (!rec.seq.length || !rec.ploc.length) {
      throw new Error("Parsed ABIF but found no PBAS/PLOC. File may use an unsupported tag layout.");
    }
    els.meta.textContent =
      `${f.name} · ${rec.seq.length} bases · ${rec.nScans} scans · channel order ${rec.fwo}`;
    els.start.max = rec.seq.length;
    els.start.value = 1;
    draw();
  } catch (e) {
    showError(String(e.message || e));
  }
});

["start", "width", "alpha", "beta"].forEach((id) => {
  els[id].addEventListener("input", draw);
});
window.addEventListener("resize", draw);
