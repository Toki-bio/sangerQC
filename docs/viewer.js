const DYE = { A: "#2ca02c", C: "#1f77b4", G: "#111111", T: "#d62728" };
const TRACE_H = 236;
const ROW_H = 14;
const PAD = { l: 56, r: 8, t: 22, b: 8 };
const REASON_FILL = {
  relative: "rgba(204,68,68,0.20)",
  spike: "rgba(230,140,20,0.22)",
  low_amp_keep: "rgba(42,157,143,0.20)",
  below_stop_level: "rgba(90,40,90,0.24)",
};
const CLASS_FILL = {
  clean: "rgba(34,120,80,0.35)",
  stochastic: "rgba(120,120,120,0.35)",
  bump_start: "rgba(160,80,20,0.40)",
  slow_decay_readable: "rgba(42,157,143,0.45)",
  stop_level: "rgba(90,40,90,0.45)",
  shared_hump: "rgba(180,60,140,0.40)",
  g_raise: "rgba(40,140,40,0.40)",
  merged_blob: "rgba(180,30,30,0.45)",
  oversaturation: "rgba(230,140,20,0.45)",
  homopolymer: "rgba(80,80,160,0.35)",
  periodicity_loss: "rgba(180,80,80,0.40)",
  primer_not_locus: "rgba(100,100,40,0.40)",
  unresolved_replicate_conflict: "rgba(0,0,0,0.35)",
};
const CLASSES = Object.keys(CLASS_FILL);

const els = {
  file: document.getElementById("file"),
  drop: document.getElementById("drop"),
  panels: document.getElementById("panels"),
  err: document.getElementById("err"),
  px: document.getElementById("px"),
  pxOut: document.getElementById("pxOut"),
  yZoom: document.getElementById("yZoom"),
  yOut: document.getElementById("yOut"),
  autoY: document.getElementById("autoY"),
  overlay: document.getElementById("overlay"),
  lock: document.getElementById("lock"),
  alpha: document.getElementById("alpha"),
  alphaOut: document.getElementById("alphaOut"),
  beta: document.getElementById("beta"),
  betaOut: document.getElementById("betaOut"),
  exportBtn: document.getElementById("exportLabels"),
  dlg: document.getElementById("labelDlg"),
  dlgRange: document.getElementById("dlgRange"),
  dlgClass: document.getElementById("dlgClass"),
  dlgOk: document.getElementById("dlgOk"),
};

const files = [];
const pendingFastas = [];
let lockBusy = false;
let pendingLabel = null;

function letterRowCount(item) {
  return item && item.chromas ? 3 : 2;
}

function cssH(item) {
  return PAD.t + TRACE_H + PAD.b + letterRowCount(item) * ROW_H + 6;
}

function showError(msg) {
  els.err.hidden = !msg;
  els.err.textContent = msg || "";
}

function fileMax(rec) {
  let m = 1;
  for (const b of "GATC") {
    const ch = rec.channels[b];
    if (!ch) continue;
    for (let i = 0; i < ch.length; i++) if (ch[i] > m) m = ch[i];
  }
  return m;
}

function visibleScans(viewport, rec, px) {
  const sl = viewport.scrollLeft;
  const w = viewport.clientWidth;
  const scan0 = Math.max(0, Math.floor((sl - PAD.l) / px) - 2);
  const scan1 = Math.min(rec.nScans - 1, Math.ceil((sl + w - PAD.l) / px) + 2);
  return [scan0, scan1];
}

function visibleMax(rec, scan0, scan1) {
  let m = 1;
  for (const b of "GATC") {
    const ch = rec.channels[b];
    if (!ch) continue;
    for (let s = scan0; s <= scan1; s++) if (ch[s] > m) m = ch[s];
  }
  return m;
}

function niceTicks(ymax) {
  if (!(ymax > 0)) return [0];
  const raw = ymax / 4;
  const exp = Math.pow(10, Math.floor(Math.log10(raw)));
  const f = raw / exp;
  const nf = f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10;
  const step = nf * exp;
  const ticks = [];
  for (let v = 0; v <= ymax * 1.001; v += step) ticks.push(v);
  if (ticks[ticks.length - 1] < ymax * 0.92) ticks.push(ymax);
  return ticks;
}

function controls() {
  return {
    px: Number(els.px.value) || 1,
    yZoom: Number(els.yZoom.value) || 1,
    autoY: els.autoY.checked,
    overlay: els.overlay.checked,
    lock: els.lock.checked,
    alpha: Number(els.alpha.value) || 0.35,
    beta: Number(els.beta.value) || 0.18,
  };
}

function syncOutputs() {
  const c = controls();
  els.pxOut.textContent = c.px.toFixed(1);
  els.yOut.textContent = c.yZoom.toFixed(2) + "×";
  els.alphaOut.textContent = c.alpha.toFixed(2);
  els.betaOut.textContent = c.beta.toFixed(2);
}

function scanToBase(rec, scan) {
  const ploc = rec.ploc;
  let lo = 0, hi = ploc.length - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (ploc[mid] < scan) lo = mid + 1;
    else hi = mid;
  }
  if (lo > 0 && Math.abs(ploc[lo - 1] - scan) < Math.abs(ploc[lo] - scan)) return lo - 1;
  return lo;
}

function runClassify(item) {
  const c = controls();
  const { pred, hq, island } = classifyV02(item.rec, c.alpha, c.beta);
  item.pred = pred;
  item.hq = hq;
  item.island = island;
  item.v02 = v02Letters(item.rec.seq, pred);
  refreshStats(item);
}

function refreshStats(item) {
  item.stats = compareStats(item.rec.seq, item.pred, item.chromas ? item.chromas.map : null);
  if (!item.panel) return;
  const meta = item.panel.querySelector(".meta");
  if (!meta) return;
  const rec = item.rec;
  const st = item.stats;
  let t =
    `${rec.seq.length} ABI · v0.2 keep ${st.vKeep} · ` +
    `HQ ${item.hq.toFixed(0)} island ${item.island[0]}–${item.island[1]}`;
  if (st.chKeep != null) {
    t +=
      ` · Chromas keep ${st.chKeep}` +
      ` · recovered ${st.recovered} · lost ${st.lost}` +
      (st.disagree ? ` · letter≠ ${st.disagree}` : "");
    if (item.chromas.rc) t += " · FASTA was RC";
  } else {
    t += " · drop Chromas FASTA to compare";
  }
  meta.textContent = t;
}

function drawOne(item) {
  const rec = item.rec;
  const c = controls();
  const canvas = item.canvas;
  const viewport = item.viewport;
  if (!canvas || !viewport) return;
  const dpr = window.devicePixelRatio || 1;
  const height = cssH(item);
  const keepScroll = viewport.scrollLeft;
  const cssW = Math.max(viewport.clientWidth, PAD.l + rec.nScans * c.px + PAD.r);
  viewport.style.height = height + "px";
  canvas.style.width = cssW + "px";
  canvas.style.height = height + "px";
  canvas.width = Math.floor(cssW * dpr);
  canvas.height = Math.floor(height * dpr);
  viewport.scrollLeft = keepScroll;
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  const sx = (scan) => PAD.l + scan * c.px;
  const [scan0, scan1] = visibleScans(viewport, rec, c.px);
  const dataMax = c.autoY ? visibleMax(rec, scan0, scan1) : item.fileMax;
  const ymax = Math.max(1, dataMax / c.yZoom);
  const plotH = TRACE_H;
  const sy = (rfu) => PAD.t + (1 - Math.min(rfu, ymax) / ymax) * plotH;

  ctx.fillStyle = "#fff";
  ctx.fillRect(0, 0, cssW, height);

  const spacing = rec.ploc.length > 1
    ? (rec.ploc[rec.ploc.length - 1] - rec.ploc[0]) / (rec.ploc.length - 1)
    : 12;
  const showLetters = c.px * spacing >= 7;
  const showBars = c.px * spacing >= 3;

  if (c.overlay && item.pred) {
    for (let i = 0; i < rec.ploc.length; i++) {
      const sc = rec.ploc[i];
      if (sc < scan0 - spacing || sc > scan1 + spacing) continue;
      const reason = item.pred[i].reason;
      const fill = REASON_FILL[reason];
      if (!fill) continue;
      const left = i > 0 ? (rec.ploc[i - 1] + sc) / 2 : sc - spacing / 2;
      const right = i < rec.ploc.length - 1 ? (sc + rec.ploc[i + 1]) / 2 : sc + spacing / 2;
      ctx.fillStyle = fill;
      ctx.fillRect(sx(left), PAD.t, Math.max(1, sx(right) - sx(left)), plotH);
    }
  }

  for (const lab of item.labels) {
    const a = rec.ploc[lab.start - 1];
    const b = rec.ploc[lab.end - 1];
    if (b < scan0 || a > scan1) continue;
    ctx.fillStyle = CLASS_FILL[lab.class] || "rgba(0,0,0,0.2)";
    ctx.fillRect(sx(a) - 2, PAD.t + plotH - 8, Math.max(4, sx(b) - sx(a) + 4), 8);
  }

  if (item._drag) {
    const a = Math.min(item._drag.a, item._drag.b);
    const b = Math.max(item._drag.a, item._drag.b);
    ctx.fillStyle = "rgba(30,90,180,0.18)";
    ctx.fillRect(sx(rec.ploc[a]), PAD.t, sx(rec.ploc[b]) - sx(rec.ploc[a]), plotH);
  }

  if (showBars) {
    ctx.strokeStyle = "rgba(80,80,80,0.28)";
    ctx.lineWidth = 0.6;
    for (let i = 0; i < rec.ploc.length; i++) {
      const sc = rec.ploc[i];
      if (sc < scan0 || sc > scan1) continue;
      const x = sx(sc);
      ctx.beginPath();
      ctx.moveTo(x, PAD.t);
      ctx.lineTo(x, PAD.t + plotH);
      ctx.stroke();
    }
  }

  ctx.setLineDash([5, 4]);
  ctx.lineWidth = 1;
  const lines = [
    [item.hq, "#555555", "HQ"],
    [c.alpha * item.hq, "#2a9d8f", "α"],
    [c.beta * item.hq, "#cc4444", "β"],
  ];
  const x0 = sx(scan0);
  const x1 = sx(scan1);
  for (const [y, color] of lines) {
    if (y <= 0) continue;
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(x0, sy(y));
    ctx.lineTo(x1, sy(y));
    ctx.stroke();
  }
  ctx.setLineDash([]);

  for (const b of "GATC") {
    const ch = rec.channels[b];
    if (!ch) continue;
    ctx.beginPath();
    ctx.strokeStyle = DYE[b];
    ctx.lineWidth = 1.05;
    let started = false;
    const step = c.px >= 0.8 ? 1 : Math.max(1, Math.round(1 / c.px));
    for (let s = scan0; s <= scan1; s += step) {
      const x = sx(s);
      const y = sy(ch[s] || 0);
      if (!started) { ctx.moveTo(x, y); started = true; }
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  ctx.font = "10px ui-monospace, Consolas, monospace";
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.fillStyle = "#666";
  for (let i = 0; i < rec.ploc.length; i++) {
    const pos = i + 1;
    if (pos % 10 !== 0 && pos !== 1) continue;
    const sc = rec.ploc[i];
    if (sc < scan0 || sc > scan1) continue;
    ctx.fillText(String(pos), sx(sc), PAD.t - 2);
  }

  const rows = [
    { label: "ABI", letters: rec.seq, kind: "abi" },
    { label: "v0.2", letters: item.v02, kind: "v02" },
  ];
  if (item.chromas) rows.push({ label: "Chr", letters: item.chromas.map, kind: "chr" });

  ctx.font = "11px ui-monospace, Consolas, monospace";
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  const row0 = PAD.t + plotH + 4;
  if (showLetters) {
    for (let r = 0; r < rows.length; r++) {
      const row = rows[r];
      const y = row0 + r * ROW_H;
      for (let i = 0; i < rec.ploc.length; i++) {
        const sc = rec.ploc[i];
        if (sc < scan0 || sc > scan1) continue;
        const ch = row.letters[i] || "";
        ctx.fillStyle = letterColor(row.kind, ch, rec.seq[i]);
        ctx.fillText(ch, sx(sc), y);
      }
    }
  }

  const sl = viewport.scrollLeft;
  const vw = viewport.clientWidth;
  ctx.fillStyle = "rgba(255,255,255,0.94)";
  ctx.fillRect(sl, 0, PAD.l - 2, height);
  ctx.strokeStyle = "rgba(0,0,0,0.12)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(sl + PAD.l - 2, PAD.t);
  ctx.lineTo(sl + PAD.l - 2, PAD.t + plotH);
  ctx.stroke();

  ctx.fillStyle = "#333";
  ctx.font = "10px ui-monospace, Consolas, monospace";
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const v of niceTicks(ymax)) {
    const y = sy(v);
    ctx.fillText(String(Math.round(v)), sl + PAD.l - 6, y);
    ctx.strokeStyle = "rgba(0,0,0,0.18)";
    ctx.beginPath();
    ctx.moveTo(sl + PAD.l - 6, y);
    ctx.lineTo(sl + PAD.l - 2, y);
    ctx.stroke();
  }
  ctx.fillStyle = "#888";
  ctx.textAlign = "center";
  ctx.save();
  ctx.translate(sl + 10, PAD.t + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("RFU", 0, 0);
  ctx.restore();

  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  const right = sl + vw - 8;
  for (const [y, color, lab] of lines) {
    if (y <= 0 || y > ymax * 1.05) continue;
    ctx.fillStyle = color;
    ctx.fillText(`${lab} ${y.toFixed(0)}`, right - 52, sy(y));
  }

  ctx.font = "10px ui-monospace, Consolas, monospace";
  ctx.textAlign = "right";
  ctx.textBaseline = "top";
  for (let r = 0; r < rows.length; r++) {
    ctx.fillStyle = "#555";
    ctx.fillText(rows[r].label, sl + PAD.l - 6, row0 + r * ROW_H);
  }
}

function letterColor(kind, ch, abi) {
  if (!ch || ch === "-") return "#bbb";
  if (kind === "abi") {
    return /[ACGT]/i.test(ch) ? "#666" : "#aa0000";
  }
  if (kind === "v02") {
    if (ch === "N") return "#993399";
    if (ch === ch.toLowerCase()) return "#1a7a72";
    return "#111";
  }
  if (ch.toUpperCase() === "N") return "#2a9d8f";
  if (abi && ch.toUpperCase() !== abi.toUpperCase() && /[ACGT]/i.test(ch) && /[ACGT]/i.test(abi)) {
    return "#aa0000";
  }
  return "#1a5f8a";
}

function drawAll() {
  syncOutputs();
  for (const item of files) drawOne(item);
}

function reclassifyAll() {
  for (const item of files) runClassify(item);
  drawAll();
}

function eventScan(item, ev) {
  const rect = item.viewport.getBoundingClientRect();
  const x = ev.clientX - rect.left + item.viewport.scrollLeft;
  const px = controls().px;
  const scan = (x - PAD.l) / px;
  return Math.max(0, Math.min(item.rec.nScans - 1, scan));
}

function syncScrollFrom(src) {
  if (!controls().lock || lockBusy) return;
  lockBusy = true;
  const maxSrc = src.viewport.scrollWidth - src.viewport.clientWidth;
  const frac = maxSrc > 0 ? src.viewport.scrollLeft / maxSrc : 0;
  for (const item of files) {
    if (item === src) continue;
    const max = item.viewport.scrollWidth - item.viewport.clientWidth;
    item.viewport.scrollLeft = frac * Math.max(0, max);
  }
  lockBusy = false;
}

function jumpHQ(item) {
  const i = Math.max(0, (item.island[0] || 1) - 1);
  const px = controls().px;
  item.viewport.scrollLeft = Math.max(0, PAD.l + item.rec.ploc[i] * px - item.viewport.clientWidth * 0.15);
}

function addPanel(item) {
  const panel = document.createElement("div");
  panel.className = "panel";
  panel.innerHTML = `
    <div class="panel-bar">
      <span class="name"></span>
      <span class="meta"></span>
      <span class="hover"></span>
      <button type="button" class="jump">HQ body</button>
      <button type="button" class="close">close</button>
    </div>
    <div class="viewport"><canvas></canvas></div>
    <ul class="labellist"></ul>`;
  panel.querySelector(".name").textContent = item.name;
  item.hoverEl = panel.querySelector(".hover");
  item.listEl = panel.querySelector(".labellist");
  panel.querySelector(".close").addEventListener("click", () => {
    const i = files.indexOf(item);
    if (i >= 0) files.splice(i, 1);
    panel.remove();
    els.drop.classList.toggle("has-files", files.length > 0);
  });
  panel.querySelector(".jump").addEventListener("click", () => jumpHQ(item));
  const viewport = panel.querySelector(".viewport");
  const canvas = panel.querySelector("canvas");
  item.viewport = viewport;
  item.canvas = canvas;
  item.panel = panel;
  item._lastPx = controls().px;
  refreshLabelList(item);
  refreshStats(item);

  let scrollTick = 0;
  viewport.addEventListener("scroll", () => {
    if (scrollTick) return;
    scrollTick = requestAnimationFrame(() => {
      scrollTick = 0;
      drawOne(item);
      syncScrollFrom(item);
    });
  });
  viewport.addEventListener("wheel", (ev) => {
    if (ev.ctrlKey) {
      ev.preventDefault();
      const next = Math.min(8, Math.max(0.25, Number(els.yZoom.value) - ev.deltaY * 0.002));
      els.yZoom.value = String(next);
      drawAll();
    } else if (ev.shiftKey) {
      ev.preventDefault();
      const scan = viewport.scrollLeft / (Number(els.px.value) || 1);
      const next = Math.min(4, Math.max(0.3, Number(els.px.value) - ev.deltaY * 0.002));
      els.px.value = String(next);
      drawAll();
      viewport.scrollLeft = scan * next;
    } else if (Math.abs(ev.deltaY) > Math.abs(ev.deltaX)) {
      ev.preventDefault();
      viewport.scrollLeft += ev.deltaY;
    }
  }, { passive: false });

  canvas.addEventListener("mousemove", (ev) => {
    const scan = eventScan(item, ev);
    const bi = scanToBase(item.rec, scan);
    const d = item.pred[bi];
    const ch = item.chromas ? item.chromas.map[bi] : "–";
    item.hoverEl.textContent =
      `pos ${d.pos}  ABI ${item.rec.seq[bi]}  v0.2 ${item.v02[bi]}  Chr ${ch}  ` +
      `${d.reason}  h=${d.top_h.toFixed(0)} RFU  per=${d.periodicity.toFixed(2)}  ` +
      `valley=${d.valley_ratio.toFixed(2)}  z=${d.spike_z.toFixed(2)}`;
    if (item._drag) {
      item._drag.b = bi;
      drawOne(item);
    }
  });
  canvas.addEventListener("mousedown", (ev) => {
    if (ev.button !== 0) return;
    const bi = scanToBase(item.rec, eventScan(item, ev));
    item._drag = { a: bi, b: bi };
  });

  els.panels.appendChild(panel);
  requestAnimationFrame(() => drawOne(item));
}

function refreshLabelList(item) {
  if (!item.listEl) return;
  item.listEl.textContent = "";
  item.labels.forEach((lab, i) => {
    const li = document.createElement("li");
    li.textContent = `${lab.start}–${lab.end}  ${lab.class}  ${lab.grade}`;
    const rm = document.createElement("button");
    rm.type = "button";
    rm.textContent = "remove";
    rm.addEventListener("click", () => {
      item.labels.splice(i, 1);
      refreshLabelList(item);
      drawOne(item);
    });
    li.appendChild(rm);
    item.listEl.appendChild(li);
  });
}

function attachFastaTo(item, faRec) {
  const al = alignFastaToPbas(item.rec.seq, faRec.seq);
  item.chromas = {
    name: faRec.name,
    seq: faRec.seq,
    map: al.map,
    score: al.score,
    rc: al.rc,
  };
  refreshStats(item);
  drawOne(item);
}

function matchFastaToOpen(faRec) {
  const hint = nameHint(faRec.name);
  if (hint) {
    const hit = files.find((f) => nameHint(f.name) === hint);
    if (hit) return hit;
  }
  return files.length === 1 ? files[0] : files[files.length - 1] || null;
}

function isAb1Name(name) {
  return /\.(ab1|abif)$/i.test(name);
}

function isFastaName(name) {
  return /\.(fa|fasta|fas|seq)$/i.test(name);
}

els.dlg.addEventListener("close", () => {
  pendingLabel = null;
});

els.dlgOk.addEventListener("click", (ev) => {
  ev.preventDefault();
  if (!pendingLabel) return;
  const { item, start, end } = pendingLabel;
  item.labels.push({
    file: item.name,
    start,
    end,
    called: item.rec.seq.slice(start - 1, end),
    class: els.dlgClass.value,
    grade: "M1",
    evidence: ["MANUAL"],
    job: "call_reliability",
    date: new Date().toISOString().slice(0, 10),
  });
  pendingLabel = null;
  els.dlg.close();
  refreshLabelList(item);
  drawOne(item);
});

els.exportBtn.addEventListener("click", () => {
  const rows = files.flatMap((f) => f.labels);
  const blob = new Blob([JSON.stringify(rows, null, 2)], { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "sangerqc_labels.json";
  a.click();
  URL.revokeObjectURL(a.href);
});

async function loadAb1(f) {
  const rec = parseABIF(await f.arrayBuffer());
  if (!rec.seq.length || !rec.ploc.length) {
    throw new Error(`${f.name}: no PBAS/PLOC`);
  }
  const item = {
    id: f.name + ":" + files.length + ":" + Date.now(),
    name: f.name,
    rec,
    fileMax: fileMax(rec),
    labels: [],
  };
  runClassify(item);
  files.push(item);
  els.drop.classList.add("has-files");
  addPanel(item);
  const hint = nameHint(f.name);
  const pi = pendingFastas.findIndex((p) => !hint || nameHint(p.name) === hint);
  if (pi >= 0) {
    attachFastaTo(item, pendingFastas.splice(pi, 1)[0]);
  }
  showError("");
  return item;
}

async function loadFasta(f) {
  const recs = parseFasta(await f.text());
  if (!recs.length) throw new Error(`${f.name}: no FASTA sequence`);
  for (const faRec of recs) {
    const item = matchFastaToOpen(faRec);
    if (!item) {
      pendingFastas.push(faRec);
      showError("FASTA held until an .ab1 is open (match by well name if possible).");
      continue;
    }
    attachFastaTo(item, faRec);
    showError("");
  }
}

async function loadFile(f) {
  if (!f) return;
  try {
    if (isAb1Name(f.name)) await loadAb1(f);
    else if (isFastaName(f.name)) await loadFasta(f);
    else await loadAb1(f);
  } catch (e) {
    showError(String(e.message || e));
  }
}

async function loadList(list) {
  const arr = Array.from(list || []);
  const ab1s = arr.filter((f) => isAb1Name(f.name));
  const fas = arr.filter((f) => isFastaName(f.name));
  const rest = arr.filter((f) => !isAb1Name(f.name) && !isFastaName(f.name));
  for (const f of ab1s.concat(rest)) await loadFile(f);
  for (const f of fas) await loadFile(f);
  els.file.value = "";
}

els.file.addEventListener("change", (ev) => loadList(ev.target.files || []));
["dragenter", "dragover"].forEach((name) => {
  els.drop.addEventListener(name, (ev) => {
    ev.preventDefault();
    els.drop.classList.add("drag");
  });
});
["dragleave", "drop"].forEach((name) => {
  els.drop.addEventListener(name, (ev) => {
    ev.preventDefault();
    els.drop.classList.remove("drag");
  });
});
els.drop.addEventListener("drop", (ev) => loadList(ev.dataTransfer.files || []));

["px", "yZoom", "alpha", "beta"].forEach((id) => {
  els[id].addEventListener("input", () => {
    if (id === "alpha" || id === "beta") {
      reclassifyAll();
      return;
    }
    if (id === "px") {
      const px = Number(els.px.value) || 1;
      for (const item of files) {
        const oldPx = item._lastPx || px;
        const scan = item.viewport.scrollLeft / oldPx;
        item._lastPx = px;
        drawOne(item);
        item.viewport.scrollLeft = scan * px;
      }
      syncOutputs();
      return;
    }
    drawAll();
  });
});
els.autoY.addEventListener("change", drawAll);
els.overlay.addEventListener("change", drawAll);
window.addEventListener("resize", drawAll);
window.addEventListener("mouseup", () => {
  const item = files.find((f) => f._drag);
  if (!item) return;
  const a0 = Math.min(item._drag.a, item._drag.b);
  const b0 = Math.max(item._drag.a, item._drag.b);
  item._drag = null;
  drawOne(item);
  if (b0 === a0) return;
  const a = a0 + 1, b = b0 + 1;
  pendingLabel = { item, start: a, end: b };
  els.dlgRange.textContent = `${item.name}  ${a}–${b}`;
  els.dlg.showModal();
});

CLASSES.forEach((c) => {
  const opt = document.createElement("option");
  opt.value = c;
  opt.textContent = c;
  els.dlgClass.appendChild(opt);
});

syncOutputs();
window.sangerQC = { loadFile, loadAb1, attachFastaTo, drawAll, files, reclassifyAll, jumpHQ };
