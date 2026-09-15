const DYE = { A: "#2ca02c", C: "#1f77b4", G: "#111111", T: "#d62728" };
const CSS_H = 280;
const PAD = { l: 48, r: 8, t: 18, b: 28 };

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
  alpha: document.getElementById("alpha"),
  alphaOut: document.getElementById("alphaOut"),
  beta: document.getElementById("beta"),
  betaOut: document.getElementById("betaOut"),
};

const files = []; // { id, name, rec, hq, heights, viewport, canvas }

function showError(msg) {
  els.err.hidden = !msg;
  els.err.textContent = msg || "";
}

function peakHeights(rec) {
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

function hqBody(heights) {
  if (!heights.length) return 1;
  const mid = heights.slice(
    Math.floor(heights.length * 0.25),
    Math.ceil(heights.length * 0.75)
  );
  const s = (mid.length ? mid : heights).slice().sort((a, b) => a - b);
  return s[Math.floor(s.length / 2)] || 1;
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

function controls() {
  return {
    px: Number(els.px.value) || 1,
    yZoom: Number(els.yZoom.value) || 1,
    autoY: els.autoY.checked,
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

function drawOne(item) {
  const rec = item.rec;
  const c = controls();
  const canvas = item.canvas;
  const viewport = item.viewport;
  const dpr = window.devicePixelRatio || 1;
  const cssW = Math.max(viewport.clientWidth, PAD.l + rec.nScans * c.px + PAD.r);
  const cssH = CSS_H;
  canvas.style.width = cssW + "px";
  canvas.style.height = cssH + "px";
  canvas.width = Math.floor(cssW * dpr);
  canvas.height = Math.floor(cssH * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  const sx = (scan) => PAD.l + scan * c.px;
  const [scan0, scan1] = visibleScans(viewport, rec, c.px);
  const dataMax = c.autoY ? visibleMax(rec, scan0, scan1) : item.fileMax;
  const ymax = Math.max(1, dataMax / c.yZoom);
  const plotH = cssH - PAD.t - PAD.b;
  const sy = (rfu) => PAD.t + (1 - Math.min(rfu, ymax) / ymax) * plotH;

  ctx.fillStyle = "#fff";
  ctx.fillRect(0, 0, cssW, cssH);

  const spacing = rec.ploc.length > 1
    ? (rec.ploc[rec.ploc.length - 1] - rec.ploc[0]) / (rec.ploc.length - 1)
    : 12;
  const showLetters = c.px * spacing >= 7;
  const showBars = c.px * spacing >= 3;

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
    [item.hq, "#555555"],
    [c.alpha * item.hq, "#2a9d8f"],
    [c.beta * item.hq, "#cc4444"],
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

  if (showLetters) {
    ctx.font = "11px ui-monospace, Consolas, monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillStyle = "#222";
    for (let i = 0; i < rec.ploc.length; i++) {
      const sc = rec.ploc[i];
      if (sc < scan0 || sc > scan1) continue;
      ctx.fillText(rec.seq[i] || "", sx(sc), PAD.t + plotH + 4);
    }
  }

  ctx.fillStyle = "#666";
  ctx.font = "11px system-ui, sans-serif";
  ctx.textAlign = "left";
  ctx.textBaseline = "top";
  const yMode = c.autoY ? "auto-Y" : "full-Y";
  ctx.fillText(
    `scans ${scan0}–${scan1}  ymax ${ymax.toFixed(0)} RFU (${yMode})  HQ ${item.hq.toFixed(0)}`,
    viewport.scrollLeft + 8,
    2
  );
}

function drawAll() {
  syncOutputs();
  for (const item of files) drawOne(item);
}

function addPanel(item) {
  const panel = document.createElement("div");
  panel.className = "panel";
  panel.innerHTML = `
    <div class="panel-bar">
      <span class="name"></span>
      <span class="meta"></span>
      <button type="button" class="close">close</button>
    </div>
    <div class="viewport"><canvas></canvas></div>`;
  panel.querySelector(".name").textContent = item.name;
  panel.querySelector(".meta").textContent =
    `${item.rec.seq.length} bases · ${item.rec.nScans} scans · ${item.rec.fwo}`;
  panel.querySelector(".close").addEventListener("click", () => {
    const i = files.indexOf(item);
    if (i >= 0) files.splice(i, 1);
    panel.remove();
    els.drop.classList.toggle("has-files", files.length > 0);
  });
  const viewport = panel.querySelector(".viewport");
  const canvas = panel.querySelector("canvas");
  item.viewport = viewport;
  item.canvas = canvas;
  item._lastPx = controls().px;
  let scrollTick = 0;
  viewport.addEventListener("scroll", () => {
    if (scrollTick) return;
    scrollTick = requestAnimationFrame(() => {
      scrollTick = 0;
      drawOne(item);
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
  els.panels.appendChild(panel);
  drawOne(item);
}

async function loadFile(f) {
  if (!f) return;
  try {
    const rec = parseABIF(await f.arrayBuffer());
    if (!rec.seq.length || !rec.ploc.length) {
      throw new Error(`${f.name}: no PBAS/PLOC`);
    }
    const heights = peakHeights(rec);
    const item = {
      id: f.name + ":" + files.length + ":" + Date.now(),
      name: f.name,
      rec,
      heights,
      hq: hqBody(heights),
      fileMax: fileMax(rec),
    };
    files.push(item);
    els.drop.classList.add("has-files");
    addPanel(item);
    showError("");
  } catch (e) {
    showError(String(e.message || e));
  }
}

async function loadList(list) {
  for (const f of list) await loadFile(f);
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
window.addEventListener("resize", drawAll);
syncOutputs();
window.sangerQC = { loadFile, drawAll, files };
