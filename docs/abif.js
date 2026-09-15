/**
 * Minimal ABIF (.ab1) reader for the browser.
 * Extracts FWO_1, PBAS, PLOC, DATA9-12. Big-endian.
 */
(function (global) {
  const ENTRY = 28;

  function u16(v, o) { return v.getUint16(o, false); }
  function i16(v, o) { return v.getInt16(o, false); }
  function i32(v, o) { return v.getInt32(o, false); }

  function tagName(v, o) {
    return String.fromCharCode(v.getUint8(o), v.getUint8(o + 1), v.getUint8(o + 2), v.getUint8(o + 3));
  }

  function readEntry(v, o) {
    return {
      name: tagName(v, o),
      number: i32(v, o + 4),
      eltype: i16(v, o + 8),
      elsize: i16(v, o + 10),
      numels: i32(v, o + 12),
      datasize: i32(v, o + 16),
      dataoffset: i32(v, o + 20),
    };
  }

  function payloadOffset(e) {
    return e.datasize <= 4 ? null : e.dataoffset;
  }

  function readPayload(view, e) {
    const n = e.numels;
    const inline = e.datasize <= 4;
    const base = inline ? null : e.dataoffset;
    const t = e.eltype;
    const out = [];
    if (t === 2 || t === 18 || t === 19 || t === 1) {
      // char / p-string / c-string / byte
      if (inline) {
        const buf = new Uint8Array(4);
        new DataView(buf.buffer).setInt32(0, e.dataoffset, false);
        return String.fromCharCode(...buf).replace(/\0/g, "");
      }
      let s = "";
      const start = t === 18 ? base + 1 : base;
      const len = t === 18 ? view.getUint8(base) : n;
      for (let i = 0; i < len && start + i < view.byteLength; i++) {
        const c = view.getUint8(start + i);
        if (t === 19 && c === 0) break;
        s += String.fromCharCode(c);
      }
      return s;
    }
    for (let i = 0; i < n; i++) {
      if (inline) {
        // at most one 16/32-bit value in the offset field
        if (t === 4 || t === 3) out.push(t === 4 ? (e.dataoffset >> 16) : (e.dataoffset >>> 16));
        else if (t === 5) out.push(e.dataoffset);
        break;
      }
      const o = base + i * e.elsize;
      if (t === 4) out.push(i16(view, o));
      else if (t === 3) out.push(u16(view, o));
      else if (t === 5) out.push(i32(view, o));
      else if (t === 7) out.push(view.getFloat32(o, false));
      else out.push(i16(view, o));
    }
    return out;
  }

  function parseABIF(arrayBuffer) {
    const view = new DataView(arrayBuffer);
    const magic = String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3));
    if (magic !== "ABIF") throw new Error("Not an ABIF file (missing ABIF magic).");
    const dirMeta = readEntry(view, 6);
    const nDir = dirMeta.numels;
    const dirOff = dirMeta.dataoffset;
    const tags = {};
    for (let i = 0; i < nDir; i++) {
      const e = readEntry(view, dirOff + i * ENTRY);
      const key = e.name + e.number;
      tags[key] = e;
    }
    const get = (name) => {
      const e = tags[name];
      if (!e) return null;
      return readPayload(view, e);
    };
    const fwo = (get("FWO_1") || "GATC").toString().replace(/\0/g, "").slice(0, 4);
    const seq = (get("PBAS2") || get("PBAS1") || "").toString();
    const ploc = get("PLOC2") || get("PLOC1") || [];
    const data = {
      DATA9: get("DATA9") || [],
      DATA10: get("DATA10") || [],
      DATA11: get("DATA11") || [],
      DATA12: get("DATA12") || [],
    };
    const channels = {};
    ["DATA9", "DATA10", "DATA11", "DATA12"].forEach((k, i) => {
      channels[fwo[i] || "?"] = data[k];
    });
    return { fwo, seq, ploc, channels, nScans: (data.DATA9 || []).length };
  }

  global.parseABIF = parseABIF;
})(typeof window !== "undefined" ? window : globalThis);
