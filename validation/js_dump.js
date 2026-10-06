// node validation/js_dump.js FILE.ab1  -> JSON per position from the browser classifier
const fs = require("fs");
const path = require("path");
for (const f of ["abif.js", "classify.js", "compare.js"]) {
  eval(fs.readFileSync(path.join(__dirname, "..", "docs", f), "utf8"));
}
const buf = fs.readFileSync(process.argv[2]);
const rec = parseABIF(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
const res = classifyV03(rec, 0.35, 0.18);
const letters = v02Letters(rec.seq, res.pred);
console.log(JSON.stringify({
  hq: res.hq, lowSnr: res.lowSnr,
  pos: res.pred.map((p, i) => [p.bad, p.reason, p.sec.cls, p.sec.iupac, letters[i]]),
}));
