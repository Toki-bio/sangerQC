"""Python (sangerqc v0.3) vs browser JS (docs/classify.js) on the same .ab1 files.

    python validation/parity.py FILE.ab1 [...]
"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sangerqc.v0_3 import classify  # noqa: E402

for f in sys.argv[1:]:
    js = json.loads(subprocess.run(["node", str(Path(__file__).parent / "js_dump.js"), f],
                                   capture_output=True, text=True, check=True).stdout)
    pred, raw = classify(f)
    n = len(raw["seq"])
    diff = {"bad": [], "reason": [], "sec_cls": [], "iupac": []}
    for i in range(n):
        p, j = pred[i + 1], js["pos"][i]
        if p["bad"] != j[0]: diff["bad"].append(i + 1)
        if p["reason"] != j[1]: diff["reason"].append((i + 1, p["reason"], j[1]))
        if p["secondary"]["cls"] != j[2]: diff["sec_cls"].append((i + 1, p["secondary"]["cls"], j[2]))
        if (p["secondary"]["iupac"] or None) != j[3]: diff["iupac"].append(i + 1)
    print(Path(f).name[:10], "n", n, "HQ py %.0f js %.0f" % (pred[1]["hq_body"], js["hq"]),
          {k: (len(v), v[:4]) for k, v in diff.items()})
