"""Compare v0.2 and v0.3 bad flags against labelled zones (label 1 = bad, 0 = clean).

    python validation/compare_labels.py LABELS.json READ.ab1 [LABELS.json READ.ab1 ...]

Labels are the calibration-loop records (pos, label); unpublished data stays outside the repo.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sangerqc.v0_2 import classify as c2  # noqa: E402
from sangerqc.v0_3 import classify as c3  # noqa: E402


def score(pred, labels):
    tp = tn = fp = fn = 0
    for r in labels:
        bad = pred[r["pos"]]["bad"]
        if r["label"] == 1:
            tp += bad; fn += not bad
        else:
            tn += not bad; fp += bad
    return tp, fn, tn, fp


def main(args):
    for lab_path, ab1 in zip(args[0::2], args[1::2]):
        labels = json.load(open(lab_path))
        p2, _ = c2(ab1)
        p3, _ = c3(ab1)
        s2, s3 = score(p2, labels), score(p3, labels)
        name = labels[0].get("file", Path(ab1).stem)
        acc = lambda s: (s[0] + s[2]) / sum(s)
        print(f"{name}: v0.2 bad caught {s2[0]}/{s2[0]+s2[1]} clean kept {s2[2]}/{s2[2]+s2[3]} acc {acc(s2):.3f} | "
              f"v0.3 bad caught {s3[0]}/{s3[0]+s3[1]} clean kept {s3[2]}/{s3[2]+s3[3]} acc {acc(s3):.3f}")
        changed = [i for i in p2 if p2[i]["bad"] != p3[i]["bad"]]
        lab = {r["pos"]: r["label"] for r in labels}
        print("   flag changes:", [(i, p3[i]["reason"], lab.get(i, "-")) for i in changed])


if __name__ == "__main__":
    main(sys.argv[1:])
