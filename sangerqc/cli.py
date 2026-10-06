"""Command line: ab1 -> FASTA + per-position decision table (v0.3).

    python -m sangerqc.cli READ1.ab1 [READ2.ab1 ...] --out DIR
        [--forward SEQ] [--reverse SEQ] [--reference REF.fasta]

With --reference and two or more reads, secondary peaks are checked for
reproduction across reads (concordance.py) and reproduced sites are coded IUPAC
even below the single-read 0.33 threshold.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO

from .concordance import reference_status, reproduced_sites
from .export import FIELDS, to_sequence
from . import __version__
from .secondary import COLOCATE_MAX, RATIO_CALL, RATIO_MIN, background_rate
from .v0_3 import classify


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ab1", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--forward")
    ap.add_argument("--reverse")
    ap.add_argument("--reference", help="FASTA; first record is the shared coordinate for concordance")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    items, paths = [], {}
    for f in a.ab1:
        pred, raw = classify(f, forward=a.forward, reverse=a.reverse)
        items.append((Path(f).stem, pred, raw))
        paths[Path(f).stem] = Path(f)

    reproduced = {name: {} for name, _, _ in items}
    sites, ref = {}, None
    if a.reference:
        ref = str(next(SeqIO.parse(a.reference, "fasta")).seq).upper()
    if ref and len(items) > 1:
        sites, reproduced = reproduced_sites(items, ref)

    summary = {}
    with open(out / "sangerqc_v0_3.fasta", "w", encoding="utf-8") as fa:
        for name, pred, raw in items:
            rs = reference_status(raw["seq"], ref) if ref else None
            seq, rows, span = to_sequence(pred, raw, reproduced=reproduced[name], ref_status=rs)
            d = dict(raw["diagnostics"])
            d["file"] = paths[name].name
            d["file_md5"] = hashlib.md5(paths[name].read_bytes()).hexdigest()
            d["kept_span"] = span
            d["length_out"] = len(seq)
            d["iupac_out"] = sum(c.upper() in "RYSWKM" for c in seq)
            d["lowercase_out"] = sum(c.islower() for c in seq)
            # control for the IUPAC calls: how often any non-called base passes the
            # co-location test at unflagged positions outside reproduced sites
            ordinary = [i for i in range(len(raw["seq"])) if not pred[i + 1]["bad"]
                        and (i + 1) not in reproduced[name]]
            d["colocation_pass_rate"] = background_rate(raw, raw["secondary_records"], ordinary)
            summary[name] = d
            fa.write(f">{name} sangerQC v0.3 kept {span}\n{seq}\n")
            with open(out / f"{name}_decisions.tsv", "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
                w.writeheader()
                w.writerows(rows)
    if sites:
        with open(out / "reproduced_sites.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["ref_pos", "pair", "read", "read_pos", "primary", "secondary", "ratio"])
            for j, s in sorted(sites.items()):
                for name, (i, p, q, r) in s["reads"].items():
                    w.writerow([j, "/".join(sorted(s["pair"])), name, i, p, q, r])
    summary["_run"] = {"sangerqc_version": __version__, "params": {
        "ratio_min": RATIO_MIN, "ratio_call": RATIO_CALL, "colocate_max": COLOCATE_MAX},
        "forward": a.forward, "reverse": a.reverse, "reference": a.reference}
    summary["_reproduced_sites"] = {j: {"pair": "/".join(sorted(s["pair"])),
                                         "also": ["/".join(sorted(p)) for p in s["also"]]}
                                     for j, s in sorted(sites.items())}
    with open(out / "summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, default=str)
    for name, d in summary.items():
        if name.startswith("_"):
            continue
        print(name, {k: d[k] for k in ("kept_span", "length_out", "iupac_out", "lowercase_out", "low_snr")})
    if sites:
        print("reproduced sites:", len(sites))


if __name__ == "__main__":
    main()
