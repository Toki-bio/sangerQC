"""Forward + reverse reads of the same template, real heterozygotes with manual genotypes.

    python validation/eval_multiread.py PLATE1_DIR PLATE2_DIR REPORT.docx [--cache FILE]

Data: two plates of nested-PCR Sanger reads (SeqStudio), 23 samples x 2 SNPs x forward/reverse.
Plate 1 = rs10156191 (C/T), plate 2 = rs1049793 (G/C) (identified from the data: the double-peak
count matches the report's heterozygote count). Layout inferred from the data, NOT given: forward block column-major, 8 rows,
sample s at row (s-1) % 8, column (s-1) // 8 within the block; the reverse block holds the same samples,
and the order of its three columns is chosen from forward/reverse agreement of the double-peak calls alone
(no report used); on plate 1 this finds columns 5 and 6 swapped.
Truth: genotypes in the report's table (called by eye in Lasergene/SnapGene = manual grade M1).
Only sample numbers and genotypes are read from the report.

Per site, methods compared for zygosity (heterozygous vs homozygous):
  caller   the basecaller's own letter is an IUPAC code
  rule     co-located second peak >= 0.33 of the primary (sangerQC v0.3 rule)
  bayes1   single read, p_mixed from the per-read model > 0.5
  joint    forward + reverse combined (multiread.py), flat prior
  joint+SNP  same with a biallelic prior (alleles = the two commonest primary bases at the site)
"""
import glob
import itertools
import math
import os
import pickle
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sangerqc import bayes  # noqa: E402
from sangerqc.concordance import map_to_reference  # noqa: E402
from sangerqc.multiread import COMP, analyze_read, joint_genotype, read_likelihood  # noqa: E402

ROWS = "ABCDEFGH"
# plate -> (reference well, reference column of the SNP, report column, name)
SITES = {
    "p1": [("B", 1, 108, "rs10156191")],
    "p2": [("F", 7, 191, "rs1049793"), ("F", 7, 148, "site148 (not in report)")],
}
FWD_FIRST = {"p1": 1, "p2": 10}     # first column of the forward block; reverse block = +/-3
REV_FIRST = {"p1": 4, "p2": 7}


def wells(plate_dir):
    out = {}
    for f in glob.glob(plate_dir + "/*.ab1"):
        m = re.match(r"([A-H])(\d+)_", os.path.basename(f))
        out[(m.group(1), int(m.group(2)))] = f
    return out


def load_truth(docx_path):
    import docx
    d = docx.Document(docx_path)
    for t in d.tables:
        if t.rows[0].cells[0].text.strip() == "№":
            tr = {}
            for r in t.rows[1:]:
                c = [x.text.strip() for x in r.cells]
                tr[int(c[0])] = {"rs1049793": c[4].replace(" ", ""), "rs10156191": c[5].replace(" ", "")}
            return tr
    raise SystemExit("genotype table not found")


def infer_reverse_order(plate, analyses, maps, j):
    """Order of the three reverse-block columns that makes forward and reverse zygosity calls agree most.
    Uses single-read calls at the plate's SNP column only; returns (best order, table of scores)."""
    def call(k):
        if k not in analyses or j not in maps[k]:
            return None
        i, _ = maps[k][j]
        a = analyses[k]
        if a.pred[i + 1]["bad"]:
            return None
        o = bayes.observe(a.raw["secondary_records"][i]["probe"])
        return bayes.posterior(bayes.log10_bf(o, a.model), a.pi) > 0.5

    cols = list(range(REV_FIRST[plate], REV_FIRST[plate] + 3))
    table = []
    for perm in itertools.permutations(cols):
        agree = n = 0
        for s in range(1, 24):
            row, off = ROWS[(s - 1) % 8], (s - 1) // 8
            f = call((plate, row, FWD_FIRST[plate] + off))
            r = call((plate, row, perm[off]))
            if f is not None and r is not None:
                n += 1
                agree += f == r
        table.append((agree, n, perm))
    table.sort(reverse=True)
    return table[0][2], table


def is_het(g):
    a, b = g.split("/")
    return a != b


def main(p1, p2, report, cache=None):
    truth = load_truth(report)
    dirs = {"p1": p1, "p2": p2}
    if cache and os.path.exists(cache):
        analyses = pickle.load(open(cache, "rb"))
    else:
        analyses = {}
        for plate in dirs:
            for (r, c), f in sorted(wells(dirs[plate]).items()):
                if len(open(f, "rb").read()) and True:
                    try:
                        ra = analyze_read(f, name=f"{plate}:{r}{c}")
                    except Exception as e:      # failed wells (5 bases) have no usable trace
                        continue
                    if len(ra.raw["seq"]) >= 150:
                        analyses[(plate, r, c)] = ra
        if cache:
            pickle.dump(analyses, open(cache, "wb"))
    print(f"{len(analyses)} usable reads")

    results = []
    for plate, sites in SITES.items():
        ref_r, ref_c, _, _ = sites[0]
        ref = analyses[(plate, ref_r, ref_c)].raw["seq"]
        maps = {k: map_to_reference(a.raw["seq"], ref) for k, a in analyses.items() if k[0] == plate}
        order, table = infer_reverse_order(plate, analyses, maps, sites[0][2] - 1)
        print(f"\n{plate}: reverse-column order by forward/reverse agreement (no report used):",
              ", ".join(f"{perm}: {a}/{n}" for a, n, perm in table[:3]), "-> chosen", order)
        for (_, _, col, sname) in sites:
            j = col - 1
            prim = Counter()
            for k, mp in maps.items():
                if j in mp:
                    i, st = mp[j]
                    p = analyses[k].raw["secondary_records"][i]["probe"]["primary"]
                    prim[p if st == "+" else COMP[p]] += 1
            alleles = [b for b, _ in prim.most_common(2)]
            print(f"\n=== {sname} ({plate}, reference column {col}); primaries at the site {dict(prim)}; alleles {alleles}")
            pi_site = float(np.mean([analyses[k].pi for k in maps]))
            for s in range(1, 24):
                row, off = ROWS[(s - 1) % 8], (s - 1) // 8
                fk = (plate, row, FWD_FIRST[plate] + off)
                rk = (plate, row, order[off])
                reads = []
                for k in (fk, rk):
                    if k in analyses and j in maps[k]:
                        i, st = maps[k][j]
                        reads.append((analyses[k], i, st, "FR"))
                if not reads:
                    continue
                t = truth.get(s, {}).get(sname)
                rec = dict(plate=plate, site=sname, sample=s, truth_het=None if t is None else is_het(t), n_reads=len(reads))
                # per-read quantities
                per = []
                for ra, i, st, _ in reads:
                    pr = ra.raw["secondary_records"][i]["probe"]
                    o = bayes.observe(pr)
                    caller = ra.raw["seq"][i]
                    lbf = bayes.log10_bf(o, ra.model)
                    per.append(dict(caller_het=caller in "RYSWKM", rule=bool(o and o[0] >= 0.33 and o[0] >= 0.20),
                                    p1=bayes.posterior(lbf, ra.pi), bad=ra.pred[i + 1]["bad"]))
                rec["per_read"] = per
                for name, kw in (("joint", {}), ("joint+SNP", dict(alleles=alleles))):
                    rec[name] = joint_genotype(reads, pi=pi_site, **kw)["p_het"]
                rec["joint_w05"] = joint_genotype(reads, pi=pi_site, weights={"FR": 0.5})["p_het"]
                for tag, rd in (("F", reads[:1]), ("R", reads[1:])):
                    if rd:
                        rec[f"only{tag}"] = joint_genotype(rd, pi=pi_site)["p_het"]
                results.append(rec)

    # ---- layout check without the truth table: forward and reverse agree on zygosity?
    print("\nLayout check (no truth used): forward vs reverse read of the same sample, single-read p_mixed > 0.5")
    for plate, sites in SITES.items():
        for (_, _, col, sname) in sites:
            pairs = [r for r in results if r["site"] == sname and len(r["per_read"]) == 2
                     and not any(p["bad"] for p in r["per_read"])]
            agree = sum((r["per_read"][0]["p1"] > .5) == (r["per_read"][1]["p1"] > .5) for r in pairs)
            print(f"  {sname}: {agree}/{len(pairs)} pairs agree")

    # ---- accuracy against the report
    print("\nZygosity vs the report (heterozygous = positive); samples with the SNP covered")
    for sname in ("rs10156191", "rs1049793"):
        rs = [r for r in results if r["site"] == sname and r["truth_het"] is not None]
        if not rs:
            continue
        pos = sum(r["truth_het"] for r in rs)
        print(f"\n{sname}: {len(rs)} samples, {pos} heterozygous in the report")

        def tab(label, pred):
            tp = sum(pred(r) and r["truth_het"] for r in rs)
            fp = sum(pred(r) and not r["truth_het"] for r in rs)
            fn = sum((not pred(r)) and r["truth_het"] for r in rs)
            tn = sum((not pred(r)) and not r["truth_het"] for r in rs)
            print(f"  {label:34s} TP {tp:2d} FP {fp:2d} FN {fn:2d} TN {tn:2d}  accuracy {(tp+tn)/len(rs):.3f}")

        # per-read level (every usable read counts)
        def rd_tab(label, key):
            tp = fp = fn = tn = 0
            for r in rs:
                for p in r["per_read"]:
                    if p["bad"]:
                        continue
                    call = p[key] if key != "p1" else p["p1"] > 0.5
                    if call and r["truth_het"]: tp += 1
                    elif call: fp += 1
                    elif r["truth_het"]: fn += 1
                    else: tn += 1
            print(f"  {label:34s} TP {tp:2d} FP {fp:2d} FN {fn:2d} TN {tn:2d}  accuracy {(tp+tn)/(tp+fp+fn+tn):.3f}  (reads)")
        rd_tab("basecaller IUPAC, per read", "caller_het")
        rd_tab("ratio rule >=0.33, per read", "rule")
        rd_tab("Bayes single read p>0.5", "p1")
        tab("joint F+R flat prior p>0.5", lambda r: r["joint"] > .5)
        tab("joint F+R biallelic prior p>0.5", lambda r: r["joint+SNP"] > .5)
        tab("  weight 0.5, flat prior", lambda r: r["joint_w05"] > .5)
        for key, lab in (("joint", "flat"), ("joint+SNP", "biallelic"), ("joint_w05", "flat, w=0.5"),
                         ("onlyF", "forward only"), ("onlyR", "reverse only")):
            ll = []
            for r in rs:
                if key in r:
                    p = min(max(r[key], 1e-4), 1 - 1e-4)
                    ll.append(-math.log(p) if r["truth_het"] else -math.log(1 - p))
            if ll:
                print(f"  log-loss {lab:14s} {np.mean(ll):.3f}  (n={len(ll)})")
        wrong = [(r["sample"], round(r["joint+SNP"], 3), [round(p["p1"], 2) for p in r["per_read"]], r["truth_het"])
                 for r in rs if (r["joint+SNP"] > .5) != r["truth_het"]]
        print("  joint+SNP errors (sample, p_het, single-read p, truth het):", wrong)
    # site 148, no truth
    rs = [r for r in results if r["site"].startswith("site148")]
    print(f"\nsite 148 (not in the report): joint p_het > 0.5 in {sum(r['joint'] > .5 for r in rs)}/{len(rs)} samples; "
          f"forward/reverse agree on {sum((r['per_read'][0]['p1'] > .5) == (r['per_read'][1]['p1'] > .5) for r in rs if len(r['per_read']) == 2)}"
          f"/{sum(len(r['per_read']) == 2 for r in rs)} pairs")
    pickle.dump(results, open((cache or "multiread_results") + ".results", "wb"))


if __name__ == "__main__":
    a = sys.argv[1:]
    cache = a[a.index("--cache") + 1] if "--cache" in a else None
    main(a[0], a[1], a[2], cache)
