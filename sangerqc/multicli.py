"""Several reads per template -> joint consensus with per-position evidence.

    python -m sangerqc.multicli --reference REF.fasta --out DIR \\
        --template T1=T1_fwd.ab1,T1_rev.ab1 --template T2=T2_a.ab1,T2_b.ab1,T2_c.ab1
    python -m sangerqc.multicli --reference REF.fasta --out DIR --sheet reads.tsv

A sheet is tab-separated with columns  template  path  [relation]  where relation (optional) is one of
FR, DP, RR, DM, SAME and overrides the automatic choice. Automatic: a read on the opposite strand of the
first read of its template is FR (independent sequence context, full weight); one on the same strand is
SAME (same primer assumed, half weight) unless --distinct-primers is given (then DP, full weight).
Reads are placed on the reference by local alignment (both strands); only reference positions are
reported (insertions relative to the reference are ignored).

Per reference position the joint genotype over the 10 unordered base pairs is computed (multiread.py).
Consensus letter: homozygous base, or the IUPAC code of a heterozygous call; lowercase when GQ < --min-gq;
N where no read is usable. Written: <template>_consensus.fasta, <template>_positions.tsv (every position
with a heterozygous call, low GQ, or a conflict) and summary.json.

Conflict = at least two reads each individually confident (> 0.9) in different genotypes: replicate
disagreement is reported, never settled by majority vote (ROADMAP 3b).
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from . import bayes
from .concordance import map_to_reference
from .multiread import RELATION_WEIGHT, analyze_read, joint_genotype, oriented, prior_over_genotypes, read_likelihood
from .secondary import IUPAC

COLS = ["ref_pos", "ref_base", "call", "out", "gq", "p_het", "reads_used", "conflict", "per_read"]


def load_templates(args):
    t = defaultdict(list)
    for spec in args.template or []:
        name, paths = spec.split("=", 1)
        for p in paths.split(","):
            t[name].append((p.strip(), None))
    if args.sheet:
        with open(args.sheet, encoding="utf-8") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if not row or row[0].startswith("#") or row[0] == "template":
                    continue
                t[row[0]].append((row[1], row[2] if len(row) > 2 and row[2] else None))
    return t


def consensus_template(name, reads, reference, pi=None, min_gq=10, alleles_at=None, distinct_primers=False):
    """reads: [(ReadAnalysis, relation or None)]. Returns (sequence, rows, info)."""
    maps = [map_to_reference(ra.raw["seq"], reference) for ra, _ in reads]
    strands = []
    for mp in maps:
        st = next(iter(mp.values()))[1] if mp else "+"
        strands.append(st)
    first = strands[0]
    rels = []
    for k, (_, rel) in enumerate(reads):
        if rel:
            rels.append(rel)
        elif k == 0:
            rels.append("FR")        # the anchor read: full weight
        elif strands[k] != first:
            rels.append("FR")
        else:
            rels.append("DP" if distinct_primers else "SAME")
    pi = pi if pi is not None else sum(ra.pi for ra, _ in reads) / len(reads)
    seq, rows = [], []
    n_het = n_conflict = n_cov = 0
    for j in range(len(reference)):
        items = []
        for k, (ra, _) in enumerate(reads):
            if j in maps[k]:
                i, st = maps[k][j]
                items.append((ra, i, st, rels[k]))
        if not items:
            seq.append("N" if rows or seq else "")
            continue
        alleles = (alleles_at or {}).get(j + 1)
        res = joint_genotype(items, pi=pi, alleles=alleles)
        if res["reads_used"] == 0:
            seq.append("N")
            continue
        n_cov += 1
        g = res["call"]
        letter = g[0] if g[0] == g[1] else IUPAC[frozenset(g)]
        # conflict: >= 2 reads individually confident in different genotypes
        confident = []
        for ra, i, st, rel in items:
            if ra.pred[i + 1]["bad"]:
                continue
            lik, _ = read_likelihood(ra, i)
            lik = oriented(lik, st)
            prior = prior_over_genotypes(pi, alleles)
            z = sum(lik[x] * prior[x] for x in lik)
            post = {x: lik[x] * prior[x] / z for x in lik}
            top = max(post, key=post.get)
            if post[top] > 0.9:
                confident.append(top)
        conflict = len(set(confident)) > 1
        out = letter.lower() if res["gq"] < min_gq else letter
        seq.append(out)
        het = g[0] != g[1]
        n_het += het
        n_conflict += conflict
        if het or conflict or res["gq"] < min_gq:
            rows.append(dict(ref_pos=j + 1, ref_base=reference[j], call=g, out=out, gq=res["gq"],
                             p_het=round(res["p_het"], 4), reads_used=res["reads_used"], conflict=int(conflict),
                             per_read=";".join(f"{n}:{p if isinstance(p, str) else round(p['lbf'], 2)}" for n, _, p in res["per_read"])))
    info = dict(reads=[ra.name for ra, _ in reads], relations=rels, strands=strands, pi=round(pi, 4),
                covered=n_cov, het_calls=n_het, conflicts=n_conflict)
    return "".join(seq), rows, info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--template", action="append", help="NAME=read1.ab1,read2.ab1")
    ap.add_argument("--sheet")
    ap.add_argument("--distinct-primers", action="store_true")
    ap.add_argument("--min-gq", type=int, default=10)
    ap.add_argument("--snp", action="append", default=[],
                    help="REFPOS:A/B  known biallelic site (1-based reference position; alleles on the reference strand)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    reference = str(next(SeqIO.parse(a.reference, "fasta")).seq).upper()
    alleles_at = {}
    for s in a.snp:
        pos, al = s.split(":")
        alleles_at[int(pos)] = tuple(al.split("/"))
    summary = {}
    cache = {}
    for name, specs in load_templates(a).items():
        reads = []
        for path, rel in specs:
            if path not in cache:
                cache[path] = analyze_read(path, name=Path(path).stem)
            reads.append((cache[path], rel))
        seq, rows, info = consensus_template(name, reads, reference, min_gq=a.min_gq, alleles_at=alleles_at,
                                             distinct_primers=a.distinct_primers)
        (out / f"{name}_consensus.fasta").write_text(f">{name} sangerQC joint consensus, {len(reads)} reads\n{seq}\n", encoding="utf-8")
        with open(out / f"{name}_positions.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        summary[name] = info
        print(name, {k: info[k] for k in ("relations", "covered", "het_calls", "conflicts")})
    (out / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
