# Changes after v0.3

## 2026-10-07
- **Fix:** primary channel at a call = tallest channel with a real apex on the call, not the tallest value in the
  slot (a neighbour's tail could win and turn a single base into a "double"). Python and JS. Found on the
  allergology plates; synthetic benchmarks unchanged, cyclops site at reference 269 moved from strong to against.
- **Fix:** miss rate of a real second peak floored at 2 % (default 5 %); was 0.1 %, fitted where it never happens.
- **New:** per-read two-base model (EM, shrunk to the global fit); genotype probabilities, Phred likelihoods and
  GQ per call (`gt`, `gq`, `pl`, `p_mixed`) after reading Clair3's source (BAYES.md).
- **New:** `multiread.py`, `multicli.py`: joint genotype and consensus from several reads of one template;
  `validation/eval_multiread.py` (needs the plates, which are not in the repo).
- **New:** `validation/eval_bayes*.py`, `bench_hetindel.py`; heterozygous-indel decoder prototype (`hetindel.py`).
- Viewer: classifier ported to JS (parity checked on 9 reads).
