"""Local Streamlit app: upload an .ab1, run v0.3, draw PLOC grid + thresholds, download FASTA + decisions."""
import sys
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sangerqc.plot import figure_strips
import csv
import io

from sangerqc.export import FIELDS, to_sequence
from sangerqc.v0_2 import ALPHA, BETA
from sangerqc.v0_3 import classify

st.set_page_config(page_title="sangerQC", layout="wide")
st.title("sangerQC")
st.caption(
    "Upload a Sanger `.ab1`. Vertical bars are the called-peak scan positions "
    "(PLOC), not a uniform base grid. Horizontal lines are this read’s own "
    "HQ-body amplitude thresholds."
)

uploaded = st.file_uploader("AB1 chromatogram", type=["ab1", "abif"])
c1, c2, c3 = st.columns(3)
p1, p2 = st.columns(2)
fwd = p1.text_input("forward primer (optional)", "").strip().upper() or None
rev = p2.text_input("reverse primer (optional)", "").strip().upper() or None
alpha = c1.number_input("alpha (relative tests off below this × HQ)", 0.05, 1.0, ALPHA, 0.01)
beta = c2.number_input("beta (stop below this × HQ)", 0.01, 1.0, BETA, 0.01)
bases = c3.slider("bases per strip", 40, 160, 80, 10)

if uploaded is None:
    st.info("Drop an `.ab1` to classify. Nothing is uploaded anywhere — this runs locally.")
    st.stop()

tmp = Path(st.session_state.get("_tmp") or "")
cache = Path(".streamlit_tmp")
cache.mkdir(exist_ok=True)
path = cache / uploaded.name
path.write_bytes(uploaded.getvalue())

pred, raw = classify(str(path), alpha=alpha, beta=beta, forward=fwd, reverse=rev)
n = len(raw["seq"])
n_bad = sum(1 for v in pred.values() if v["bad"])
n_low = sum(1 for v in pred.values() if v["reason"] == "low_amp_keep")
one = pred[1]
st.write(
    f"{n} bases · HQ body {one['hq_body']:.0f} RFU on island {one['hq_island']} · "
    f"{n_bad} flagged ({100*n_bad/n:.1f}%) · {n_low} low-amp kept"
)

diag = raw["diagnostics"]
if diag["low_snr"]:
    st.warning(f"Low SNR: HQ {diag['hq']:.0f} RFU < 2 x noise floor {diag['noise_floor']:.0f}; "
               "the alpha/beta gates do not mean much on this read.")
if diag.get("primer_notes"):
    st.write("Primers: " + "; ".join(diag["primer_notes"]))

fig = figure_strips(raw, pred, bases_per_strip=bases, alpha=alpha, beta=beta)
st.pyplot(fig, clear_figure=True)

seq, rows, span = to_sequence(pred, raw)
st.subheader("v0.3 sequence")
st.caption("IUPAC = co-located second peak >= 0.33; lowercase = flagged or 5' low-amplitude start; "
           "single read, no reference: ends are not checked against a reference here (use the CLI with --reference).")
st.code(f">{uploaded.name} kept {span}\n{seq}")
buf = io.StringIO()
w = csv.DictWriter(buf, fieldnames=FIELDS, delimiter="\t")
w.writeheader()
w.writerows(rows)
st.download_button("decisions.tsv", buf.getvalue(), file_name=f"{Path(uploaded.name).stem}_decisions.tsv")
