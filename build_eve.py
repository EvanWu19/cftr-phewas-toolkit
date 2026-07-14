"""Build a REAL CFTR EVE extract from the official EVE release zip.

Source: 05_analysis/EVE_all_data.zip  (EVE, Frazer et al. 2021 release, Aug 2021),
        member file  variant_files/CFTR_HUMAN.csv.
Output: data/eve_cftr_2021-08.csv  (small, per-variant EVE scores for CFTR)

Derives a 1-letter ``protein_variant`` key (e.g. 'M1A') = wt_aa + position + mt_aa,
matching the AlphaMissense / gnomAD / CFTR2 missense key, and keeps the EVE score
(EVE_scores_ASM, 0-1, higher = more pathogenic) plus EVE's own class call.
Re-run only when a newer EVE release is dropped in.

NOTE: EVE data is redistributed under EVE's (evemodel.org) terms — confirm the
license before committing this extract to a public repo. The 9.6 GB source zip
stays external (never committed).
"""
from pathlib import Path
import zipfile
import io
import pandas as pd

PKG = Path(__file__).resolve().parent
ZIP = PKG / "data" / "EVE_all_data.zip"
MEMBER = "variant_files/CFTR_HUMAN.csv"
OUT = PKG / "data" / "eve_cftr_2021-08.csv"
OUT.parent.mkdir(exist_ok=True)

with zipfile.ZipFile(ZIP) as z:
    with z.open(MEMBER) as fh:
        raw = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"),
                          usecols=["wt_aa", "position", "mt_aa",
                                   "EVE_scores_ASM", "EVE_classes_75_pct_retained_ASM"],
                          low_memory=False)

print("rows in CFTR_HUMAN.csv:", len(raw))
# keep only rows EVE actually scored
df = raw.dropna(subset=["EVE_scores_ASM"]).copy()
df["protein_variant"] = df["wt_aa"].astype(str) + df["position"].astype(int).astype(str) + df["mt_aa"].astype(str)
df = df.rename(columns={"EVE_scores_ASM": "eve_score",
                        "EVE_classes_75_pct_retained_ASM": "eve_class"})
df["source"] = "REAL"
out = df[["protein_variant", "wt_aa", "position", "mt_aa", "eve_score", "eve_class", "source"]]
out = out.sort_values("position").reset_index(drop=True)
out.to_csv(OUT, index=False)

print("scored EVE variants written:", len(out), "->", OUT.relative_to(PKG))
print("eve_score range:", round(out["eve_score"].min(), 4), "->", round(out["eve_score"].max(), 4))
print("unique residues covered:", out["position"].nunique())
print("class counts:\n", out["eve_class"].value_counts(dropna=False).to_string())
