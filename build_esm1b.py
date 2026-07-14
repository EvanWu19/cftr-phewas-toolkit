"""Build a REAL CFTR ESM1b extract from the ntranoslab esm_variants release.

Source: data/ALL_hum_isoforms_ESM1b_LLR.zip :: content/.../P13569_LLR.csv
        (ESM1b LLR matrix for CFTR canonical isoform; Brandes et al. 2023).
Output: data/esm1b_cftr.csv  (per-variant ESM1b LLR, protein_variant keyed)

The source file is an LLR matrix: columns "<wt> <pos>" (e.g. 'M 1'), rows = mutant
amino acid, cell = log-likelihood ratio (more NEGATIVE = more damaging). We melt it
to long form with a 1-letter protein_variant key (e.g. 'M1A'), matching the
AlphaMissense/EVE/CFTR2 key. Canonical P13569 (per plan R14).
"""
from pathlib import Path
import zipfile
import io
import pandas as pd

PKG = Path(__file__).resolve().parent
ZIP = PKG / "data" / "ALL_hum_isoforms_ESM1b_LLR.zip"
MEMBER = "content/ALL_hum_isoforms_ESM1b_LLR/P13569_LLR.csv"
OUT = PKG / "data" / "esm1b_cftr.csv"
OUT.parent.mkdir(exist_ok=True)

with zipfile.ZipFile(ZIP) as z:
    with z.open(MEMBER) as fh:
        m = pd.read_csv(io.TextIOWrapper(fh, encoding="utf-8"), index_col=0)

rows = []
for col in m.columns:                      # col like "M 1"
    wt, pos = col.split(" ")
    pos = int(pos)
    for mut in m.index:                    # mutant amino acid
        if mut == wt:
            continue
        val = m.at[mut, col]
        if pd.isna(val):
            continue
        rows.append((f"{wt}{pos}{mut}", wt, pos, mut, round(float(val), 4)))

df = pd.DataFrame(rows, columns=["protein_variant", "wt_aa", "position", "mt_aa", "esm1b_score"])
df = df.sort_values("position").reset_index(drop=True)
df["source"] = "REAL"
df.to_csv(OUT, index=False)
print("REAL ESM1b CFTR variants written:", len(df), "->", OUT.relative_to(PKG))
print("residues covered:", df["position"].nunique())
print("esm1b_score range:", df["esm1b_score"].min(), "->", df["esm1b_score"].max())
print("pathogenic (<= -7.5):", int((df["esm1b_score"] <= -7.5).sum()))
