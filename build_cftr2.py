"""Build a REAL, committable CFTR2 extract from the official CFTR2 release xlsx.

Source: docs/CFTR2_30January2026.xlsx  (public CFTR2 variant list, cftr2.org).
Output: data/cftr2_2026-01-30.csv

Keeps only variant-level, public columns (no patient-level data). Derives a
1-letter missense key (`protein_variant`, e.g. 'G551D') from the protein name so
CFTR2 class can be joined onto the AlphaMissense/gnomAD missense tables.
Re-run this only when a newer CFTR2 release is dropped into docs/.
"""
from pathlib import Path
import re
import openpyxl
import pandas as pd

PKG = Path(__file__).resolve().parent
XLSX = PKG.parents[1] / "docs" / "CFTR2_30January2026.xlsx"
OUT = PKG / "data" / "cftr2_2026-01-30.csv"
OUT.parent.mkdir(exist_ok=True)

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
# p.Gly551Asp  -> G551D   (simple single-residue missense only)
MIS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def missense_key(protein_name: str):
    if not protein_name:
        return ""
    m = MIS.match(protein_name.strip())
    if not m:
        return ""
    a, pos, b = m.group(1), m.group(2), m.group(3)
    if a in AA3TO1 and b in AA3TO1:
        return f"{AA3TO1[a]}{pos}{AA3TO1[b]}"
    return ""


wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
ws = wb["CFTR2 variants by legacy name"]
rows = []
for r in ws.iter_rows(min_row=13, values_only=True):
    if r[0] is None:
        continue
    legacy, protein, cdna, alt, alleles, af, prev, cur, changed = r[:9]
    rows.append({
        "protein_variant": missense_key(protein or ""),
        "legacy_name": legacy,
        "protein_name": protein,
        "cdna_name": cdna,
        "cftr2_alleles": alleles,
        "cftr2_af": af,
        "cftr2_class": cur,
    })
df = pd.DataFrame(rows)

# merge GRCh38 genomic coordinates from sheet 2 (on cDNA name)
ws2 = wb["Genomic coordinates"]
g = pd.DataFrame(ws2.iter_rows(min_row=2, values_only=True),
                 columns=[c for c in next(ws2.iter_rows(min_row=1, max_row=1, values_only=True))])
gcols = {"Variant cDNA name": "cdna_name", "grch38_chr": "grch38_chr",
         "grch38_pos": "grch38_pos", "grch38_ref": "grch38_ref", "grch38_alt": "grch38_alt"}
g = g[[c for c in gcols]].rename(columns=gcols).drop_duplicates("cdna_name")
df = df.merge(g, on="cdna_name", how="left")

df.to_csv(OUT, index=False)
print("wrote", OUT.relative_to(PKG), "rows:", len(df))
print("with a 1-letter missense key:", (df["protein_variant"] != "").sum())
print("class counts:\n", df["cftr2_class"].value_counts(dropna=False).to_string())
