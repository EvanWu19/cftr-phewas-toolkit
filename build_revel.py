"""Build a REAL CFTR REVEL extract from the genome-wide REVEL v1.3 table.

Source: data/revel-v1.3_all_chromosomes.zip :: revel_with_transcript_ids
        (one 6.5 GB CSV; Ioannidis et al. 2016; free for non-commercial use).
Output: data/revel_cftr_v1.3.csv

NOTE ON THE ODD FILENAME: inside the zip the data is a single ~6.5 GB CSV literally
named `revel_with_transcript_ids` — with **no `.csv` extension**. That is just how
the REVEL authors package it (sites.google.com/site/revelgenomics); it is not a
corrupt or unknown file. It has a normal CSV header/rows.

REVEL is keyed by genomic coordinate (no protein position in the file). We stream
the CSV (chr-grouped), keep chr7 rows in the CFTR GRCh38 window, and stop once we
leave the chr7 block. Join onto observed variants by (chrom,pos,ref,alt).
Columns: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid.

⚠ REVEL is SUPERVISED (trained ~2015 on HGMD/ClinVar-lineage pathogenic + ESP
neutral), so benchmarking it against ClinVar/CFTR2 labels added after that risks
temporal leakage — see notebook 12.
"""
from pathlib import Path
import zipfile
import io
import csv
import pandas as pd

PKG = Path(__file__).resolve().parent
ZIP = PKG / "data" / "revel-v1.3_all_chromosomes.zip"
MEMBER = "revel_with_transcript_ids"
OUT = PKG / "data" / "revel_cftr_v1.3.csv"
OUT.parent.mkdir(exist_ok=True)
START, END = 117_470_000, 117_670_000

rows = []
seen7 = False
scanned = 0
with zipfile.ZipFile(ZIP) as z:
    with z.open(MEMBER) as fh:
        reader = csv.reader(io.TextIOWrapper(fh, encoding="utf-8", newline=""))
        header = next(reader)
        for f in reader:
            scanned += 1
            chrom = f[0]
            if chrom == "7":
                seen7 = True
                g = f[2]                       # grch38_pos
                if g == "." or g == "":
                    continue
                p = int(g)
                if START <= p <= END:
                    rows.append((chrom, p, f[3], f[4], f[5], f[6], round(float(f[7]), 4)))
            elif seen7:
                break                          # left the contiguous chr7 block

df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "aaref", "aaalt", "revel_score"])
df = df.sort_values("pos").reset_index(drop=True)
df["source"] = "REAL"
df.to_csv(OUT, index=False)
print(f"scanned {scanned:,} REVEL rows to reach/pass chr7")
print("REAL REVEL CFTR variants written:", len(df), "->", OUT.relative_to(PKG))
if len(df):
    print("pos range:", df["pos"].min(), "-", df["pos"].max())
    print("revel_score range:", df["revel_score"].min(), "->", df["revel_score"].max())
