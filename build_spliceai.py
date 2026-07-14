"""Build a REAL CFTR SpliceAI extract from the precomputed Illumina VCF.

Source: 05_analysis/spliceai_scores.masked.snv.hg38.vcf.gz  (Illumina precomputed
        SpliceAI v1.3, GRCh38, masked SNV; ~28.6 GB, from BaseSpace share
        https://basespace.illumina.com/s/otSPW8hnhaZR — external, never committed).
Output: data/spliceai_cftr_2021_v1.3.csv

pysam won't build on Windows, so we parse the .tbi index by hand and use
Bio.bgzf to SEEK to the CFTR region — reading ~200 kb, not 28.6 GB. Parses the
SpliceAI INFO field and computes DS_max. Keyed by genomic coordinate.

LICENSE: SpliceAI scores are CC BY-NC 4.0 — source VCF stays external; publishing
the derived extract needs attribution + a non-commercial notice (plan R43).
"""
import gzip
import struct
from pathlib import Path
from Bio import bgzf
import pandas as pd

PKG = Path(__file__).resolve().parent
VCF = PKG / "data" / "spliceai_scores.masked.snv.hg38.vcf.gz"
TBI = Path(str(VCF) + ".tbi")
OUT = PKG / "data" / "spliceai_cftr_2021_v1.3.csv"
OUT.parent.mkdir(exist_ok=True)

# CFTR, GRCh38 (safe superset of chr7:117,480,025-117,668,665)
REGION_START, REGION_END = 117_470_000, 117_670_000
LIDX_SHIFT = 14   # tabix linear-index window = 16 kb


def tbi_start_offset(tbi_path):
    """Parse the .tbi and return (contig_name, bgzf virtual offset) to start reading."""
    raw = gzip.open(tbi_path, "rb").read()
    off = [0]

    def take(fmt):
        sz = struct.calcsize(fmt)
        v = struct.unpack_from(fmt, raw, off[0]); off[0] += sz
        return v

    (magic,) = take("<4s")
    assert magic == b"TBI\x01", magic
    n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = take("<8i")
    names = raw[off[0]:off[0] + l_nm].split(b"\x00")[:-1]
    off[0] += l_nm
    names = [n.decode() for n in names]
    target = next((names.index(c) for c in ("7", "chr7") if c in names), None)
    if target is None:
        raise SystemExit(f"no chr7 contig in tbi; sample names: {names[:6]}")

    for r in range(n_ref):
        (n_bin,) = take("<i")
        for _ in range(n_bin):
            _bin, n_chunk = take("<Ii")
            off[0] += n_chunk * 16          # skip chunks (2 x uint64)
        (n_intv,) = take("<i")
        intv = struct.unpack_from("<%dQ" % n_intv, raw, off[0]); off[0] += n_intv * 8
        if r == target:
            li = min(REGION_START >> LIDX_SHIFT, n_intv - 1)
            # first non-zero offset at/after the window
            voff = next((intv[i] for i in range(li, n_intv) if intv[i]), 0)
            return names[target], voff
    raise SystemExit("target ref not reached")


contig, voff = tbi_start_offset(TBI)
print(f"contig={contig}  start virtual-offset={voff}")

rows = []
reader = bgzf.BgzfReader(str(VCF), "rb")
reader.seek(voff)
for bline in reader:
    line = bline.decode() if isinstance(bline, bytes) else bline
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if f[0] != contig:
        break
    pos = int(f[1])
    if pos < REGION_START:
        continue
    if pos > REGION_END:
        break
    info = f[7]
    if "SpliceAI=" not in info:
        continue
    field = next(x for x in info.split(";") if x.startswith("SpliceAI="))[len("SpliceAI="):]
    parts = field.split("|")   # ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_...
    if len(parts) < 6 or parts[1] != "CFTR":
        continue
    try:
        ds = [float(parts[i]) for i in (2, 3, 4, 5)]
    except ValueError:
        continue
    rows.append({"chrom": f[0], "pos": pos, "ref": f[3], "alt": f[4], "symbol": parts[1],
                 "DS_AG": ds[0], "DS_AL": ds[1], "DS_DG": ds[2], "DS_DL": ds[3],
                 "spliceai_ds_max": round(max(ds), 4), "source": "REAL"})
reader.close()

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print("REAL SpliceAI CFTR SNVs written:", len(df), "->", OUT.relative_to(PKG))
if len(df):
    print("pos range:", df["pos"].min(), "-", df["pos"].max())
    print("ds_max range:", df["spliceai_ds_max"].min(), "->", df["spliceai_ds_max"].max())
    print("HIGH (>=0.5):", int((df["spliceai_ds_max"] >= 0.5).sum()),
          "| MODERATE (0.2-0.5):", int(((df["spliceai_ds_max"] >= 0.2) & (df["spliceai_ds_max"] < 0.5)).sum()))
