"""Build a REAL CFTR SpliceAI extract from the precomputed Illumina VCFs.

Sources (external, never committed — from the BaseSpace share
https://basespace.illumina.com/s/otSPW8hnhaZR, Illumina precomputed SpliceAI v1.3,
GRCh38). Illumina ships a 2x2 of {masked, raw} x {snv, indel}; this builder reads
**one SNV file and one indel file** and concatenates them:

    data/spliceai_scores.masked.snv.hg38.vcf.gz     ~28.6 GB   -> SNVs   (preferred)
    data/spliceai_scores.raw.snv.hg38.vcf.gz        ~28.8 GB   -> SNVs   (fallback)
    data/spliceai_scores.masked.indel.hg38.vcf.gz   (often a broken 0-byte download)
    data/spliceai_scores.raw.indel.hg38.vcf.gz      ~69.3 GB   -> indels (fallback)

Output: data/spliceai_cftr_2021_v1.3.csv  (~2.08 M rows: ~566 k SNVs + ~1.51 M indels)

⚠ MIXED masked/raw IS LIKELY, AND IS RECORDED. `masked.indel` is very often a 0-byte
download, so the usual result is **masked SNVs + raw indels**. Masked zeroes out
biologically implausible directions (a gain at an existing site, a loss where no site
exists); raw does not, so raw scores are >= masked for the same variant. Every row
therefore carries `score_type` ('masked' / 'raw') and `variant_class` ('snv' / 'indel'),
and the builder warns when the two are mixed. Measured on the CFTR region, the seam is
small but not nothing: masked and raw agree exactly on 95.8% of the 566,106 SNVs
(mean |diff| 0.0019), yet **113 SNVs cross the 0.5 HIGH cut** between the two. Do not
compare a masked score against a raw one and call the difference biology.

WHY INDELS MATTER HERE: the missense predictors cannot score an indel at all, so before
this file was included every CFTR2 deletion/insertion was invisible to SpliceAI —
including F508del, which SpliceAI scores DS_max = 0.01 (correctly: it is a folding
defect, not a splice one). See notebook 13.

pysam won't build on Windows, so we parse the .tbi index by hand and use Bio.bgzf to
SEEK to the CFTR region — reading ~200 kb of each file, not 28-69 GB. Parses the SpliceAI
INFO field and computes DS_max. Keyed by genomic coordinate.

LICENSE: SpliceAI scores are CC BY-NC 4.0 — source VCFs stay external; publishing the
derived extract needs attribution + a non-commercial notice (plan R43).
"""
import gzip
import struct
from pathlib import Path
from Bio import bgzf
import pandas as pd

PKG = Path(__file__).resolve().parent
DATA = PKG / "data"
OUT = DATA / "spliceai_cftr_2021_v1.3.csv"
OUT.parent.mkdir(exist_ok=True)

# CFTR, GRCh38 (safe superset of chr7:117,480,025-117,668,665)
REGION_START, REGION_END = 117_470_000, 117_670_000
LIDX_SHIFT = 14   # tabix linear-index window = 16 kb


def pick(*names):
    """First of `names` that exists in data/ and is non-empty.

    The size check is the point: a 0-byte `masked.indel` is a *common* failed download,
    and silently reading it would produce an empty indel set that looks like "SpliceAI
    cannot score indels" rather than "the file never downloaded".
    """
    for n in names:
        p = DATA / n
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


SNV_VCF = pick("spliceai_scores.masked.snv.hg38.vcf.gz",
               "spliceai_scores.raw.snv.hg38.vcf.gz")
INDEL_VCF = pick("spliceai_scores.masked.indel.hg38.vcf.gz",
                 "spliceai_scores.raw.indel.hg38.vcf.gz")


def tbi_start_offset(tbi_path, region_start=REGION_START):
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
            li = min(region_start >> LIDX_SHIFT, n_intv - 1)
            # first non-zero offset at/after the window
            voff = next((intv[i] for i in range(li, n_intv) if intv[i]), 0)
            return names[target], voff
    raise SystemExit("target ref not reached")


def scan(vcf_path):
    """Seek to the CFTR region in one precomputed VCF and yield parsed score rows."""
    score_type = "masked" if ".masked." in vcf_path.name else "raw"
    variant_class = "indel" if ".indel." in vcf_path.name else "snv"
    contig, voff = tbi_start_offset(Path(str(vcf_path) + ".tbi"))
    print(f"  {vcf_path.name}  (contig={contig}, {score_type}/{variant_class})")

    rows = []
    reader = bgzf.BgzfReader(str(vcf_path), "rb")
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
                     "spliceai_ds_max": round(max(ds), 4),
                     "variant_class": variant_class, "score_type": score_type,
                     "source": "REAL"})
    reader.close()
    print(f"    -> {len(rows):,} CFTR {variant_class} records")
    return rows


rows = []
print("reading:")
if SNV_VCF is None:
    print("  ! no SNV VCF found in data/ — extract will have no SNVs")
else:
    rows += scan(SNV_VCF)
if INDEL_VCF is None:
    print("  ! no non-empty indel VCF found in data/ — extract will have no indels.\n"
          "    NB: a 0-byte spliceai_scores.masked.indel.hg38.vcf.gz is a common failed\n"
          "    download; re-fetch it, or use spliceai_scores.raw.indel.hg38.vcf.gz.")
else:
    rows += scan(INDEL_VCF)

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"\nREAL SpliceAI CFTR records written: {len(df):,} -> {OUT.relative_to(PKG)}")
if len(df):
    print(df.groupby(["variant_class", "score_type"]).size().to_string())
    if df["score_type"].nunique() > 1:
        print("\n  ⚠ MIXED masked/raw extract — see the module docstring. Rows carry\n"
              "    `score_type`; do not compare a masked score against a raw one.")
    print("\npos range:", df["pos"].min(), "-", df["pos"].max())
    print("ds_max range:", df["spliceai_ds_max"].min(), "->", df["spliceai_ds_max"].max())
    print("HIGH (>=0.5):", int((df["spliceai_ds_max"] >= 0.5).sum()),
          "| MODERATE (0.2-0.5):", int(((df["spliceai_ds_max"] >= 0.2) & (df["spliceai_ds_max"] < 0.5)).sum()))
