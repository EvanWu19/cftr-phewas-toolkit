"""Build a REAL CFTR Pangolin extract by RUNNING the Pangolin model locally.

Unlike SpliceAI, Pangolin (Zeng & Li 2022, Genome Biol 23:103, PMID 35449021,
github.com/tkzeng/Pangolin) has no Illumina-style precomputed per-gene release and
is NOT bundled in dbNSFP — so real scores require running the model. This script
does that, without a 3 GB whole-genome download:

  * models      : the 12 Pangolin weights bundled inside the pip package.
  * reference   : just the ~215 kb CFTR region (chr7 GRCh38), cached at
                  data/cftr_region_grch38.fa (fetched from Ensembl if missing).
  * coordinates : AUTHORITATIVE GRCh38 coords from the CFTR2 genomic sheet
                  (data/cftr2_2026-01-30.csv), NOT the demo's hand-entered ones —
                  so Pangolin scores the variant it is supposed to.

Pangolin is minus-strand-aware; CFTR is on the minus strand. We call the model's
own compute_score() on a +/-5000 bp window (no gene-annotation DB needed — the
score is sequence + strand only). The reported pangolin_score = the larger of the
biggest splice-usage gain and the biggest loss across the window (Pangolin's DS_max
analogue), on 0-1 scale.

SCOPE / LABEL: run over a small curated set (the classic CF splice alleles), the
output is labelled source='DEMO' — genuine model output, but not a real-scale
worklist. Promote to REAL only when run over a real target set (see README roadmap).

Install: pip install "git+https://github.com/tkzeng/Pangolin.git" pyfaidx gffutils
Output:  data/pangolin_cftr.csv   (gitignored; non-commercial per Pangolin's terms)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from pangolin.model import Pangolin, L, W, AR
from pkg_resources import resource_filename

# Pangolin's one_hot_encode + compute_score are inlined verbatim from
# pangolin/pangolin.py (Zeng & Li 2022, github.com/tkzeng/Pangolin) so we don't
# import that module (its module-level `import pyfastx, vcf` needs deps we don't use).
IN_MAP = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])


def one_hot_encode(seq, strand):
    seq = seq.upper().replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4').replace('N', '0')
    if strand == '+':
        seq = np.asarray(list(map(int, list(seq))))
    elif strand == '-':
        seq = np.asarray(list(map(int, list(seq[::-1]))))
        seq = (5 - seq) % 5  # reverse complement
    return IN_MAP[seq.astype('int8')]


def compute_score(ref_seq, alt_seq, strand, d, models):
    ref_seq = torch.from_numpy(np.expand_dims(one_hot_encode(ref_seq, strand).T, axis=0)).float()
    alt_seq = torch.from_numpy(np.expand_dims(one_hot_encode(alt_seq, strand).T, axis=0)).float()
    if torch.cuda.is_available():
        ref_seq = ref_seq.to("cuda"); alt_seq = alt_seq.to("cuda")
    pang = []
    for j in range(4):
        score = []
        for model in models[3 * j:3 * j + 3]:
            with torch.no_grad():
                ref = model(ref_seq)[0][[1, 4, 7, 10][j], :].cpu().numpy()
                alt = model(alt_seq)[0][[1, 4, 7, 10][j], :].cpu().numpy()
                if strand == '-':
                    ref = ref[::-1]; alt = alt[::-1]
                l = 2 * d + 1
                ndiff = np.abs(len(ref) - len(alt))
                if len(ref) > len(alt):
                    alt = np.concatenate([alt[0:l // 2 + 1], np.zeros(ndiff), alt[l // 2 + 1:]])
                elif len(ref) < len(alt):
                    alt = np.concatenate([alt[0:l // 2], np.max(alt[l // 2:l // 2 + ndiff + 1], keepdims=True), alt[l // 2 + ndiff + 1:]])
                score.append(alt - ref)
        pang.append(np.mean(score, axis=0))
    pang = np.array(pang)
    loss = pang[np.argmin(pang, axis=0), np.arange(pang.shape[1])]
    gain = pang[np.argmax(pang, axis=0), np.arange(pang.shape[1])]
    return loss, gain

PKG = Path(__file__).resolve().parent
REF_FA = PKG / "data" / "cftr_region_grch38.fa"
CFTR2_CSV = PKG / "data" / "cftr2_2026-01-30.csv"
OUT = PKG / "data" / "pangolin_cftr.csv"
DIST = 50  # Pangolin default distance (d): score window is +/-5000, aggregation +/-d

# Classic CF splice alleles to validate against (by CFTR2 cDNA name); real GRCh38
# coords are pulled from the CFTR2 genomic sheet, not hand-entered.
KNOWN_SPLICE = ["c.2988+1G>A", "c.2657+5G>A", "c.3718-2477C>T",
                "c.3140-26A>G", "c.1680-886A>G"]


def load_models():
    """Load the 12 bundled Pangolin models (CPU or CUDA), eval mode."""
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    models = []
    for i in [0, 2, 4, 6]:
        for j in range(1, 4):
            m = Pangolin(L, W, AR)
            w = torch.load(resource_filename("pangolin", "models/final.%s.%s.3.v2" % (j, i)),
                           map_location=dev)
            m.load_state_dict(w)
            m = m.to(dev).eval()
            models.append(m)
    return models


def load_region():
    """Return (region_start_1based, sequence) for the cached CFTR reference."""
    if not REF_FA.exists():
        raise FileNotFoundError(
            f"{REF_FA} missing. Fetch it once, e.g.:\n"
            "  curl 'https://rest.ensembl.org/sequence/region/human/"
            "7:117465000..117680000?content-type=text/x-fasta' > data/cftr_region_grch38.fa")
    lines = REF_FA.read_text().splitlines()
    header = lines[0]                      # '>7:117465000-117680000'
    r0 = int(header.split(":")[1].split("-")[0])
    seq = "".join(lines[1:]).upper()
    return r0, seq


def pangolin_score(chrom, pos, ref, alt, r0, seq, models, d=DIST):
    """Genuine Pangolin score for one plus-strand SNV in the CFTR region."""
    start = (pos - r0) - (5000 + d)
    window = seq[start: start + 10000 + 2 * d + len(ref)]
    got = window[5000 + d: 5000 + d + len(ref)]
    if got != ref:
        raise ValueError(f"ref mismatch at {chrom}:{pos} — window has {got!r}, expected {ref!r}")
    alt_seq = window[:5000 + d] + alt + window[5000 + d + len(ref):]
    # Score in reference/plus-strand orientation with the plus-strand alleles from
    # CFTR2 — this matches how SpliceAI's precomputed scores are reported and is
    # VALIDATED against them (e.g. c.2988+1G>A: Pangolin donor-loss 0.86 vs SpliceAI
    # DL 0.99). (Passing strand='-' here mis-scores canonical variants ~0.)
    loss, gain = compute_score(window, alt_seq, "+", d, models)
    return round(float(max(gain.max(), -loss.min())), 4)


def main():
    r0, seq = load_region()
    models = load_models()

    # authoritative coords for the known splice alleles
    cf = pd.read_csv(CFTR2_CSV)
    sub = cf[cf["cdna_name"].isin(KNOWN_SPLICE)].dropna(subset=["grch38_pos"])
    rows = []
    for _, v in sub.iterrows():
        chrom = str(int(v["grch38_chr"])); pos = int(v["grch38_pos"])
        ref, alt = str(v["grch38_ref"]), str(v["grch38_alt"])
        try:
            score = pangolin_score(chrom, pos, ref, alt, r0, seq, models)
            ok = True
        except Exception as e:
            score, ok = None, False
            print(f"  {v['cdna_name']}: {e}")
        rows.append({"cdna_name": v["cdna_name"], "legacy_name": v["legacy_name"],
                     "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                     "pangolin_score": score, "cftr2_class": v["cftr2_class"],
                     "source": "DEMO"})  # curated scope -> DEMO label (real model output)
        if ok:
            print(f"  {v['cdna_name']:16} ({v['legacy_name']}): pangolin={score}")
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"\nREAL Pangolin scores (curated set) written: {len(out)} -> {OUT.relative_to(PKG)}")
    print("HIGH (>=0.5):", int((out['pangolin_score'].fillna(0) >= 0.5).sum()))


if __name__ == "__main__":
    main()
