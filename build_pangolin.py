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

Pangolin is strand-aware; CFTR is on the PLUS strand, so we score in plus-strand
orientation (strand='+', see pangolin_score()). We call the model's
own compute_score() on a +/-5000 bp window (no gene-annotation DB needed — the
score is sequence + strand only). The reported pangolin_score = the larger of the
biggest splice-usage gain and the biggest loss across the window (Pangolin's DS_max
analogue), on 0-1 scale.

SCOPE / LABEL — two modes, and the label follows the scope:

  --scope cftr2    (default)  every CFTR2 variant carrying GRCh38 coordinates
                              (~1,893 of 2,097), SNVs *and* indels → source='REAL'.
                              This is a real-scale worklist, so it is labelled REAL.
  --scope curated             just the 5 classic CF splice alleles → source='DEMO'.
                              Genuine model output, but its *coverage* is a teaching
                              subset, not a worklist — hence the DEMO label.

Pangolin scores indels as well as SNVs (compute_score pads the delta track by the
length difference), which is why the cftr2 scope reaches variant classes the
precomputed masked-SNV SpliceAI release cannot. Events larger than MAX_EVENT bp are
skipped: the ±50 bp aggregation window cannot say anything meaningful about a
multi-kb deletion, and padding one would be a fabricated number.

Every target is accounted for — rows that cannot be scored are still written, with
pangolin_score empty and `skip_reason` saying why, so coverage is auditable rather
than silently short.

Install: pip install "git+https://github.com/tkzeng/Pangolin.git" pyfaidx gffutils
Output:  data/pangolin_cftr.csv   (gitignored; non-commercial per Pangolin's terms)
Runtime: ~4 min for the full CFTR2 scope on a CUDA GPU; CPU is many times slower.
"""
import argparse
import re
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
# Largest ref/alt event we will score. Pangolin aggregates over +/-DIST bp, so a
# multi-kb deletion (CFTR2 has one, the 24 kb 'Deletion of Intron 1') is outside what
# the score means — skipped with a reason rather than reported as a number.
MAX_EVENT = 100
_ACGT = re.compile(r"^[ACGT]+$")

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
    """Genuine Pangolin score for one plus-strand SNV or small indel in the CFTR region."""
    start = (pos - r0) - (5000 + d)
    end = start + 10000 + 2 * d + len(ref)
    # Negative/over-long slices would silently return a short (wrong) window, so the
    # bounds are checked rather than left to Python's forgiving slice semantics.
    if start < 0 or end > len(seq):
        raise ValueError(f"outside the cached reference window (needs {chrom}:{pos-5050}-{pos+5050})")
    window = seq[start:end]
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


def targets(scope: str) -> pd.DataFrame:
    """CFTR2 rows to score, with the un-scoreable ones already flagged.

    Returns every requested variant — including the ones we cannot score — so the
    caller can write a complete, auditable table instead of a silently short one.
    """
    cf = pd.read_csv(CFTR2_CSV)
    if scope == "curated":
        cf = cf[cf["cdna_name"].isin(KNOWN_SPLICE)]
    cf = cf.copy()
    ref, alt = cf["grch38_ref"].astype(str), cf["grch38_alt"].astype(str)
    cf["skip_reason"] = None
    cf.loc[cf["grch38_pos"].isna(), "skip_reason"] = "no GRCh38 coordinates in CFTR2"
    bad = cf["skip_reason"].isna() & ~(ref.str.match(_ACGT) & alt.str.match(_ACGT))
    cf.loc[bad, "skip_reason"] = "allele is not plain ACGT"
    big = cf["skip_reason"].isna() & ((ref.str.len() > MAX_EVENT) | (alt.str.len() > MAX_EVENT))
    cf.loc[big, "skip_reason"] = f"event larger than {MAX_EVENT} bp"
    return cf


def main(scope: str = "cftr2") -> None:
    r0, seq = load_region()
    models = load_models()
    tgt = targets(scope)
    label = "REAL" if scope == "cftr2" else "DEMO"   # label follows coverage, not model
    print(f"scope={scope} → {len(tgt):,} CFTR2 variants, source label '{label}'")

    rows, done = [], 0
    for _, v in tgt.iterrows():
        score, reason = None, v["skip_reason"]
        chrom = str(int(v["grch38_chr"])) if pd.notna(v.get("grch38_chr")) else "7"
        pos = int(v["grch38_pos"]) if pd.notna(v["grch38_pos"]) else None
        ref, alt = str(v["grch38_ref"]), str(v["grch38_alt"])
        if reason is None:
            try:
                score = pangolin_score(chrom, pos, ref, alt, r0, seq, models)
                done += 1
                if scope == "curated" or done % 200 == 0:
                    print(f"  [{done:5,}] {str(v['cdna_name'])[:24]:24} "
                          f"({str(v['legacy_name'])[:16]:16}) pangolin={score}")
            except Exception as e:                    # ref mismatch / out of window
                reason = str(e).split(" — ")[0] if "—" in str(e) else str(e)
        rows.append({"cdna_name": v["cdna_name"], "legacy_name": v["legacy_name"],
                     "chrom": chrom, "pos": pos, "ref": ref if pos else None,
                     "alt": alt if pos else None, "pangolin_score": score,
                     "cftr2_class": v["cftr2_class"], "source": label,
                     "skip_reason": reason})

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    scored = out["pangolin_score"].notna()
    print(f"\nPangolin ({label}) written: {int(scored.sum()):,} scored / {len(out):,} targets "
          f"-> {OUT.relative_to(PKG)}")
    print(f"  HIGH (>=0.5): {int((out['pangolin_score'].fillna(0) >= 0.5).sum()):,}"
          f" | MODERATE (0.2-0.5): "
          f"{int(((out['pangolin_score'] >= 0.2) & (out['pangolin_score'] < 0.5)).sum()):,}")
    if (~scored).any():
        print("\n  not scored, by reason:")
        print(out.loc[~scored, "skip_reason"].value_counts().to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scope", choices=["cftr2", "curated"], default="cftr2",
                    help="cftr2 = every CFTR2 variant with coords (REAL); "
                         "curated = the 5 classic splice alleles (DEMO)")
    main(ap.parse_args().scope)
