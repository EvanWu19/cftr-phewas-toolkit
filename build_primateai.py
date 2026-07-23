"""Build a REAL CFTR PrimateAI extract from the dbNSFP v5.0a parquet.

Source: data/dbNSFP5.0a_variant.clin_var_re-annot_pdb_variants_plddt_rsa.parquet
        (dbNSFP v5.0a, ClinVar-re-annotated subset; bundles PrimateAI_score etc.)
Output: data/primateai_cftr.csv

COVERAGE: this dbNSFP parquet is the ClinVar-re-annotated subset, so PrimateAI
covers ~1,976 CFTR missense variants in the extract (of which ~1,271 overlap gnomAD's
observed 2,466), NOT genome-wide saturation. Honestly labelled as such.

PROVENANCE: the original PrimateAI (Sundaram 2018 — semi-supervised on common human +
6-primate-species variants as a benign proxy, ~2017 data) shipped as precomputed
genome-wide scores from Illumina (BaseSpace/Zenodo), which are gated/large. dbNSFP
v5.0a bundles the PrimateAI_score column, so we pull it from that (smaller, accessible)
parquet instead — trading genome-wide coverage for convenience. dbNSFP is CC BY-NC-ND
(non-commercial); PrimateAI non-commercial.
"""
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd

PKG = Path(__file__).resolve().parent
PARQ = PKG / "data" / "dbNSFP5.0a_variant.clin_var_re-annot_pdb_variants_plddt_rsa.parquet"
OUT = PKG / "data" / "primateai_cftr.csv"
OUT.parent.mkdir(exist_ok=True)
START, END = 117_470_000, 117_670_000
AA1 = set("ACDEFGHIKLMNPQRSTVWY")

cols = ['#chr', 'pos(1-based)', 'ref', 'alt', 'aaref', 'aaalt', 'aapos',
        'PrimateAI_score', 'PrimateAI_pred']
df = pq.read_table(PARQ, columns=cols, filters=[('#chr', '==', '7')]).to_pandas()
df['pos'] = df['pos(1-based)'].astype('int64')
cf = df[(df['pos'] >= START) & (df['pos'] <= END)].copy()

def num(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None

cf['primate_ai_score'] = cf['PrimateAI_score'].map(num)
cf = cf.dropna(subset=['primate_ai_score'])
# missense only + 1-letter protein_variant key
cf = cf[cf['aaref'].isin(AA1) & cf['aaalt'].isin(AA1)]
cf['protein_variant'] = cf['aaref'] + cf['aapos'].astype(str) + cf['aaalt']
cf['chrom'] = cf['#chr'].astype(str)
out = cf[['chrom', 'pos', 'ref', 'alt', 'protein_variant',
          'primate_ai_score', 'PrimateAI_pred']].rename(columns={'PrimateAI_pred': 'primate_ai_pred'})
out = out.drop_duplicates('protein_variant').sort_values('pos').reset_index(drop=True)
out['source'] = 'REAL'
out.to_csv(OUT, index=False)
print("REAL PrimateAI CFTR variants written:", len(out), "->", OUT.relative_to(PKG))
print("score range:", out['primate_ai_score'].min(), "->", out['primate_ai_score'].max())
print("pathogenic (>= 0.803):", int((out['primate_ai_score'] >= 0.803).sum()))
