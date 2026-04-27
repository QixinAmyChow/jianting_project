"""
Full pipeline: reorganize GEO flat files → analysis → figure generation.
"""
import os, sys, re, shutil
sys.path.insert(0, os.path.dirname(__file__))

OUT     = "bm_analysis_out"
FIG_DIR = os.path.join(OUT, "figures")
DAT_DIR = os.path.join(OUT, "data")
RAW_DIR = os.path.join(OUT, "raw_geo")
EXT_DIR = os.path.join(RAW_DIR, "extracted")
SAMP_DIR = os.path.join(RAW_DIR, "samples")

for d in [FIG_DIR, DAT_DIR, SAMP_DIR]:
    os.makedirs(d, exist_ok=True)

# ── Step 1: Reorganise flat GEO files into per-sample subdirs ────────────────
print("Step 1: Reorganising flat GEO files …")

files = os.listdir(EXT_DIR)
# Map: prefix → {barcodes, features, matrix}
sample_map = {}
for f in files:
    if f.endswith("_barcodes.tsv.gz"):
        pfx = f[:-len("_barcodes.tsv.gz")]
        sample_map.setdefault(pfx, {})["barcodes"] = f
    elif f.endswith("_features.tsv.gz"):
        pfx = f[:-len("_features.tsv.gz")]
        sample_map.setdefault(pfx, {})["features"] = f
    elif f.endswith("_matrix.mtx.gz"):
        pfx = f[:-len("_matrix.mtx.gz")]
        sample_map.setdefault(pfx, {})["matrix"] = f

print(f"  Found {len(sample_map)} samples.")

for pfx, fmap in sample_map.items():
    sdir = os.path.join(SAMP_DIR, pfx)
    os.makedirs(sdir, exist_ok=True)
    rename = {"barcodes": "barcodes.tsv.gz",
              "features": "features.tsv.gz",
              "matrix":   "matrix.mtx.gz"}
    for key, std_name in rename.items():
        src = os.path.join(EXT_DIR, fmap[key])
        dst = os.path.join(sdir, std_name)
        if not os.path.exists(dst):
            os.symlink(os.path.abspath(src), dst)

print("  Done reorganising.")

# ── Step 2: Parse metadata ───────────────────────────────────────────────────
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
sc.settings.n_jobs = 4
sc.settings.figdir = FIG_DIR

from bm_macrophage_analysis import (
    parse_geo_metadata, preprocess,
    annotate_and_extract_macrophages, subcluster_macrophages,
    harmonise_metadata,
    fig_umap_overview, fig_cancer_type_bar, fig_er_sp1,
    fig_metastasis_enrichment, fig_treatment_response,
    fig_clinical_metadata, fig_er_sp1_coexpr,
    MACRO_MARKERS
)

meta_gz = os.path.join(RAW_DIR, "GSE266330_series_matrix.txt.gz")
print("Step 2: Parsing metadata …")
meta_df = parse_geo_metadata(meta_gz)
print(f"  Metadata shape: {meta_df.shape}")

# ── Step 3: Load count matrices ──────────────────────────────────────────────
processed_path = os.path.join(DAT_DIR, "mac_subclustered.h5ad")
full_path       = os.path.join(DAT_DIR, "full_preprocessed.h5ad")

if os.path.exists(processed_path):
    print("Loading cached macrophage AnnData …")
    mac = sc.read_h5ad(processed_path)
else:
    if os.path.exists(full_path):
        print("Loading cached full AnnData …")
        adata = sc.read_h5ad(full_path)
    else:
        print("Step 3: Loading count matrices …")
        adatas = []
        sample_dirs = sorted(
            [os.path.join(SAMP_DIR, d) for d in os.listdir(SAMP_DIR)
             if os.path.isdir(os.path.join(SAMP_DIR, d))]
        )
        print(f"  Loading {len(sample_dirs)} samples …")
        for sdir in sample_dirs:
            sample_id = os.path.basename(sdir)
            try:
                a = sc.read_10x_mtx(sdir, var_names="gene_symbols", cache=True)
                a.var_names_make_unique()
                a.obs["sample_id"] = sample_id
                adatas.append(a)
            except Exception as e:
                print(f"    Skipping {sample_id}: {e}")

        print(f"  Concatenating {len(adatas)} samples …")
        adata = ad.concat(
            adatas,
            label="sample_id",
            keys=[a.obs["sample_id"][0] for a in adatas],
            merge="same"
        )
        adata.obs_names_make_unique()

        # Attach metadata
        if not meta_df.empty:
            meta_sub = meta_df.reindex(adata.obs["sample_id"])
            meta_sub.index = adata.obs.index
            for col in meta_df.columns:
                adata.obs[col] = meta_sub[col].values

        print(f"  Full dataset: {adata.n_obs} cells, {adata.n_vars} genes")

        print("Step 4: Pre-processing …")
        adata = preprocess(adata)
        adata.write_h5ad(full_path)

    print("Step 5: Isolating macrophages …")
    mac = annotate_and_extract_macrophages(adata)

    print("Step 6: Sub-clustering macrophages …")
    mac = subcluster_macrophages(mac)

    print("Step 7: Harmonising metadata …")
    mac.obs = harmonise_metadata(mac.obs)
    mac.write_h5ad(processed_path)
    print(f"  Saved: {processed_path}")

# ── Step 8: Generate figures ──────────────────────────────────────────────────
print("\nGenerating figures …")
fig_umap_overview(mac)
fig_cancer_type_bar(mac)
fig_er_sp1(mac)
fig_metastasis_enrichment(mac)
fig_treatment_response(mac)
fig_clinical_metadata(mac)
top_cluster = fig_er_sp1_coexpr(mac)

print(f"\nAll figures saved to: {FIG_DIR}")
print(f"ER+/SP1-high cluster: {top_cluster}")
print(f"Macrophage sub-clusters: {mac.obs['mac_cluster'].nunique()}")
print("DONE")
