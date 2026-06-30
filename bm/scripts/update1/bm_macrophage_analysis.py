"""
BM Project — Human Data Analysis
Source: Liu et al. 2025, Cell Genomics (GSE266330)
Goal:   Sub-cluster macrophages; identify ER+/SP1-high clusters;
        characterize by cancer origin, metastasis status, treatment response,
        gender, and tissue site.

Usage:
    python3 bm_macrophage_analysis.py

Outputs (saved to ./bm_analysis_out/):
    figures/  — all plot panels (PDF + PNG)
    data/     — processed AnnData objects
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
import anndata as ad

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
sc.settings.n_jobs = 4

# ── Paths ────────────────────────────────────────────────────────────────────
OUT     = "bm_analysis_out"
FIG_DIR = os.path.join(OUT, "figures")
DAT_DIR = os.path.join(OUT, "data")
RAW_DIR = os.path.join(OUT, "raw_geo")
for d in [FIG_DIR, DAT_DIR, RAW_DIR]:
    os.makedirs(d, exist_ok=True)

sc.settings.figdir = FIG_DIR


# ═══════════════════════════════════════════════════════════════════════════
# 1. DOWNLOAD DATA FROM GEO
# ═══════════════════════════════════════════════════════════════════════════

def download_geo(geo_id="GSE266330"):
    """
    Download count matrices from GEO.
    GSE266330 provides per-sample matrices (barcodes/features/matrix).
    Falls back to h5 if available.
    """
    import urllib.request, gzip, shutil, tarfile

    # Try to fetch the series matrix for metadata
    meta_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{geo_id[:6]}nnn/{geo_id}/matrix/"
        f"{geo_id}_series_matrix.txt.gz"
    )
    meta_path = os.path.join(RAW_DIR, f"{geo_id}_series_matrix.txt.gz")
    if not os.path.exists(meta_path):
        print(f"  Downloading series matrix from GEO …")
        try:
            urllib.request.urlretrieve(meta_url, meta_path)
        except Exception as e:
            print(f"  WARNING: Could not download series matrix: {e}")

    # Download supplementary tar (contains per-sample 10x dirs or h5 files)
    supp_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{geo_id[:6]}nnn/{geo_id}/suppl/"
        f"{geo_id}_RAW.tar"
    )
    tar_path = os.path.join(RAW_DIR, f"{geo_id}_RAW.tar")
    if not os.path.exists(tar_path):
        print(f"  Downloading {geo_id}_RAW.tar (~may be large) …")
        urllib.request.urlretrieve(supp_url, tar_path)
        print("  Download complete.")

    extract_dir = os.path.join(RAW_DIR, "extracted")
    if not os.path.exists(extract_dir):
        print("  Extracting RAW tar …")
        os.makedirs(extract_dir)
        with tarfile.open(tar_path) as tar:
            tar.extractall(extract_dir)

    return extract_dir, meta_path


def parse_geo_metadata(meta_gz_path):
    """
    Parse GEO series matrix for sample-level clinical metadata.
    Returns a DataFrame indexed by sample/GSM ID.
    """
    if not os.path.exists(meta_gz_path):
        return pd.DataFrame()

    import gzip, re
    rows = {}
    with gzip.open(meta_gz_path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Sample_"):
                key, _, val = line.partition("\t")
                key = key.lstrip("!")
                vals = [v.strip('"') for v in val.split("\t")]
                rows[key] = vals
    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).T
    # transpose so samples are rows
    if "Sample_geo_accession" in df.index:
        df.columns = df.loc["Sample_geo_accession"].tolist()
        df = df.drop(index="Sample_geo_accession").T
    return df


def load_10x_samples(extract_dir, meta_df=None):
    """
    Read all 10x sample directories (or .h5 files) found under extract_dir.
    Returns a single concatenated AnnData with obs column 'sample_id'.
    """
    import gzip, re

    adatas = []

    # look for h5 files first
    h5_files = []
    for root, dirs, files in os.walk(extract_dir):
        for f in files:
            if f.endswith(".h5") or f.endswith(".h5ad"):
                h5_files.append(os.path.join(root, f))

    if h5_files:
        print(f"  Found {len(h5_files)} h5 files.")
        for path in sorted(h5_files):
            sample_id = os.path.basename(path).split(".")[0]
            try:
                adata = sc.read_10x_h5(path) if path.endswith(".h5") else sc.read_h5ad(path)
                adata.var_names_make_unique()
                adata.obs["sample_id"] = sample_id
                adatas.append(adata)
            except Exception as e:
                print(f"    Skipping {path}: {e}")

    else:
        # look for subdirectories with barcodes/features/matrix
        sample_dirs = []
        for root, dirs, files in os.walk(extract_dir):
            gz_files = [f for f in files if f.endswith(".gz")]
            has_barcodes = any("barcodes" in f for f in gz_files)
            has_matrix   = any("matrix" in f for f in gz_files)
            if has_barcodes and has_matrix:
                sample_dirs.append(root)

        print(f"  Found {len(sample_dirs)} 10x sample directories.")
        for d in sorted(sample_dirs):
            sample_id = os.path.basename(d)
            try:
                adata = sc.read_10x_mtx(d, var_names="gene_symbols", cache=True)
                adata.var_names_make_unique()
                adata.obs["sample_id"] = sample_id
                adatas.append(adata)
            except Exception as e:
                print(f"    Skipping {d}: {e}")

    if not adatas:
        raise RuntimeError(
            "No count data found in extracted GEO directory.\n"
            f"  Please check: {extract_dir}\n"
            "  Expected: .h5 files OR subdirs with barcodes/features/matrix.mtx.gz"
        )

    print(f"  Concatenating {len(adatas)} samples …")
    adata = ad.concat(adatas, label="sample_id", keys=[a.obs["sample_id"][0] for a in adatas],
                      merge="same")
    adata.obs_names_make_unique()

    # Attach clinical metadata if available
    if meta_df is not None and not meta_df.empty:
        shared_cols = meta_df.columns.tolist()
        meta_sub = meta_df.reindex(adata.obs["sample_id"])
        meta_sub.index = adata.obs.index
        for col in shared_cols:
            adata.obs[col] = meta_sub[col].values

    return adata


# ═══════════════════════════════════════════════════════════════════════════
# 2. PRE-PROCESSING (full dataset)
# ═══════════════════════════════════════════════════════════════════════════

def preprocess(adata, min_genes=200, min_cells=3, max_pct_mt=25):
    print("  QC filtering …")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()

    print(f"  After QC: {adata.n_obs} cells, {adata.n_vars} genes")

    # Normalise & log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep raw counts for DEG

    # HVG
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample_id")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
    sc.tl.umap(adata)

    # Broad clustering for cell-type annotation
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_broad")
    return adata


# ═══════════════════════════════════════════════════════════════════════════
# 3. MACROPHAGE ISOLATION
# ═══════════════════════════════════════════════════════════════════════════

# Canonical macrophage markers (human)
MACRO_MARKERS = ["CD68", "CSF1R", "MRC1", "CD163", "ADGRE1", "ITGAM",
                 "FCGR3A", "MSR1", "C1QA", "C1QB", "AIF1"]

def annotate_and_extract_macrophages(adata):
    """
    Score each cell for macrophage signature, then extract high-scoring cells.
    Also keeps cells from clusters where macrophage score is highest.
    """
    # Only use markers present in the data
    markers_present = [g for g in MACRO_MARKERS if g in adata.var_names]
    print(f"  Macrophage markers found: {markers_present}")

    sc.tl.score_genes(adata, gene_list=markers_present, score_name="macro_score")

    # Per cluster mean score → annotate macrophage-enriched clusters
    cluster_scores = adata.obs.groupby("leiden_broad")["macro_score"].mean()
    macro_clusters = cluster_scores[cluster_scores > cluster_scores.quantile(0.6)].index.tolist()
    print(f"  Macrophage-enriched broad clusters: {macro_clusters}")

    mask = (
        adata.obs["leiden_broad"].isin(macro_clusters) |
        (adata.obs["macro_score"] > adata.obs["macro_score"].quantile(0.85))
    )
    mac = adata[mask].copy()
    print(f"  Macrophages extracted: {mac.n_obs} cells")
    return mac


# ═══════════════════════════════════════════════════════════════════════════
# 4. MACROPHAGE SUB-CLUSTERING
# ═══════════════════════════════════════════════════════════════════════════

def subcluster_macrophages(mac, resolution=0.4):
    print("  Re-processing macrophage subset …")
    # Re-run HVG/PCA/UMAP on the subset
    sc.pp.highly_variable_genes(mac, n_top_genes=2000, batch_key="sample_id")
    sc.pp.scale(mac, max_value=10)
    sc.tl.pca(mac, svd_solver="arpack")
    sc.pp.neighbors(mac, n_neighbors=20, n_pcs=20)
    sc.tl.umap(mac)
    sc.tl.leiden(mac, resolution=resolution, key_added="mac_cluster")
    n_clusters = mac.obs["mac_cluster"].nunique()
    print(f"  Macrophage sub-clusters: {n_clusters}")
    return mac


# ═══════════════════════════════════════════════════════════════════════════
# 5. METADATA HARMONISATION
# ═══════════════════════════════════════════════════════════════════════════

# Map common GEO characteristic fields → standardised column names.
# Adjust these keys to match the actual GSE266330 series matrix fields.
META_MAP = {
    "cancer_type"        : ["cancer type", "tumor type", "tissue type", "cancer_type",
                             "Sample_characteristics_ch1_cancer type"],
    "condition"          : ["condition", "disease state", "metastasis", "tumor/normal",
                             "Sample_characteristics_ch1_disease state"],
    "treatment"          : ["treatment", "prior treatment", "Sample_characteristics_ch1_treatment"],
    "treatment_response" : ["response", "treatment response", "Sample_characteristics_ch1_response"],
    "gender"             : ["gender", "sex", "Sample_characteristics_ch1_Sex",
                             "Sample_characteristics_ch1_gender"],
    "tissue_origin"      : ["tissue", "primary site", "tissue origin",
                             "Sample_characteristics_ch1_tissue"],
}

def harmonise_metadata(obs):
    """
    Try to find standardised columns in obs; add them if not already present.
    """
    obs = obs.copy()
    for std_col, candidates in META_MAP.items():
        if std_col not in obs.columns:
            for cand in candidates:
                # case-insensitive partial match
                matches = [c for c in obs.columns
                           if cand.lower() in c.lower()]
                if matches:
                    obs[std_col] = obs[matches[0]]
                    break
            else:
                obs[std_col] = "unknown"
    return obs


# ═══════════════════════════════════════════════════════════════════════════
# 6. FIGURES
# ═══════════════════════════════════════════════════════════════════════════

PALETTE = sc.pl.palettes.default_20

def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"),
                    bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {name}")


# ── Fig 1: Sub-cluster UMAP + cancer type ────────────────────────────────────
def fig_umap_overview(mac):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    sc.pl.umap(mac, color="mac_cluster", legend_loc="on data",
               title="Macrophage sub-clusters", ax=axes[0, 0], show=False,
               palette=PALETTE)

    ct_col = "cancer_type" if "cancer_type" in mac.obs else "sample_id"
    sc.pl.umap(mac, color=ct_col, title="Cancer type of origin",
               ax=axes[0, 1], show=False)

    # Condition (metastasis vs control)
    cond_col = "condition" if "condition" in mac.obs else "mac_cluster"
    sc.pl.umap(mac, color=cond_col, title="Condition",
               ax=axes[1, 0], show=False)

    sc.pl.umap(mac, color="sample_id", title="Sample",
               ax=axes[1, 1], show=False)

    fig.suptitle("Macrophage Sub-clustering — All Samples (GSE266330)",
                 fontsize=14, y=1.02)
    save(fig, "01_umap_overview")


# ── Fig 2: Stacked bar — cancer type per cluster ─────────────────────────────
def fig_cancer_type_bar(mac):
    ct_col = "cancer_type" if "cancer_type" in mac.obs.columns else None
    if ct_col is None:
        print("  Skipping cancer-type bar: no cancer_type column")
        return

    df = (mac.obs.groupby(["mac_cluster", ct_col])
               .size()
               .reset_index(name="n"))
    totals = df.groupby("mac_cluster")["n"].transform("sum")
    df["frac"] = df["n"] / totals

    pivot = df.pivot(index="mac_cluster", columns=ct_col, values="frac").fillna(0)

    fig, ax = plt.subplots(figsize=(max(6, len(pivot)*0.8), 5))
    pivot.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", edgecolor="none")
    ax.set_xlabel("Macrophage sub-cluster"); ax.set_ylabel("Fraction")
    ax.set_title("Cancer type composition per macrophage cluster")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    plt.xticks(rotation=0)
    save(fig, "02_cancer_type_bar")


# ── Fig 3: ER (ESR1) + SP1 feature plots + dot plot ─────────────────────────
TARGET_GENES = ["ESR1", "SP1"]

def fig_er_sp1(mac):
    # Use raw counts for expression
    adata_plot = mac.raw.to_adata() if mac.raw is not None else mac

    genes_present = [g for g in TARGET_GENES if g in adata_plot.var_names]
    if not genes_present:
        print(f"  WARNING: genes {TARGET_GENES} not found in data.")
        return

    n = len(genes_present)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes_present):
        sc.pl.umap(adata_plot, color=gene, ax=ax, show=False,
                   title=f"{gene} expression", color_map="Reds", vmin=0)

    save(fig, "03a_feature_plots_ESR1_SP1")

    # Dot plot across clusters
    fig2, ax2 = plt.subplots(figsize=(max(6, mac.obs["mac_cluster"].nunique() * 0.7), 4))
    sc.pl.dotplot(mac, var_names=genes_present, groupby="mac_cluster",
                  ax=ax2, show=False, title="ESR1 & SP1 across macrophage sub-clusters",
                  use_raw=True)
    save(fig2, "03b_dotplot_ESR1_SP1")


# ── Fig 4: Metastasis vs. control enrichment ─────────────────────────────────
def fig_metastasis_enrichment(mac):
    cond_col = "condition" if "condition" in mac.obs.columns else None
    if cond_col is None:
        print("  Skipping metastasis enrichment: no condition column")
        return

    # Proportion of each cluster per condition
    df = (mac.obs.groupby([cond_col, "mac_cluster"])
               .size()
               .reset_index(name="n"))
    totals = df.groupby(cond_col)["n"].transform("sum")
    df["frac"] = df["n"] / totals

    pivot = df.pivot(index=cond_col, columns="mac_cluster", values="frac").fillna(0)

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns)*0.8), 4))
    pivot.T.plot(kind="bar", ax=ax, edgecolor="none")
    ax.set_xlabel("Macrophage sub-cluster"); ax.set_ylabel("Fraction of condition")
    ax.set_title("Macrophage cluster proportion — Metastasis vs. Control")
    ax.legend(title=cond_col, bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.xticks(rotation=0)
    save(fig, "04_metastasis_vs_control")

    # DEG heatmap per cluster (met vs ctrl)
    # Only run if enough cells per group
    try:
        sc.tl.rank_genes_groups(mac, groupby="mac_cluster", use_raw=True,
                                method="wilcoxon", key_added="rank_mac")
        top_n = 5
        top_genes = []
        for cl in mac.obs["mac_cluster"].cat.categories:
            try:
                result = sc.get.rank_genes_groups_df(mac, group=str(cl), key="rank_mac")
                top_genes += result.head(top_n)["names"].tolist()
            except Exception:
                pass
        top_genes = list(dict.fromkeys(top_genes))  # deduplicate, preserve order

        if top_genes:
            fig3, ax3 = plt.subplots(figsize=(12, max(4, len(top_genes) * 0.25)))
            sc.pl.heatmap(mac, var_names=top_genes, groupby="mac_cluster",
                          use_raw=True, ax=ax3, show=False,
                          cmap="RdBu_r", vcenter=0,
                          figsize=(12, max(4, len(top_genes) * 0.25)))
            ax3.set_title("Top DEGs per macrophage cluster")
            save(fig3, "04b_deg_heatmap")
    except Exception as e:
        print(f"  DEG heatmap skipped: {e}")


# ── Fig 5: Treatment response ────────────────────────────────────────────────
def fig_treatment_response(mac):
    resp_col = "treatment_response" if "treatment_response" in mac.obs.columns else None
    if resp_col is None or mac.obs[resp_col].eq("unknown").all():
        print("  Skipping treatment response: metadata not available in this dataset")
        # Make an informative placeholder note
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5,
                "Treatment response metadata\nnot available in GSE266330.\n"
                "Manual annotation or external\nclinical data linkage required.",
                ha="center", va="center", fontsize=13,
                bbox=dict(boxstyle="round", facecolor="#f0f4f8", edgecolor="#aabbcc"))
        ax.axis("off")
        ax.set_title("Fig 5 — Treatment Response", fontsize=12)
        save(fig, "05_treatment_response_placeholder")
        return

    df = (mac.obs.groupby([resp_col, "mac_cluster"])
               .size().reset_index(name="n"))
    totals = df.groupby(resp_col)["n"].transform("sum")
    df["frac"] = df["n"] / totals
    pivot = df.pivot(index=resp_col, columns="mac_cluster", values="frac").fillna(0)

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns)), 4))
    pivot.T.plot(kind="bar", ax=ax, edgecolor="none")
    ax.set_xlabel("Macrophage sub-cluster"); ax.set_ylabel("Fraction")
    ax.set_title("Macrophage cluster ~ Treatment Response")
    ax.legend(title=resp_col, bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.xticks(rotation=0)
    save(fig, "05_treatment_response")


# ── Fig 6: Gender + tissue origin ────────────────────────────────────────────
def fig_clinical_metadata(mac):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, col, title in zip(
        axes,
        ["gender", "tissue_origin"],
        ["Gender per macrophage cluster", "Tissue origin per macrophage cluster"]
    ):
        if col not in mac.obs.columns or mac.obs[col].eq("unknown").all():
            ax.text(0.5, 0.5, f"'{col}' not in metadata",
                    ha="center", va="center", fontsize=12)
            ax.axis("off"); ax.set_title(title)
            continue

        df = (mac.obs.groupby([col, "mac_cluster"])
                   .size().reset_index(name="n"))
        totals = df.groupby(col)["n"].transform("sum")
        df["frac"] = df["n"] / totals
        pivot = df.pivot(index=col, columns="mac_cluster", values="frac").fillna(0)
        pivot.T.plot(kind="bar", ax=ax, edgecolor="none", colormap="Set2")
        ax.set_xlabel("Macrophage sub-cluster"); ax.set_ylabel("Fraction")
        ax.set_title(title)
        ax.legend(title=col, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
        ax.tick_params(axis="x", rotation=0)

    fig.suptitle("Clinical Metadata per Macrophage Sub-cluster", fontsize=13)
    save(fig, "06_gender_tissue_origin")


# ── Fig 7: Summary — co-expression score ER+SP1-high cluster ─────────────────
def fig_er_sp1_coexpr(mac):
    adata_plot = mac.raw.to_adata() if mac.raw is not None else mac
    genes = [g for g in ["ESR1", "SP1"] if g in adata_plot.var_names]
    if len(genes) < 2:
        print("  Skipping co-expression score: need both ESR1 and SP1")
        return

    sc.tl.score_genes(adata_plot, gene_list=genes, score_name="ER_SP1_score")
    mac.obs["ER_SP1_score"] = adata_plot.obs["ER_SP1_score"].values

    # Violin per cluster
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.violin(mac, keys="ER_SP1_score", groupby="mac_cluster",
                 ax=axes[0], show=False, rotation=0)
    axes[0].set_title("ER + SP1 co-expression score per cluster")

    sc.pl.umap(mac, color="ER_SP1_score", ax=axes[1], show=False,
               color_map="RdYlBu_r", title="ER + SP1 co-expression (UMAP)")

    # Annotate the top cluster
    top_cl = mac.obs.groupby("mac_cluster")["ER_SP1_score"].mean().idxmax()
    mac.obs["ER_SP1_top"] = (mac.obs["mac_cluster"] == top_cl).map(
        {True: f"Cluster {top_cl} (ER+/SP1-high)", False: "Other"})
    print(f"  ER+/SP1-high cluster: {top_cl}")

    save(fig, "07_ER_SP1_coexpression")
    return top_cl


# ═══════════════════════════════════════════════════════════════════════════
# 7. MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    processed_path = os.path.join(DAT_DIR, "mac_subclustered.h5ad")
    full_path       = os.path.join(DAT_DIR, "full_preprocessed.h5ad")

    if os.path.exists(processed_path):
        print("Loading cached macrophage AnnData …")
        mac = sc.read_h5ad(processed_path)
    else:
        # ── Download ──────────────────────────────────────────────────────
        print("Step 1: Downloading GSE266330 from GEO …")
        extract_dir, meta_gz = download_geo("GSE266330")

        # ── Parse metadata ────────────────────────────────────────────────
        print("Step 2: Parsing clinical metadata …")
        meta_df = parse_geo_metadata(meta_gz)
        print(f"  Metadata shape: {meta_df.shape}")

        # ── Load count data ───────────────────────────────────────────────
        print("Step 3: Loading count matrices …")
        adata = load_10x_samples(extract_dir, meta_df)
        print(f"  Full dataset: {adata.n_obs} cells, {adata.n_vars} genes")

        # ── Pre-process ───────────────────────────────────────────────────
        print("Step 4: Pre-processing …")
        adata = preprocess(adata)
        adata.write_h5ad(full_path)

        # ── Macrophage isolation ──────────────────────────────────────────
        print("Step 5: Isolating macrophages …")
        mac = annotate_and_extract_macrophages(adata)

        # ── Sub-clustering ────────────────────────────────────────────────
        print("Step 6: Sub-clustering macrophages …")
        mac = subcluster_macrophages(mac)

        # ── Harmonise metadata columns ────────────────────────────────────
        print("Step 7: Harmonising metadata …")
        mac.obs = harmonise_metadata(mac.obs)
        mac.write_h5ad(processed_path)
        print(f"  Saved: {processed_path}")

    # ── Figures ───────────────────────────────────────────────────────────
    print("\nGenerating figures …")
    fig_umap_overview(mac)
    fig_cancer_type_bar(mac)
    fig_er_sp1(mac)
    fig_metastasis_enrichment(mac)
    fig_treatment_response(mac)
    fig_clinical_metadata(mac)
    top_cluster = fig_er_sp1_coexpr(mac)

    print(f"\nAll done. Figures in: {FIG_DIR}")
    print(f"  ER+/SP1-high cluster: {top_cluster}")
    print(f"  Macrophage AnnData: {processed_path}")


if __name__ == "__main__":
    main()
