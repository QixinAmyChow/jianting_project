"""
Load macrophage MTX from R extraction, sub-cluster, generate all figures.
"""
import os, warnings
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
sc.settings.n_jobs = 4

MAC_DIR = "bm_analysis_out/data/macrophages"
FIG_DIR = "bm_analysis_out/figures"
DAT_DIR = "bm_analysis_out/data"
os.makedirs(FIG_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

PALETTE = sc.pl.palettes.default_20

def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"),
                    bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {name}")


# ── 1. Load data ──────────────────────────────────────────────────────────────
cached = os.path.join(DAT_DIR, "mac_subclustered.h5ad")

if os.path.exists(cached):
    print("Loading cached AnnData …")
    mac = sc.read_h5ad(cached)
else:
    print("Loading MTX …")
    mat = scipy.io.mmread(os.path.join(MAC_DIR, "mac_counts.mtx")).T  # cells x genes
    genes    = open(os.path.join(MAC_DIR, "genes.txt")).read().splitlines()
    barcodes = open(os.path.join(MAC_DIR, "barcodes.txt")).read().splitlines()
    meta     = pd.read_csv(os.path.join(MAC_DIR, "mac_metadata.csv"), index_col=0)

    mac = ad.AnnData(X=scipy.sparse.csr_matrix(mat),
                     obs=meta.loc[barcodes] if set(barcodes).issubset(meta.index)
                         else meta.reindex(barcodes),
                     var=pd.DataFrame(index=genes))
    mac.obs_names = barcodes
    mac.var_names = genes
    mac.var_names_make_unique()
    print(f"  Loaded: {mac.n_obs} cells x {mac.n_vars} genes")

    # ── 2. HVG + PCA + UMAP ──────────────────────────────────────────────────
    print("HVG / PCA / UMAP …")
    sc.pp.highly_variable_genes(mac, n_top_genes=2000, flavor="seurat")
    sc.tl.pca(mac, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(mac, n_neighbors=20, n_pcs=20)
    sc.tl.umap(mac)
    sc.tl.leiden(mac, resolution=0.4, key_added="mac_cluster")
    n_cl = mac.obs["mac_cluster"].nunique()
    print(f"  Macrophage sub-clusters: {n_cl}")

    # Harmonise metadata columns
    META_MAP = {
        "cancer_type" : ["cancer", "cancer type", "tumor type"],
        "condition"   : ["condition", "disease state", "tissue.origin"],
        "gender"      : ["gender", "sex"],
        "tissue_origin": ["tissue.origin", "tissue origin", "tissue"],
    }
    for std, cands in META_MAP.items():
        if std not in mac.obs.columns:
            for c in cands:
                hits = [col for col in mac.obs.columns if c.lower() in col.lower()]
                if hits:
                    mac.obs[std] = mac.obs[hits[0]]
                    break
            else:
                mac.obs[std] = "unknown"

    # Coerce all object columns to string to avoid h5py write errors
    for col in mac.obs.columns:
        if mac.obs[col].dtype == object:
            mac.obs[col] = mac.obs[col].astype(str)
    mac.write_h5ad(cached)
    print(f"  Saved: {cached}")

print(f"\nClusters: {sorted(mac.obs['mac_cluster'].unique())}")

# ── 3. Figures ────────────────────────────────────────────────────────────────

# Fig 1: UMAP overview
print("Fig 1: UMAP overview …")
fig, axes = plt.subplots(2, 2, figsize=(14, 14))
sc.pl.umap(mac, color="mac_cluster", legend_loc="on data",
           title="Macrophage sub-clusters", ax=axes[0, 0], show=False, palette=PALETTE)
ct_col = "cancer_type" if "cancer_type" in mac.obs else "sample_id"
sc.pl.umap(mac, color=ct_col, title="Cancer type", ax=axes[0, 1], show=False)
sc.pl.umap(mac, color="sample_id", title="Sample", ax=axes[1, 0], show=False,
           legend_fontsize=6)
leg = axes[1, 0].get_legend()
if leg:
    leg.set_bbox_to_anchor((0.5, -0.02))
    leg.set_loc("upper center")
    leg._ncols = 4
    leg.get_frame().set_visible(False)
cond_col = "condition" if "condition" in mac.obs else "mac_cluster"
sc.pl.umap(mac, color=cond_col, title="Condition", ax=axes[1, 1], show=False)
fig.subplots_adjust(hspace=0.45, wspace=0.35)
fig.suptitle("Macrophage Sub-clustering — GSE266330 (31 samples)", fontsize=13, y=1.02)
save(fig, "01_umap_overview")

# Fig 2: Stacked bar — cancer type per cluster
print("Fig 2: Cancer type bar …")
ct_col = "cancer_type" if "cancer_type" in mac.obs.columns else None
if ct_col:
    df = mac.obs.groupby(["mac_cluster", ct_col]).size().reset_index(name="n")
    df["frac"] = df["n"] / df.groupby("mac_cluster")["n"].transform("sum")
    pivot = df.pivot(index="mac_cluster", columns=ct_col, values="frac").fillna(0)
    fig, ax = plt.subplots(figsize=(max(6, len(pivot)*0.9), 5))
    pivot.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", edgecolor="none")
    ax.set_xlabel("Macrophage sub-cluster"); ax.set_ylabel("Fraction")
    ax.set_title("Cancer type per macrophage cluster")
    ax.legend(bbox_to_anchor=(1.01,1), loc="upper left", fontsize=9)
    plt.xticks(rotation=0)
    save(fig, "02_cancer_type_bar")

# Fig 3: ESR1 + SP1 feature plots + dot plot
print("Fig 3: ESR1 / SP1 …")
TARGET = [g for g in ["ESR1", "SP1"] if g in mac.var_names]
if TARGET:
    fig, axes = plt.subplots(1, len(TARGET), figsize=(6*len(TARGET), 5))
    if len(TARGET) == 1: axes = [axes]
    for ax, gene in zip(axes, TARGET):
        sc.pl.umap(mac, color=gene, ax=ax, show=False,
                   title=f"{gene} expression", color_map="Reds", vmin=0)
    save(fig, "03a_feature_plots_ESR1_SP1")

    fig2, ax2 = plt.subplots(figsize=(max(6, mac.obs["mac_cluster"].nunique()*0.8), 4))
    sc.pl.dotplot(mac, var_names=TARGET, groupby="mac_cluster",
                  ax=ax2, show=False, title="ESR1 & SP1 across clusters")
    save(fig2, "03b_dotplot_ESR1_SP1")
else:
    print("  WARNING: ESR1/SP1 not found in var_names")

# Fig 4: Gender + tissue origin
print("Fig 4: Clinical metadata …")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, col, title in zip(axes,
    ["gender", "tissue_origin"],
    ["Gender per cluster", "Tissue origin per cluster"]):
    if col not in mac.obs or mac.obs[col].eq("unknown").all():
        ax.text(0.5, 0.5, f"'{col}' not available", ha="center", va="center")
        ax.axis("off"); ax.set_title(title); continue
    df = mac.obs.groupby([col, "mac_cluster"]).size().reset_index(name="n")
    df["frac"] = df["n"] / df.groupby(col)["n"].transform("sum")
    pivot = df.pivot(index=col, columns="mac_cluster", values="frac").fillna(0)
    pivot.T.plot(kind="bar", ax=ax, edgecolor="none", colormap="Set2")
    ax.set_xlabel("Cluster"); ax.set_ylabel("Fraction")
    ax.set_title(title)
    ax.legend(title=col, bbox_to_anchor=(1.01,1), loc="upper left", fontsize=9)
    ax.tick_params(axis="x", rotation=0)
fig.suptitle("Clinical Metadata per Macrophage Sub-cluster", fontsize=13)
save(fig, "04_gender_tissue_origin")

# Fig 5: ER+SP1 co-expression score
print("Fig 5: ER+SP1 co-expression …")
genes = [g for g in ["ESR1","SP1"] if g in mac.var_names]
if len(genes) == 2:
    sc.tl.score_genes(mac, gene_list=genes, score_name="ER_SP1_score")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.violin(mac, keys="ER_SP1_score", groupby="mac_cluster",
                 ax=axes[0], show=False, rotation=0)
    axes[0].set_title("ER + SP1 co-expression score per cluster")
    sc.pl.umap(mac, color="ER_SP1_score", ax=axes[1], show=False,
               color_map="RdYlBu_r", title="ER + SP1 co-expression (UMAP)")
    top_cl = mac.obs.groupby("mac_cluster")["ER_SP1_score"].mean().idxmax()
    print(f"  ER+/SP1-high cluster: {top_cl}")
    save(fig, "05_ER_SP1_coexpression")

    # DEG heatmap
    sc.tl.rank_genes_groups(mac, groupby="mac_cluster", method="wilcoxon",
                            key_added="rank_mac")
    top_genes = []
    for cl in mac.obs["mac_cluster"].cat.categories:
        try:
            res = sc.get.rank_genes_groups_df(mac, group=str(cl), key="rank_mac")
            top_genes += res.head(3)["names"].tolist()
        except Exception:
            pass
    top_genes = list(dict.fromkeys(top_genes))
    if top_genes:
        n_genes = len(top_genes)
        n_clusters = mac.obs["mac_cluster"].nunique()
        fig_h = max(6, n_genes * 0.25)
        fig_w = max(10, n_clusters * 0.8)
        sc.pl.heatmap(mac, var_names=top_genes, groupby="mac_cluster",
                      cmap="RdBu_r", vcenter=0, show=False,
                      figsize=(fig_w, fig_h),
                      save="_deg_heatmap.png")
        os.rename(os.path.join(FIG_DIR, "heatmap_deg_heatmap.png"),
                  os.path.join(FIG_DIR, "06_deg_heatmap.png"))
        print("  Saved: 06_deg_heatmap")
else:
    top_cl = None
    print("  Skipping ER+SP1 score (genes not found)")

print(f"\nAll figures saved to: {FIG_DIR}")
print(f"ER+/SP1-high cluster: {top_cl}")
print("DONE")
