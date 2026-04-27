"""
BC ER+/ER- Macrophage Sub-clustering Analysis
Source : GSE266330 — 47.integrated.h5ad (from integrated_Seurat_objects.zip)
Goal   : Same 3-figure analysis as BM project update slides, but BC samples
         are split into BC_ER+ and BC_ER- based on macrophage ESR1 expression.

Figures (saved to ./bm_analysis_out/figures/):
    Fig 1 — Macrophage sub-cluster UMAP
    Fig 2 — Cancer type composition per cluster  (BC → BC_ER+ / BC_ER-)
    Fig 3 — ESR1 & SP1 feature plots + dot plot
"""

import os, zipfile, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 2
sc.settings.n_jobs   = 4

ZIP_PATH  = "bm_analysis_out/raw_geo/integrated_Seurat_objects.zip"
H5AD_NAME = "integrated_Seurat_objects/47.integrated.h5ad"
FIG_DIR   = "bm_analysis_out/figures"
DAT_DIR   = "bm_analysis_out/data"
CACHED    = os.path.join(DAT_DIR, "mac_bc_er_subclustered.h5ad")

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DAT_DIR, exist_ok=True)
sc.settings.figdir = FIG_DIR

PALETTE = sc.pl.palettes.default_20


def save(fig, name):
    for ext in ["pdf", "png"]:
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"),
                    bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {name}")


# ── Colour map for cancer types (BC split gets two shades of pink/red) ─────────
CANCER_COLORS = {
    "BC_ER+"    : "#d62728",   # brick red
    "BC_ER-"    : "#ff9896",   # soft pink
    "BC"        : "#e377c2",   # fallback if split not applied
    "BDC"       : "#ff7f0e",
    "CC"        : "#2ca02c",
    "EC"        : "#9467bd",
    "KC"        : "#8c564b",
    "LC"        : "#17becf",
    "PC"        : "#bcbd22",
    "TC"        : "#7f7f7f",
    "ctrl"      : "#aec7e8",
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1.  LOAD / EXTRACT h5ad
# ═══════════════════════════════════════════════════════════════════════════════

def extract_h5ad():
    """Extract 47.integrated.h5ad from the ZIP if not already done."""
    out_path = "/tmp/47.integrated.h5ad"
    if os.path.exists(out_path):
        print(f"  h5ad already extracted: {out_path}")
        return out_path

    print(f"  Extracting {H5AD_NAME} from ZIP …")
    with zipfile.ZipFile(ZIP_PATH) as zf:
        with zf.open(H5AD_NAME) as src, open(out_path, "wb") as dst:
            chunk = 1 << 23   # 8 MB chunks
            while True:
                buf = src.read(chunk)
                if not buf:
                    break
                dst.write(buf)
    print(f"  Extracted: {out_path}")
    return out_path


# ═══════════════════════════════════════════════════════════════════════════════
# 2.  ISOLATE MACROPHAGES
# ═══════════════════════════════════════════════════════════════════════════════

MACRO_MARKERS = ["CD68", "CSF1R", "MRC1", "CD163", "FCGR3A",
                 "MSR1", "C1QA", "C1QB", "AIF1", "ITGAM"]

CELLTYPE_CANDIDATES = [
    "celltype", "cell_type", "CellType", "Celltype",
    "major_celltype", "majorCelltype", "annotation",
    "leiden", "seurat_clusters",
]
MACRO_KEYWORDS = ["macro", "Macro", "MACRO", "myeloid", "Myeloid",
                  "MΦ", "Mφ", "monocyte", "Monocyte"]


def find_celltype_col(obs):
    for c in CELLTYPE_CANDIDATES:
        if c in obs.columns:
            return c
    return None


def isolate_macrophages(adata):
    ct_col = find_celltype_col(adata.obs)

    if ct_col:
        vals = adata.obs[ct_col].astype(str)
        mask = vals.str.contains("|".join(MACRO_KEYWORDS), case=False, na=False)
        if mask.sum() > 100:
            print(f"  Using cell-type annotation '{ct_col}': {mask.sum()} macrophages")
            return adata[mask].copy()

    # Fallback: score with canonical markers
    present = [g for g in MACRO_MARKERS if g in adata.var_names]
    print(f"  No clean annotation found; scoring with markers: {present}")
    sc.tl.score_genes(adata, gene_list=present, score_name="_macro_score")
    thresh = adata.obs["_macro_score"].quantile(0.80)
    mask   = adata.obs["_macro_score"] > thresh
    print(f"  Score threshold ({thresh:.3f}): {mask.sum()} cells")
    return adata[mask].copy()


# ═══════════════════════════════════════════════════════════════════════════════
# 3.  SUB-CLUSTER MACROPHAGES
# ═══════════════════════════════════════════════════════════════════════════════

def subcluster(mac, resolution=0.4):
    print("  Re-running HVG / PCA / UMAP on macrophage subset …")
    batch_key = "sample_id" if "sample_id" in mac.obs.columns else None
    sc.pp.highly_variable_genes(mac, n_top_genes=2000,
                                flavor="seurat", batch_key=batch_key)
    sc.tl.pca(mac, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(mac, n_neighbors=20, n_pcs=20)
    sc.tl.umap(mac)
    sc.tl.leiden(mac, resolution=resolution, key_added="mac_cluster")
    print(f"  Sub-clusters: {mac.obs['mac_cluster'].nunique()}")
    return mac


# ═══════════════════════════════════════════════════════════════════════════════
# 4.  HARMONISE METADATA  +  BC ER+/ER- SPLIT
# ═══════════════════════════════════════════════════════════════════════════════

META_MAP = {
    "cancer_type"  : ["cancer_type", "cancer type", "tumor_type", "tumor type",
                      "cancertype", "CancerType", "Sample_cancer_type"],
    "sample_id"    : ["sample_id", "orig.ident", "SampleID", "sample"],
    "condition"    : ["condition", "disease state", "tissue.origin"],
    "gender"       : ["gender", "sex", "Sex"],
    "tissue_origin": ["tissue.origin", "tissue origin", "tissue"],
}


def harmonise_metadata(obs):
    obs = obs.copy()
    for std, cands in META_MAP.items():
        if std not in obs.columns:
            for c in cands:
                hits = [col for col in obs.columns if c.lower() in col.lower()]
                if hits:
                    obs[std] = obs[hits[0]]
                    break
            else:
                obs[std] = "unknown"
    return obs


def add_bc_er_split(mac):
    """
    Within BC macrophages, split by ESR1 expression into ER+ and ER-.
    Threshold: median ESR1 expression among BC cells.
    Adds 'cancer_type_er' column (same as cancer_type for non-BC; BC → BC_ER+/BC_ER-).
    """
    if "cancer_type" not in mac.obs.columns:
        mac.obs["cancer_type_er"] = "unknown"
        return mac

    # Grab ESR1 expression
    if "ESR1" not in mac.var_names:
        print("  WARNING: ESR1 not found — BC will not be split")
        mac.obs["cancer_type_er"] = mac.obs["cancer_type"].astype(str)
        return mac

    esr1 = np.asarray(mac[:, "ESR1"].X.todense()).flatten() \
           if hasattr(mac[:, "ESR1"].X, "todense") \
           else mac[:, "ESR1"].X.flatten()

    bc_mask  = mac.obs["cancer_type"].astype(str).str.upper().str.contains("BC")
    bc_esr1  = esr1[bc_mask]
    threshold = np.median(bc_esr1[bc_esr1 > 0]) if (bc_esr1 > 0).any() else np.median(bc_esr1)
    print(f"  BC cells: {bc_mask.sum()}  |  ESR1 threshold (median positive): {threshold:.4f}")

    ct_er = mac.obs["cancer_type"].astype(str).copy()
    ct_er[bc_mask & (esr1 >= threshold)] = "BC_ER+"
    ct_er[bc_mask & (esr1 <  threshold)] = "BC_ER-"
    mac.obs["cancer_type_er"] = ct_er.values

    n_pos = (ct_er == "BC_ER+").sum()
    n_neg = (ct_er == "BC_ER-").sum()
    print(f"  BC_ER+: {n_pos}  |  BC_ER-: {n_neg}")
    return mac


# ═══════════════════════════════════════════════════════════════════════════════
# 5.  FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

def fig1_umap(mac):
    """UMAP coloured by sub-cluster + cancer type."""
    print("Fig 1: UMAP overview …")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sc.pl.umap(mac, color="mac_cluster", legend_loc="on data",
               title="Macrophage sub-clusters", ax=axes[0], show=False,
               palette=PALETTE)

    ct_col = "cancer_type_er" if "cancer_type_er" in mac.obs else \
             "cancer_type"    if "cancer_type"    in mac.obs else "sample_id"
    cats = mac.obs[ct_col].astype(str).unique().tolist()
    pal  = [CANCER_COLORS.get(c, "#999999") for c in
            sorted(cats, key=lambda x: (x not in CANCER_COLORS, x))]
    sc.pl.umap(mac, color=ct_col, title="Cancer type (BC split by ER status)",
               ax=axes[1], show=False, palette=pal)

    fig.suptitle("Macrophage Sub-clustering — GSE266330 (47 samples)",
                 fontsize=13, y=1.02)
    fig.subplots_adjust(wspace=0.4)
    save(fig, "01_umap_overview")


def fig2_cancer_type_bar(mac):
    """Stacked bar: fraction of each cancer type per cluster (BC split by ER)."""
    print("Fig 2: Cancer type bar …")
    ct_col = "cancer_type_er" if "cancer_type_er" in mac.obs.columns else "cancer_type"
    if ct_col not in mac.obs.columns:
        print("  Skipping: no cancer type column")
        return

    df = (mac.obs.groupby(["mac_cluster", ct_col], observed=True)
              .size().reset_index(name="n"))
    df["frac"] = df["n"] / df.groupby("mac_cluster")["n"].transform("sum")
    pivot = df.pivot(index="mac_cluster", columns=ct_col, values="frac").fillna(0)

    # Order columns: BC_ER+ first, then BC_ER-, then the rest alphabetically
    col_order = (
        [c for c in ["BC_ER+", "BC_ER-"] if c in pivot.columns] +
        sorted([c for c in pivot.columns if c not in ("BC_ER+", "BC_ER-")])
    )
    pivot = pivot[col_order]

    colors = [CANCER_COLORS.get(c, "#999999") for c in col_order]

    n_clusters = len(pivot)
    fig, ax = plt.subplots(figsize=(max(7, n_clusters * 0.85), 5))
    pivot.plot(kind="bar", stacked=True, ax=ax, color=colors, edgecolor="none")
    ax.set_xlabel("Macrophage sub-cluster", fontsize=11)
    ax.set_ylabel("Fraction", fontsize=11)
    ax.set_title("Cancer type composition per macrophage cluster\n(BC split by ESR1 expression)",
                 fontsize=11)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9,
              title="Cancer type")
    ax.tick_params(axis="x", rotation=0)
    save(fig, "02_cancer_type_bar_bc_er")


def fig3_esr1_sp1(mac):
    """ESR1 & SP1 feature plots on UMAP + dot plot across clusters."""
    print("Fig 3: ESR1 / SP1 …")
    genes = [g for g in ["ESR1", "SP1"] if g in mac.var_names]
    if not genes:
        print("  WARNING: neither ESR1 nor SP1 found in var_names")
        return

    # Feature plots
    n = len(genes)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, gene in zip(axes, genes):
        sc.pl.umap(mac, color=gene, ax=ax, show=False,
                   title=f"{gene} expression", color_map="Reds", vmin=0)
    fig.suptitle("ESR1 & SP1 expression — Macrophage sub-clusters", fontsize=12, y=1.02)
    save(fig, "03a_feature_plots_ESR1_SP1")

    # Dot plot per cluster
    n_cl = mac.obs["mac_cluster"].nunique()
    fig2, ax2 = plt.subplots(figsize=(max(6, n_cl * 0.8), 4))
    sc.pl.dotplot(mac, var_names=genes, groupby="mac_cluster",
                  ax=ax2, show=False,
                  title="ESR1 & SP1 across macrophage sub-clusters")
    save(fig2, "03b_dotplot_ESR1_SP1")

    # Identify ER+/SP1-high cluster
    sc.tl.score_genes(mac, gene_list=genes, score_name="ER_SP1_score")
    top_cl = mac.obs.groupby("mac_cluster")["ER_SP1_score"].mean().idxmax()
    print(f"  ER+/SP1-high cluster: {top_cl}")


# ═══════════════════════════════════════════════════════════════════════════════
# 6.  MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    if os.path.exists(CACHED):
        print("Loading cached macrophage AnnData …")
        mac = sc.read_h5ad(CACHED)
    else:
        # Extract h5ad from ZIP
        h5ad_path = extract_h5ad()

        print("Loading 47.integrated.h5ad …")
        adata = sc.read_h5ad(h5ad_path)
        print(f"  Full dataset: {adata.n_obs} cells, {adata.n_vars} genes")
        print(f"  obs columns: {list(adata.obs.columns)}")

        print("Isolating macrophages …")
        # Use celltype_C directly to avoid Unicode encoding ambiguity
        mphi_cat = None
        if "celltype_C" in adata.obs.columns:
            mphi_cat = next(
                (c for c in adata.obs["celltype_C"].cat.categories
                 if c.startswith("M") and len(c) == 2),
                None
            )
        if mphi_cat is not None:
            mask = adata.obs["celltype_C"] == mphi_cat
            mac  = adata[mask].copy()
            print(f"  Using celltype_C == {repr(mphi_cat)}: {mac.n_obs} cells")
        else:
            mac = isolate_macrophages(adata)

        # Direct column mapping for GSE266330
        for dst, src_col in [("cancer_type",   "cancer"),
                              ("sample_id",     "Seq.ID"),
                              ("condition",     "tissue.origin"),
                              ("gender",        "gender"),
                              ("tissue_origin", "tissue.origin")]:
            if src_col in mac.obs.columns:
                mac.obs[dst] = mac.obs[src_col].astype(str)
        print(f"  cancer_type: {mac.obs['cancer_type'].value_counts().to_dict()}")

        del adata
        print(f"  Macrophages: {mac.n_obs} cells")

        print("Sub-clustering …")
        mac = subcluster(mac)

        print("Harmonising metadata …")
        mac.obs = harmonise_metadata(mac.obs)

        print("Splitting BC by ER status …")
        mac = add_bc_er_split(mac)

        # Coerce object columns for h5ad write
        for col in mac.obs.columns:
            if mac.obs[col].dtype == object:
                mac.obs[col] = mac.obs[col].astype(str)

        mac.write_h5ad(CACHED)
        print(f"  Saved: {CACHED}")

    print(f"\nClusters: {sorted(mac.obs['mac_cluster'].unique())}")
    if "cancer_type_er" in mac.obs.columns:
        print(f"Cancer type (with BC split):\n{mac.obs['cancer_type_er'].value_counts().to_string()}")

    print("\nGenerating figures …")
    fig1_umap(mac)
    fig2_cancer_type_bar(mac)
    fig3_esr1_sp1(mac)

    print(f"\nAll done. Figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()
