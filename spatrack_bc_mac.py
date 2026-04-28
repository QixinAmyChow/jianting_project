"""
SpaTrack trajectory analysis — BC macrophages
Input:  bm_analysis_out/spatrack/  (exported by analysis_update4.R Section 5)
Output: bm_analysis_out/figures_update_4/  (trajectory figures 17-21)

Run after analysis_update4.R:
  conda run -n scanpy_stable python spatrack_bc_mac.py
"""

import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

import spaTrack
from spaTrack.single_time import (
    get_ot_matrix,
    set_start_cells,
    get_ptime,
    get_velocity,
    get_velocity_grid,
    filter_gene,
)

# ── Paths ─────────────────────────────────────────────────────────────────────
SPATRACK_DIR = "bm_analysis_out/spatrack"
OUT_DIR      = "bm_analysis_out/figures_update_4"
os.makedirs(OUT_DIR, exist_ok=True)

DPI = 150

ARCH_COLORS = {
    "Mono":    "#ff7f0e",
    "MφOC":    "#d62728",
    "TregTex": "#9467bd",
    "ctrl":    "#aec7e8",
}

def save_fig(fig, name, w=8, h=6):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT_DIR, f"{name}.{ext}"),
                    dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {name}")

# ── Load exported data ────────────────────────────────────────────────────────
print("Loading exported BC macrophage data ...")

counts   = scipy.io.mmread(os.path.join(SPATRACK_DIR, "counts.mtx")).T.tocsr()
genes    = pd.read_csv(os.path.join(SPATRACK_DIR, "genes.txt"), header=None)[0].values
barcodes = pd.read_csv(os.path.join(SPATRACK_DIR, "barcodes.txt"), header=None)[0].values
meta     = pd.read_csv(os.path.join(SPATRACK_DIR, "metadata.csv"), index_col=0)
umap_df  = pd.read_csv(os.path.join(SPATRACK_DIR, "umap.csv"),     index_col=0)
pca_df   = pd.read_csv(os.path.join(SPATRACK_DIR, "pca.csv"),      index_col=0)

with open(os.path.join(SPATRACK_DIR, "start_cluster.txt")) as f:
    start_cluster = f.read().strip()

print(f"  Cells: {counts.shape[0]} | Genes: {counts.shape[1]}")
print(f"  Start cluster: {start_cluster}")

# Build AnnData — cells × genes
adata = ad.AnnData(
    X   = counts,
    obs = meta.loc[barcodes],
    var = pd.DataFrame(index=genes),
)
adata.obs_names = barcodes
adata.obsm["X_umap"] = umap_df.loc[barcodes].values.astype(float)
adata.obsm["X_pca"]  = pca_df.loc[barcodes].values.astype(float)

# SpaTrack needs adata.obs['cluster'] as strings
adata.obs["cluster"] = adata.obs["mac_subcluster"].astype(str)

print(f"  Cluster distribution:\n{adata.obs['cluster'].value_counts().sort_index()}")

# ── Step 1: Filter genes (SpaTrack recommended pre-processing) ────────────────
print("\n=== Step 1: Filter genes ===")
# SpaTrack's filter_gene keeps highly variable genes; data is already HVG-filtered
# so we skip redundant filtering and just normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(f"  Genes after log-norm: {adata.n_vars}")

# ── Step 2: Compute OT transition matrix ─────────────────────────────────────
print("\n=== Step 2: Optimal transport matrix ===")
print("  Computing OT (single-cell mode) ...")
Gs = get_ot_matrix(adata, data_type="single-cell", random_state=42)
from scipy.sparse import csr_matrix
adata.obsp["trans"] = csr_matrix(Gs)
print(f"  OT matrix shape: {Gs.shape}")

# ── Step 3: Set start cells ───────────────────────────────────────────────────
print("\n=== Step 3: Start cells ===")
if start_cluster in adata.obs["cluster"].values:
    start_cells = set_start_cells(adata, select_way="cell_type",
                                  cell_type=start_cluster)
else:
    # Fall back to cluster with most cells
    start_cluster = adata.obs["cluster"].value_counts().idxmax()
    print(f"  Start cluster not found — using {start_cluster}")
    start_cells = set_start_cells(adata, select_way="cell_type",
                                  cell_type=start_cluster)
print(f"  Start cells selected: {len(start_cells)}")

# ── Step 4: Compute pseudotime ────────────────────────────────────────────────
print("\n=== Step 4: Pseudotime ===")
ptime = get_ptime(adata, start_cells)
adata.obs["ptime"] = ptime
print(f"  Pseudotime range: [{ptime.min():.4f}, {ptime.max():.4f}]")

# Save pseudotime
pt_out = adata.obs[["cluster", "cancer_type_er", "archetype", "ptime"]].copy()
pt_out.to_csv(os.path.join(OUT_DIR, "bc_mac_spatrack_pseudotime.csv"))
print("  Saved: bc_mac_spatrack_pseudotime.csv")

# ── Step 5: Compute velocity ──────────────────────────────────────────────────
print("\n=== Step 5: Velocity ===")
try:
    P_grid, V_grid = get_velocity(adata, basis="umap",
                                  n_neigh_pos=10, n_neigh_gene=0,
                                  grid_num=40, smooth=0.5, density=1.0)
    velocity_ok = True
    print("  Velocity computed.")
except Exception as e:
    print(f"  Velocity ERROR: {e}")
    velocity_ok = False

# ── Plotting ──────────────────────────────────────────────────────────────────
print("\n=== Plots ===")
umap     = adata.obsm["X_umap"]
clusters = adata.obs["cluster"].values
archs    = adata.obs["archetype"].values
ptimes   = adata.obs["ptime"].values

# ── Plot 17: Pseudotime on UMAP ───────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
sc_pt = ax.scatter(umap[:, 0], umap[:, 1], c=ptimes, cmap="plasma",
                   s=8, alpha=0.8, linewidths=0)
plt.colorbar(sc_pt, ax=ax, label="SpaTrack Pseudotime")
ax.set_title(f"BC Macrophage Pseudotime (SpaTrack)\nStart: Cluster {start_cluster}",
             fontsize=13)
ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
ax.set_xticks([]); ax.set_yticks([])
save_fig(fig, "17_bc_spatrack_pseudotime")

# ── Plot 18: Cluster labels on UMAP ──────────────────────────────────────────
uniq_cl = sorted(set(clusters))
cl_cmap = plt.get_cmap("tab20", len(uniq_cl))
cl_lut  = {c: cl_cmap(i) for i, c in enumerate(uniq_cl)}

fig, ax = plt.subplots(figsize=(8, 6))
for cl in uniq_cl:
    mask = clusters == cl
    ax.scatter(umap[mask, 0], umap[mask, 1], s=8, alpha=0.8,
               color=cl_lut[cl], label=f"Cl {cl}", linewidths=0)
ax.legend(markerscale=2, fontsize=9, title="mac_subcluster",
          bbox_to_anchor=(1.01, 1), loc="upper left")
ax.set_title("BC Macrophage Clusters (SpaTrack input)", fontsize=13)
ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
save_fig(fig, "18_bc_spatrack_clusters")

# ── Plot 19: Archetype on UMAP ────────────────────────────────────────────────
uniq_arch = [a for a in ARCH_COLORS if a in archs]
fig, ax = plt.subplots(figsize=(8, 6))
for arch in uniq_arch:
    mask = archs == arch
    ax.scatter(umap[mask, 0], umap[mask, 1], s=8, alpha=0.8,
               color=ARCH_COLORS[arch], label=arch, linewidths=0)
other_mask = ~np.isin(archs, list(ARCH_COLORS.keys()))
if other_mask.any():
    ax.scatter(umap[other_mask, 0], umap[other_mask, 1], s=8, alpha=0.5,
               color="grey", label="other", linewidths=0)
ax.legend(markerscale=2, fontsize=10, title="Archetype",
          bbox_to_anchor=(1.01, 1), loc="upper left")
ax.set_title("BC Macrophage Archetype", fontsize=13)
ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
save_fig(fig, "19_bc_spatrack_archetype")

# ── Plot 20: Velocity stream on UMAP ─────────────────────────────────────────
if velocity_ok:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(umap[:, 0], umap[:, 1], c=ptimes, cmap="plasma",
               s=6, alpha=0.6, linewidths=0)
    ax.streamplot(P_grid[:, :, 0], P_grid[:, :, 1],
                  V_grid[:, :, 0], V_grid[:, :, 1],
                  color="black", linewidth=0.8, arrowsize=1.2,
                  density=1.2, broken_streamlines=True)
    ax.set_title(f"BC Macrophage Velocity (SpaTrack)\nStart: Cluster {start_cluster}",
                 fontsize=13)
    ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
    ax.set_xticks([]); ax.set_yticks([])
    save_fig(fig, "20_bc_spatrack_velocity_stream")

# ── Plot 21: Pseudotime per cluster (violin) ──────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
cl_order = sorted(set(clusters), key=lambda x: int(x) if x.isdigit() else x)
pt_by_cl = [ptimes[clusters == cl] for cl in cl_order]
parts = ax.violinplot(pt_by_cl, positions=range(len(cl_order)),
                      showmedians=True, showextrema=False)
for i, (pc, cl) in enumerate(zip(parts["bodies"], cl_order)):
    pc.set_facecolor(cl_lut[cl])
    pc.set_alpha(0.8)
parts["cmedians"].set_color("black")
ax.set_xticks(range(len(cl_order)))
ax.set_xticklabels([f"Cl {c}" for c in cl_order], rotation=45, ha="right")
ax.set_ylabel("Pseudotime")
ax.set_title("SpaTrack Pseudotime per BC Macrophage Cluster", fontsize=13)
ax.set_xlabel("mac_subcluster")
fig.tight_layout()
save_fig(fig, "21_bc_spatrack_pseudotime_violin")

# ── Plot 22: ER+/ER- pseudotime distribution ──────────────────────────────────
er_status = adata.obs["cancer_type_er"].values
er_groups = ["BC_ER+", "BC_ER-", "BC_ER?"]
er_present = [g for g in er_groups if g in er_status]

if len(er_present) >= 2:
    fig, ax = plt.subplots(figsize=(7, 5))
    er_colors = {"BC_ER+": "#d62728", "BC_ER-": "#1f77b4", "BC_ER?": "#fa9fb5"}
    pt_by_er = [ptimes[er_status == g] for g in er_present]
    parts2 = ax.violinplot(pt_by_er, positions=range(len(er_present)),
                           showmedians=True, showextrema=False)
    for pc, g in zip(parts2["bodies"], er_present):
        pc.set_facecolor(er_colors.get(g, "grey"))
        pc.set_alpha(0.85)
    parts2["cmedians"].set_color("black")
    ax.set_xticks(range(len(er_present)))
    ax.set_xticklabels(er_present)
    ax.set_ylabel("Pseudotime")
    ax.set_title("SpaTrack Pseudotime by ER Status (BC Macrophages)", fontsize=13)
    fig.tight_layout()
    save_fig(fig, "22_bc_spatrack_pseudotime_er_status")

print("\n===== spatrack_bc_mac.py complete =====")
print(f"All figures in: {OUT_DIR}")
print("Figures 17-22: SpaTrack pseudotime + velocity + violins")
