# Jianting BM scRNA-seq Project

## User context
Bioinformatician/computational support for Jianting Sheng at Houston Methodist Neal Cancer Center. Works on single-cell RNA-seq analysis (Seurat in R, Scanpy/Python), GSEA, and presentation generation (python-pptx). Comfortable with both R and Python pipelines for scRNA-seq.

## Project overview
Bone marrow (BM) macrophage sub-clustering analysis for Jianting Sheng, based on Liu et al. 2025 (Cell Genomics) dataset, GEO: GSE266330 — 42 human BM samples, 8 cancer types, ~180K cells. Manuscript in preparation (`8-6-2025 BM manuscript SW.docx`).

**Key finding:** At resolution=0.3 (nn=10, min.dist=0.05), Cluster 1 = 94.8% BC_ER+ cells (1791/1889), making it the dominant ER+ macrophage sub-population. BC08 = 54% of all BC_ER+ cells (single-patient dominance caveat).

## Analysis iterations completed
- Update 1: initial UMAP, cancer type bar, ESR1/SP1 feature plots — `BM project-human data_update.pptx`
- Update 2: added DEG heatmap/dotplot, GSEA (Hallmark + cancer) — `BM project-human data_update_2.pptx`
- Update 3 (most complete, Apr 21 2026): parameter sweep grid, sanity checks, broader GSEA (GO:BP + Reactome + C7), Cluster 1 hypothesis grouping — `BM project-human data_update_3.pptx` + `BM project-human data_update_4.pptx`

## Update 3 GSEA findings for Cluster 1 (BC_ER+)
- Hallmark/C6: no significant terms (all padj > 0.15); trending DOWN: IL-6/JAK/STAT3, IFN-g, KRAS-breast
- Broader (GO:BP + Reactome + C7): 21 deduplicated terms — mac polarization (PPARG-KO, WBP7-KO, MDSC DOWN), cytokine/IFN-g DOWN, TGF-b/AKT DOWN, T cell activation DOWN, Wnt/osteogenic DOWN, ECM weakly UP
- Interpretation: BC_ER+ macrophages appear transcriptionally suppressed / immunosuppressive TME

## Key scripts (last modified Apr 21 2026)
- `gsea_cl1_compare_save.R` — persists GSEA selections for slide 12 and slide 14 to `figures_update_3/gsea_slide12_14_selections.rds`
- `gsea_families_shared.R` — shared semantic deduplication families
- `gsea_cl1_combined.R` — Cluster 1 combined GSEA figure
- `gsea_alldb_stacking.R` — stacked GSEA figure (05e) all databases
- `embed_figures_update3.py` — assembles final pptx from `figures_update_3/`
- `regen_clustering_figs_update3.R` — regenerates figs 01–04 from `mac_update3_best.rds`

## Key data file
`bm_analysis_out/figures_update_3/mac_update3_best.rds` — Seurat object, mac subset, best clustering params applied (excluded from git, regenerable)

## Next steps
- DEG: Cluster 1 vs BC_ER- clusters
- Validate SP1/ESR1 co-expression in Cluster 1
- Link to S2C2 crosstalk model (ER+ Mφ → CD8+ T cell axis)
- Validate in mouse scRNA (Jianting W2/W3)

When suggesting next analyses, prioritize Cluster 1 (BC_ER+) characterization and the S2C2 crosstalk angle. Be aware of the BC08 single-patient dominance caveat.
