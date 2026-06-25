#!/bin/bash
# Reyfman et al. 2019 (AJRCCM 199:1517-1536)
# GEO: GSE122960 — 8 donors + 4 IPF + 5 other ILD (17 samples)
# Only filtered Cell Ranger HDF5 files are downloaded (skip raw)

set -euo pipefail

OUTDIR="/home/qixin/jianting_project/ipf/GSE122960"
mkdir -p "$OUTDIR"

BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122960/suppl"

echo "Downloading GSE122960_RAW.tar (~720 MB)..."
wget -c -P "$OUTDIR" "$BASE/GSE122960_RAW.tar"

echo "Extracting filtered H5 files only..."
cd "$OUTDIR"
tar -xf GSE122960_RAW.tar --wildcards "*filtered_gene_bc_matrices_h5.h5"

echo "Listing extracted files:"
ls -lh *filtered*.h5 2>/dev/null | head -30

echo "Done. $(ls *filtered*.h5 2>/dev/null | wc -l) filtered H5 files extracted."
