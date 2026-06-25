#!/bin/bash
# Downsampled IPF Cell Atlas data (GSE136831, Adams et al. 2020 Science Advances)
# Source: Zenodo record 6878945 (Kaminski lab, Yale)
# .Robj files — load with load() in R

set -euo pipefail

OUTDIR="/home/qixin/jianting_project/ipf/GSE136831_zenodo"
mkdir -p "$OUTDIR"

echo "Downloading control.down.Robj (~315 MB)..."
wget -c -O "$OUTDIR/control.down.Robj" \
    "https://zenodo.org/records/6878945/files/control.down.Robj?download=1"

echo "Downloading ipf.down.Robj (~581 MB)..."
wget -c -O "$OUTDIR/ipf.down.Robj" \
    "https://zenodo.org/records/6878945/files/ipf.down.Robj?download=1"

echo "Done. Files:"
ls -lh "$OUTDIR"
