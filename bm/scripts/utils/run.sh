#!/bin/bash
# BM Macrophage Pipeline — full run
# Usage: bash run.sh
# Logs: pipeline_update3.log  build_slides.log

set -e
cd "$(dirname "$0")"

echo "=== Step 1: R pipeline (figures) ==="
Rscript pipeline_update3.R 2>&1 | tee pipeline_update3.log

echo ""
echo "=== Step 2: Build slides ==="
python3 build_slides.py 2>&1 | tee build_slides.log

echo ""
echo "Done. Output: BM project-human data_update_3.pptx"
