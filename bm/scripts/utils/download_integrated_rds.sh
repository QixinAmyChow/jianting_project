#!/bin/bash
# Download integrated_Seurat_objects.zip from Zenodo (Liu et al. 2025, Cell Genomics)
# Zenodo record: 10.5281/zenodo.14270976 (latest version: 16937964)

OUTDIR="/home/qixin/jianting_project/bm/data/raw/geo"
URL="https://zenodo.org/api/records/16937964/files/integrated_Seurat_objects.zip/content"
OUTFILE="$OUTDIR/integrated_Seurat_objects.zip"
LOGFILE="/home/qixin/jianting_project/download_integrated_rds.log"

mkdir -p "$OUTDIR"

echo "=== Download started: $(date) ===" | tee -a "$LOGFILE"
echo "Target: $OUTFILE" | tee -a "$LOGFILE"
echo "Expected size: ~69.3 GB" | tee -a "$LOGFILE"

wget --continue \
     --progress=dot:giga \
     --timeout=60 \
     --tries=10 \
     --waitretry=30 \
     -O "$OUTFILE" \
     "$URL" 2>&1 | tee -a "$LOGFILE"

EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 0 ]; then
    echo "=== Download complete: $(date) ===" | tee -a "$LOGFILE"
    echo "Final size: $(du -sh "$OUTFILE" | cut -f1)  $OUTFILE" | tee -a "$LOGFILE"
    echo "Unzipping into $OUTDIR/integrated_Seurat_objects/ ..." | tee -a "$LOGFILE"
    unzip -o "$OUTFILE" -d "$OUTDIR/" 2>&1 | tee -a "$LOGFILE"
    echo "=== Unzip complete: $(date) ===" | tee -a "$LOGFILE"
else
    echo "=== Download FAILED (exit code $EXIT_CODE): $(date) ===" | tee -a "$LOGFILE"
    exit $EXIT_CODE
fi
