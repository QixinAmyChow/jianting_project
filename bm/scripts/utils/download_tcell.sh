#!/bin/bash
#SBATCH --job-name=zenodo_tcell
#SBATCH --partition=longrunq
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home/crnlqz/krenciklab/jianting_project/download_tcell_%j.log
#SBATCH --error=/home/crnlqz/krenciklab/jianting_project/download_tcell_%j.err

set -e
cd /home/crnlqz/krenciklab/jianting_project

RAW_DIR="bm/data/raw/geo"
ZIP="${RAW_DIR}/integrated_Seurat_objects.zip"
TCELL_DEST="${RAW_DIR}/integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes/T_cell.rds"

# Skip if already extracted
if [ -f "$TCELL_DEST" ]; then
  echo "T_cell.rds already present at ${TCELL_DEST} — nothing to do."
  exit 0
fi

echo "=== Step 1: Download integrated_Seurat_objects.zip (65 GB) ==="
mkdir -p "$RAW_DIR"
wget -c \
  "https://zenodo.org/records/16937964/files/integrated_Seurat_objects.zip?download=1" \
  -O "${ZIP}" \
  --tries=0 --timeout=120 --waitretry=60
echo "Download complete. Size: $(du -sh ${ZIP})"

echo ""
echo "=== Step 2: Find T cell file inside ZIP ==="
TCELL_INTERNAL=$(unzip -l "${ZIP}" | grep -i "t_cell\|t-cell\|tcell" | grep "\.rds" | awk '{print $NF}' | head -1)
if [ -z "$TCELL_INTERNAL" ]; then
  echo "ERROR: could not find T cell RDS inside ZIP. Listing all RDS files:"
  unzip -l "${ZIP}" | grep "\.rds"
  exit 1
fi
echo "Found: ${TCELL_INTERNAL}"

echo ""
echo "=== Step 3: Extract T cell RDS ==="
unzip -o "${ZIP}" "${TCELL_INTERNAL}" -d "${RAW_DIR}/"
echo "Extracted to: ${RAW_DIR}/${TCELL_INTERNAL}"

echo ""
echo "=== Step 4: Verify ==="
EXTRACTED="${RAW_DIR}/${TCELL_INTERNAL}"
if [ -f "$EXTRACTED" ]; then
  echo "OK: $(du -sh ${EXTRACTED})"
else
  echo "ERROR: extraction failed"
  exit 1
fi

echo ""
echo "=== Step 5: Remove ZIP to free disk ==="
rm -f "${ZIP}"
echo "ZIP removed."

echo ""
echo "=== Step 6: Re-run analysis_update4.R to complete T cell section ==="
Rscript analysis_update4.R 2>&1 | tee analysis_update4_tcell_rerun.log

echo ""
echo "Done. T_cell.rds at: ${EXTRACTED}"
