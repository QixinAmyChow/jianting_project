#!/bin/bash
#SBATCH --job-name=zenodo_ctrl
#SBATCH --partition=longrunq
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/home/crnlqz/krenciklab/jianting_project/download_control_%j.log
#SBATCH --error=/home/crnlqz/krenciklab/jianting_project/download_control_%j.err

set -e
cd /home/crnlqz/krenciklab/jianting_project

RAW_DIR="bm_analysis_out/raw_geo"
ZIP="${RAW_DIR}/integrated_Seurat_objects.zip"
DEST_DIR="${RAW_DIR}/integrated_Seurat_objects/47.integrated_object_subset_by_major_celltypes"

echo "=== Step 1: Download integrated_Seurat_objects.zip (65 GB) ==="
mkdir -p "$RAW_DIR"
wget -c \
  "https://zenodo.org/records/16937964/files/integrated_Seurat_objects.zip?download=1" \
  -O "${ZIP}" \
  --tries=0 --timeout=120 --waitretry=60
echo "Download complete. Size: $(du -sh ${ZIP})"

echo ""
echo "=== Step 2: Find 'all control cells' RDS inside ZIP ==="
CTRL_INTERNAL=$(unzip -l "${ZIP}" | grep -i "control\|ctrl\|all.*ctrl\|ctrl.*all" | grep "\.rds" | awk '{print $NF}' | head -1)
if [ -z "$CTRL_INTERNAL" ]; then
  echo "ERROR: could not find control cells RDS by name. Listing all RDS files in ZIP:"
  unzip -l "${ZIP}" | grep "\.rds"
  exit 1
fi
echo "Found: ${CTRL_INTERNAL}"

echo ""
echo "=== Step 3: Extract control cells RDS ==="
mkdir -p "$DEST_DIR"
unzip -o "${ZIP}" "${CTRL_INTERNAL}" -d "${RAW_DIR}/"
echo "Extracted."

echo ""
echo "=== Step 4: Verify ==="
EXTRACTED="${RAW_DIR}/${CTRL_INTERNAL}"
if [ -f "$EXTRACTED" ]; then
  echo "OK: $(du -sh ${EXTRACTED})"
else
  echo "ERROR: extraction failed — file not found at ${EXTRACTED}"
  exit 1
fi

echo ""
echo "=== Step 5: Remove ZIP to free disk ==="
rm -f "${ZIP}"
echo "ZIP removed."

echo ""
echo "Done. Control cells RDS at: ${EXTRACTED}"
