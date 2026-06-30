#!/bin/bash
#SBATCH --job-name=zenodo_download
#SBATCH --partition=longrunq
#SBATCH --time=30-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/home/crnlqz/krenciklab/jianting_project/download_seurat_%j.log
#SBATCH --error=/home/crnlqz/krenciklab/jianting_project/download_seurat_%j.err

OUT=/home/crnlqz/krenciklab/jianting_project/bm_analysis_out/raw_geo/integrated_Seurat_objects.zip

wget -c \
  "https://zenodo.org/records/16937964/files/integrated_Seurat_objects.zip?download=1" \
  -O "${OUT}" \
  --tries=0 --timeout=120 --waitretry=60

echo "Exit code: $?"
echo "Final size: $(du -sh ${OUT})"
