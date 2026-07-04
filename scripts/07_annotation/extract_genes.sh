#!/bin/bash

# Step 07: Extract gene IDs from one chromosome in a GFF file.
# Original section: N* Extract gene names from chromosome VI
#
# Run with:
#   bash scripts/07_annotation/extract_genes.sh

set -euo pipefail

gff_file="NCBI_Gac_v5_pitx1.gff"
chromosome="chrVI"
output_file="AnnotationChrVI_onlygenes.txt"

awk -v chr="${chromosome}" '$1 == chr && $3 == "gene"' "${gff_file}" | \
  grep -o 'ID=gene-[^;]*' | \
  cut -d'=' -f2 > "${output_file}"

echo "Gene names extracted and saved to ${output_file}"
