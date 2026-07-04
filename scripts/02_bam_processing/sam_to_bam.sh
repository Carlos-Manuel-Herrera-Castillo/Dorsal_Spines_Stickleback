#!/bin/bash

# Step 02a: Convert a SAM alignment file to BAM.
# Original section: C* .sam to .bam conversion using samtools
#
# Run with:
#   bash scripts/02_bam_processing/sam_to_bam.sh

set -euo pipefail

INPUT_SAM="D_S_High_part_aa_novoal.sam"
OUTPUT_BAM="D_S_High_part_aa_novoal.bam"

samtools view -S -b "${INPUT_SAM}" > "${OUTPUT_BAM}"

echo "Finished SAM to BAM conversion"
echo "Input SAM: ${INPUT_SAM}"
echo "Output BAM: ${OUTPUT_BAM}"
