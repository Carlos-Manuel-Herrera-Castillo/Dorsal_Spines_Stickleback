#!/bin/bash

# Step 02b: Merge multiple BAM files into one BAM file.
# Original section: D* Merge .bam files
#
# Run with:
#   bash scripts/02_bam_processing/merge_bam.sh

set -euo pipefail

OUTPUT_BAM="D_S_High_Merged.bam"
INPUT_BAM_GLOB="/my/folder/path/D_S_High_part_a*_novoal.bam"

samtools merge "${OUTPUT_BAM}" ${INPUT_BAM_GLOB}

echo "Finished merging BAM files"
echo "Output BAM: ${OUTPUT_BAM}"
