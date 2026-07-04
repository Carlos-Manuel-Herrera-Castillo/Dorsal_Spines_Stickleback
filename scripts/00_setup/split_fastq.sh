#!/bin/bash

# Step 00a: Split a large FASTQ file into smaller parts.
# Original section: A* Split big .fastq files
#
# Run from the project folder with:
#   bash scripts/00_setup/split_fastq.sh

set -euo pipefail

INPUT_FASTQ="my/folder/path/D_S_High.R1.fastq"
OUTPUT_PREFIX="my/folder/path/D_S_High.R1_part_"

# The number of lines should be divisible by 4 because each FASTQ read uses 4 lines.
LINES_PER_FILE=360300000

split -l "${LINES_PER_FILE}" "${INPUT_FASTQ}" "${OUTPUT_PREFIX}"

echo "Finished splitting ${INPUT_FASTQ}"
echo "Output files start with: ${OUTPUT_PREFIX}"
