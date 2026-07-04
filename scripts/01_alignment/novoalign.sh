#!/bin/bash

# Step 01a: Align paired-end DNA reads with Novoalign.
# Original section: B* Alignment with novoalign
#
# Run from the folder containing the novoalign executable with:
#   bash scripts/01_alignment/novoalign.sh

set -euo pipefail

NOVOALIGN="./novoalign"
NOVOALIGN_INDEX="/path/to/reference/genome/Gac_v5_pitx1_separateChrUn.ndx"

READ1="/path/to/read/1/D_S_High.R1_part_aa.fastq"
READ2="/path/to/read/2/D_S_High.R2_part_aa.fastq"
OUTPUT_SAM="/path/for/output/D_S_High_part_aa_novoal.sam"

"${NOVOALIGN}" \
  -d "${NOVOALIGN_INDEX}" \
  -f "${READ1}" "${READ2}" \
  -F STDFQ \
  -t320 \
  -g 40 \
  -x 12 \
  -r N \
  -e 200 \
  -i PE 200,250 \
  -o SAM \
  -o FullNW \
  > "${OUTPUT_SAM}"

echo "Finished Novoalign alignment"
echo "Output SAM: ${OUTPUT_SAM}"
