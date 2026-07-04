#!/bin/bash

# Step 04c: Split one all-chromosome FreqSum file into separate chromosome files.
# Original section: G* FreqSum splitting by chromosome
#
# Run with:
#   bash scripts/04_freqsum/split_by_chr.sh

set -euo pipefail

input_file="/my/folder/path/merged_D_S_Low.FreqSumAllChr.txt"
output_dir="/my/folder/path/DS_Separate_Chr_Novoalign"
output_prefix="D_S_Low_BWA.FreqSum"

mkdir -p "${output_dir}"

chromosomes=(
  "chrI" "chrII" "chrIII" "chrIV" "chrV"
  "chrVI" "chrVII" "chrVIII" "chrIX" "chrX"
  "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV"
  "chrXVI" "chrXVII" "chrXVIII" "chrXIX" "chrXX"
  "chrXXI" "pitx1" "chrY" "chrM"
)

for chr in "${chromosomes[@]}"; do
  output_file="${output_dir}/${output_prefix}${chr}.txt"
  echo "Writing ${output_file}"

  echo -e "chr\tpos\tA\tC\tG\tT\t-" > "${output_file}"
  grep "^${chr}[[:space:]]" "${input_file}" >> "${output_file}" || true
done

echo "Finished splitting FreqSum file by chromosome"
