#!/bin/bash

# Step 01b: Generate the STAR genome index.
# Original section: L* Preparing STAR
#
# Run on the cluster with:
#   sbatch scripts/01_alignment/star_index.sh

#SBATCH --job-name=STAR_GenomeGeneration_NCBI_Gac_v5_pitx1_ChrUnSep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=5:00:00
#SBATCH --qos=6hours
#SBATCH --array=1-1
#SBATCH --output=/my/folder/path/STAR_GenomeGeneration_V5Pitx1ChrUnSep.o
#SBATCH --error=/my/folder/path/STAR_GenomeGeneration_V5Pitx1ChrUnSep.e
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=email@email.com
#SBATCH --chdir=/my/folder/path

set -euo pipefail

module load STAR

THREADS=8
GENOME_DIR="/my/folder/path"
GENOME_FASTA="/my/folder/path/Gac_v5_pitx1_separateChrUn.fa"
ANNOTATION_GTF="/my/folder/path/NCBI_Gac_v5_pitx1.gtf"

STAR \
  --runMode genomeGenerate \
  --genomeSAindexNbases 13 \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${ANNOTATION_GTF}" \
  --sjdbGTFfeatureExon exon \
  --sjdbOverhang 100

echo "Finished STAR genome index"
echo "Genome index directory: ${GENOME_DIR}"
