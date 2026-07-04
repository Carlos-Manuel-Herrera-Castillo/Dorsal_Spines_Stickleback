#!/bin/bash

# Step 01c: Align trimmed RNA-seq reads with STAR.
# Original section: M* STAR alignment
#
# Run on the cluster with:
#   sbatch scripts/01_alignment/star_align.sh

#SBATCH --job-name=STAR_aligning_NCBI_Gac_v5_pitx1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=18G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --array=1-95
#SBATCH --output=/my/folder/path/STAR_mapping_50_Trimmomatic_NCBI_Gac_v5_pitx1_output_%A_%a.o
#SBATCH --error=/my/folder/path/STAR_mapping_50_Trimmomatic_NCBI_Gac_v5_pitx1_error_%A_%a.e
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=email@email.com
#SBATCH --chdir=/my/folder/path

set -euo pipefail

module load STAR/2.7.3a-foss-2018b

WORKDIR="/my/folder/path"
SAMPLE_IDS_FILE="Sample_ID.txt"
FASTQ_LIST_FILE="list_of_fastq.txt"
STAR_SAMPLE_OUTPUT_BASE="/scicore/home/bernerd/herrer0002/RNAseq/STAR"
GENOME_DIR="/my/folder/path"
THREADS=8

cd "${WORKDIR}"

DATE=$(date '+%d-%m-%Y %H:%M:%S')
echo "Start of the job ${DATE}"

MYID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_IDS_FILE}")
echo "Processing sample: ${MYID}"

SAMPLE_OUTDIR="${STAR_SAMPLE_OUTPUT_BASE}/${MYID}"
mkdir -p "${SAMPLE_OUTDIR}"

# list_of_fastq.txt should contain paired forward and reverse trimmed FASTQ files.
read1=$(grep "${MYID}" "${FASTQ_LIST_FILE}" | grep "1.paired")
read2=$(grep "${MYID}" "${FASTQ_LIST_FILE}" | grep "2.paired")

cd "${SAMPLE_OUTDIR}"

STAR \
  --outFilterMultimapNmax 1 \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${read1}" "${read2}" \
  --outSAMtype BAM SortedByCoordinate \
  --outReadsUnmapped Fastx \
  --runThreadN "${THREADS}" \
  --limitBAMsortRAM 40000000000 \
  --outFileNamePrefix "${WORKDIR}/${MYID}/${MYID}.STAR"

echo "Finished STAR alignment for sample: ${MYID}"
