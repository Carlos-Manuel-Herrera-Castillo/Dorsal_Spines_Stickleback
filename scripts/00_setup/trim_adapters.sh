#!/bin/bash

# Step 00b: Trim RNA-seq adapters with Trimmomatic.
# Original section: K* Trimming the adaptors
#
# Run on the cluster with:
#   sbatch scripts/00_setup/trim_adapters.sh

#SBATCH --job-name=Trimmomatic
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=06:00:00
#SBATCH --qos=6hours
#SBATCH --array=1-95
#SBATCH --output=/my/folder/path/Trimmomatic_output_%A_%a.o
#SBATCH --error=/my/folder/path/Trimmomatic_error_%A_%a.e
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=email@email.com

set -euo pipefail

module load Trimmomatic/0.39-Java-1.8

WORKDIR="/my/folder/path"
SAMPLE_IDS_FILE="Sample_ID.txt"
FASTQ_LIST_FILE="fastq_files.txt"

TRIMMOMATIC_JAR="/cluster/apps/Trimmomatic/0.39-Java-1.8/trimmomatic-0.39.jar"
ADAPTERS="/cluster/apps/Trimmomatic/0.39-Java-1.8/adapters/TruSeq3-PE.fa"

THREADS=4

cd "${WORKDIR}"

DATE=$(date '+%d-%m-%Y %H:%M:%S')
echo "Start of the job ${DATE}"

MYID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_IDS_FILE}")
echo "Processing sample: ${MYID}"

SAMPLE_OUTDIR="${WORKDIR}/${MYID}"
mkdir -p "${SAMPLE_OUTDIR}"

# fastq_files.txt should contain the locations of the files to be processed.
read1=$(grep "${MYID}" "${FASTQ_LIST_FILE}" | grep "R1")
read2=$(grep "${MYID}" "${FASTQ_LIST_FILE}" | grep "R2")

cd "${SAMPLE_OUTDIR}"

java -jar "${TRIMMOMATIC_JAR}" PE \
  -version \
  -threads "${THREADS}" \
  -phred33 \
  "${read1}" \
  "${read2}" \
  "${SAMPLE_OUTDIR}/${MYID}.1.paired.fastq" \
  "${SAMPLE_OUTDIR}/${MYID}.1.unpaired.fastq" \
  "${SAMPLE_OUTDIR}/${MYID}.2.paired.fastq" \
  "${SAMPLE_OUTDIR}/${MYID}.2.unpaired.fastq" \
  "ILLUMINACLIP:${ADAPTERS}:2:30:10:2" \
  SLIDINGWINDOW:4:15 \
  MINLEN:30

echo "Finished trimming sample: ${MYID}"
