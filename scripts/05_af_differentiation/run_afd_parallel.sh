#!/bin/bash

# Step 05a: Run AFD calculations for all chromosomes in parallel.
# Original section: H* AFD (parallelization)
#
# Run on the cluster with:
#   sbatch scripts/05_af_differentiation/run_afd_parallel.sh

#SBATCH --job-name=ADF_ParallelJob
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --qos=1day
#SBATCH --output=/my/folder/path/AFD_table_DS.o
#SBATCH --error=/my/folder/path/AFD_table_DS.e
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=email@email.com

set -euo pipefail

module load R

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPUTE_AFD="${SCRIPT_DIR}/compute_afd.R"

chromosomes=(
  "I" "II" "III" "IV" "V" "VI" "VII" "VIII" "IX" "X"
  "XI" "XII" "XIII" "XIV" "XV" "XVI" "XVII" "XVIII"
  "XIX" "XX" "XXI" "pitx1" "Y" "M"
)

locovta_values=(
  80 80 80 80 80
  80 80 80 80 80
  80 80 80 80 80
  80 80 80 50 80
  80 50 10 300
)

upcovta_values=(
  300 300 300 300 300
  300 300 300 300 300
  300 300 300 300 300
  300 300 300 250 300
  300 300 80 8800
)

locovtb_values=(
  80 80 80 80 80
  80 80 80 80 80
  80 80 80 80 80
  80 80 80 50 80
  80 50 10 300
)

upcovtb_values=(
  300 300 300 300 300
  300 300 300 300 300
  300 300 300 300 300
  300 300 300 250 300
  300 300 80 8800
)

num_cpus="${SLURM_CPUS_PER_TASK:-20}"
job_count=0

for i in "${!chromosomes[@]}"; do
  chr="${chromosomes[i]}"
  locovta="${locovta_values[i]}"
  upcovta="${upcovta_values[i]}"
  locovtb="${locovtb_values[i]}"
  upcovtb="${upcovtb_values[i]}"

  echo "Starting chromosome ${chr}"
  Rscript "${COMPUTE_AFD}" "${chr}" "${locovta}" "${upcovta}" "${locovtb}" "${upcovtb}" &

  job_count=$((job_count + 1))

  if (( job_count % num_cpus == 0 )); then
    wait
  fi
done

wait

echo "Finished all AFD jobs"
