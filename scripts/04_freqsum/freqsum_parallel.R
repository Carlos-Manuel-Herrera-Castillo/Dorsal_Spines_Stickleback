#!/usr/bin/env Rscript

# Step 04a: Convert pileup format to FreqSum format in parallel.
# Original section: F* FreqSum
#
# Run with:
#   Rscript scripts/04_freqsum/freqsum_parallel.R
#
# On a SLURM cluster, this script uses SLURM_CPUS_PER_TASK if available.

rm(list = ls())

.libPaths("/my/library/path/x86_64-pc-linux-gnu-library/4.4")

library(data.table)
library(parallel)

input_file <- "/my/folder/path/D_S_High_Merged.pileupAllChr.txt"
output_file <- "/my/folder/path/merged_D_S_Low.FreqSumAllChr.txt"

num_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

message("Reading pileup file: ", input_file)
data <- fread(
  input_file,
  header = FALSE,
  col.names = c("chr", "pos", "base", "count")
)

message("Using CPUs: ", num_cpus)
data_chunks <- split(data, cut(seq(nrow(data)), num_cpus, labels = FALSE))

process_chunk <- function(chunk) {
  chunk <- as.data.table(chunk)
  dcast(chunk, chr + pos ~ base, value.var = "count", fill = 0)
}

chunk_results <- mclapply(data_chunks, process_chunk, mc.cores = num_cpus)
transformed_data <- rbindlist(chunk_results)

unique_positions <- unique(data[, .(chr, pos)])
transformed_data <- merge(
  unique_positions,
  transformed_data,
  by = c("chr", "pos"),
  all.x = TRUE
)

transformed_data[is.na(transformed_data)] <- 0

fwrite(transformed_data, output_file, sep = "\t")

message("Finished FreqSum conversion")
message("Output file: ", output_file)
