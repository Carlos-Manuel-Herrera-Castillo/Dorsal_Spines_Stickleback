#!/usr/bin/env Rscript

# Step 04b: Fast FreqSum conversion for one chromosome.
# Original section: F* alternative speed-optimized FreqSum caller
#
# Run with:
#   Rscript scripts/04_freqsum/freqsum_fast.R
#
# This version is intended for a desktop machine and processes one chromosome.

rm(list = ls())

library(data.table)

input_file <- "C:/Users/bernerd/switchdrive/Institution/carlos/research/GWAS_DB/Basicpool.pileup.chrVI.txt"
output_file <- "C:/Users/daniel/switchdrive/Institution/carlos/research/GWAS_DB/Basicpool_freqSum_chrVI.txt"

# Alternative input/output from the original code:
# input_file <- "C:/Users/bernerd/switchdrive/Institution/carlos/research/GWAS_DB/Dpool.pileup.chrVI.txt"
# output_file <- "C:/Users/daniel/switchdrive/Institution/carlos/research/GWAS_DB/Dpool_freqSum_chrVI.txt"

message("Reading pileup file: ", input_file)
dt <- fread(
  input_file,
  header = FALSE,
  col.names = c("chromosome", "position", "nucleotide", "count")
)

dt[, count := as.integer(count)]

dt_wide <- dcast(
  dt,
  chromosome + position ~ nucleotide,
  value.var = "count",
  fun.aggregate = sum,
  fill = 0
)

nucleotides <- c("A", "C", "G", "T")
missing <- setdiff(nucleotides, names(dt_wide))

if (length(missing) > 0) {
  dt_wide[, (missing) := 0]
}

setcolorder(dt_wide, c("chromosome", "position", nucleotides))

fwrite(dt_wide, output_file, sep = "\t")

message("Finished fast FreqSum conversion")
message("Output file: ", output_file)
