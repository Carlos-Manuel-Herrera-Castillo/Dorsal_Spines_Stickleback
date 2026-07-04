# Step 03: Generate pileup counts from a BAM file.
# Original section: E* Pileup (in R)
#
# Run with:
#   Rscript scripts/03_pileup/pileup.R

rm(list = ls())

library(Rsamtools)
library(GenomicRanges)

working_dir <- "/my/folder/path/"
sample_name <- "D_S_High_Merged"
chromosomes_file <- "/my/folder/path/chrs.txt"

setwd(working_dir)

bam <- paste0(sample_name, ".bam")
outfile <- paste0(sample_name, ".pileupAllChr.txt")

pileup_params <- PileupParam(
  max_depth = 5000,
  distinguish_strands = FALSE,
  min_nucleotide_depth = 0,
  min_base_quality = 20,
  ignore_query_Ns = TRUE,
  include_deletions = TRUE
)

chromosomes <- readLines(chromosomes_file)

for (chr in chromosomes) {
  message("Processing chromosome: ", chr)

  scan_params <- ScanBamParam(
    which = GRanges(chr, IRanges(start = 1, end = 40000000))
  )

  pileup_data <- pileup(
    bam,
    scanBamParam = scan_params,
    pileupParam = pileup_params
  )

  pileup_data <- pileup_data[, -5]

  write.table(
    pileup_data,
    outfile,
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    append = TRUE
  )
}

message("Finished pileup")
message("Output file: ", outfile)
