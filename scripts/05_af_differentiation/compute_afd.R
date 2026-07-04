#!/usr/bin/env Rscript

# Step 05b: Compute allele frequency difference (AFD) for one chromosome.
# Original section: I* AFD (R script named "compute_afd.R")
#
# Run manually with:
#   Rscript scripts/05_af_differentiation/compute_afd.R VI 80 300 80 300

rm(list = ls())

.libPaths("/library/path/R/x86_64-pc-linux-gnu-library/4.4")

library(data.table)

options(scipen = 50000000)

working_dir <- "/my/folder/path/"
setwd(working_dir)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: compute_afd.R <chr> <locovta> <upcovta> <locovtb> <upcovtb>")
}

chr <- args[1]
locovta <- as.numeric(args[2])
upcovta <- as.numeric(args[3])
locovtb <- as.numeric(args[4])
upcovtb <- as.numeric(args[5])

comp <- c("D_S_High_", "D_S_Low_")
popcovt <- 50
maft <- 0.1

outfile <- paste0(
  comp[1], ".", comp[2],
  ".afd.chr", chr,
  ".locovtHigh", locovta,
  ".upcovtHigh", upcovta,
  ".locovtHigh", locovtb,
  ".upcovtHigh", upcovtb,
  ".popcovt", popcovt,
  ".maft", maft,
  ".txt"
)

write.table(
  paste("chr", "pos", "majAl", "minAl", "afd", sep = " "),
  outfile,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

message("Reading FreqSum files for chromosome: ", chr)
a <- data.frame(fread(paste0(comp[1], ".FreqSumchr", chr, ".txt"), h = TRUE, stringsAsFactors = FALSE))
b <- data.frame(fread(paste0(comp[2], ".FreqSumchr", chr, ".txt"), h = TRUE, stringsAsFactors = FALSE))

a <- a[which(a$pos >= 1), ]
b <- b[which(b$pos >= 1), ]

a$AlleleCount <- rowSums(a[, 3:7])
b$AlleleCount <- rowSums(b[, 3:7])

uloc <- intersect(a$pos, b$pos)
a <- a[match(uloc, a$pos), ]
b <- b[match(uloc, b$pos), ]

pool <- a[, 3:7] + b[, 3:7]

maf <- function(x) {
  sorted <- sort(x, decreasing = TRUE)
  sorted[2] / sum(sorted[1:2]) >= maft
}

message("Filtering by MAF")
mafs <- apply(pool, MARGIN = 1, FUN = maf)
a <- a[mafs, ]
b <- b[mafs, ]

cova <- function(x) {
  total <- sum(sort(x, decreasing = TRUE)[1:2])
  total >= locovta && total <= upcovta
}

covb <- function(x) {
  total <- sum(sort(x, decreasing = TRUE)[1:2])
  total >= locovtb && total <= upcovtb
}

message("Filtering sample A by coverage")
covsa <- apply(a[, 3:7], MARGIN = 1, FUN = cova)
a <- a[covsa, ]

message("Filtering sample B by coverage")
covsb <- apply(b[, 3:7], MARGIN = 1, FUN = covb)
b <- b[covsb, ]

uloc <- intersect(a$pos, b$pos)
a <- a[match(uloc, a$pos), ]
b <- b[match(uloc, b$pos), ]

pool <- as.matrix(a[, 3:7] + b[, 3:7])

message("Writing AFD output: ", outfile)
for (i in seq_len(nrow(a))) {
  majmin <- names(sort(pool[i, ], decreasing = TRUE)[1:2])

  ai <- a[i, ]
  aimin <- ai[which(names(ai) == majmin[2])]
  aimaj <- ai[which(names(ai) == majmin[1])]

  if (sum(aimin, aimaj) >= popcovt) {
    bi <- b[i, ]
    bimin <- bi[which(names(bi) == majmin[2])]
    bimaj <- bi[which(names(bi) == majmin[1])]

    if (sum(bimin, bimaj) >= popcovt) {
      afd <- round(aimin / sum(aimin, aimaj) - bimin / sum(bimin, bimaj), 4)
      result <- paste(ai["chr"], ai["pos"], majmin[1], majmin[2], afd, sep = " ")

      write.table(
        result,
        outfile,
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        append = TRUE
      )
    }
  }
}

message("Finished AFD for chromosome: ", chr)
