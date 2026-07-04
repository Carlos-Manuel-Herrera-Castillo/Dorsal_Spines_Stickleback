#!/usr/bin/env Rscript

# Step 05c: Fast AFD calculation for chromosome VI.
# Original section: I* speed-optimized AFD caller
#
# Run with:
#   Rscript scripts/05_af_differentiation/afd_fast.R

rm(list = ls())

library(data.table)
library(Rfast)

maft <- 0.25
bcovt <- c(35, 110)
dcovt <- c(100, 240)

working_dir <- "C:/Users/bernerd/switchdrive/Institution/carlos/research/GWAS_DB"
setwd(working_dir)

outfile <- paste0("afd.bpool.dpool.chrVI.maft", maft, ".txt")

b <- fread("Basicpool_freqSum_chrVI.txt")
d <- fread("Dpool_freqSum_chrVI.txt")

header <- paste(
  "chr", "pos",
  "bpool_A", "bpool_C", "bpool_G", "bpool_T", "bpool_DEL", "bpool_totBS",
  "dpool_A", "dpool_C", "dpool_G", "dpool_T", "dpool_DEL", "dpool_totBS",
  "afd",
  sep = " "
)

write.table(header, outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)

b[, depth := A + C + G + T]
d[, depth := A + C + G + T]

bf <- b[depth >= bcovt[1] & depth <= bcovt[2]]
df <- d[depth >= dcovt[1] & depth <= dcovt[2]]

shared_positions <- intersect(bf$position, df$position)

bs <- bf[position %in% shared_positions][order(position)]
ds <- df[position %in% shared_positions][order(position)]

message("Basic pool retained proportion: ", nrow(bf) / nrow(b))
message("D pool retained proportion: ", nrow(df) / nrow(d))

rm("b", "d", "bf", "df")
gc()

bfr <- bs[, .(A = A / depth, C = C / depth, G = G / depth, T = T / depth)]
dfr <- ds[, .(A = A / depth, C = C / depth, G = G / depth, T = T / depth)]

paf <- (bfr + dfr) / 2
mat <- as.matrix(paf)

srt <- Rfast::rowSort(mat, descending = TRUE)
maf <- srt[, 2]
idx <- which(maf >= maft)

bs <- bs[idx, ]
ds <- ds[idx, ]
bfr <- bfr[idx, ]
dfr <- dfr[idx, ]
paf <- paf[idx, ]

message("Writing fast AFD output: ", outfile)
for (i in seq_len(nrow(paf))) {
  majal <- as.integer(which.max(paf[i, ]))
  afd <- round(abs(as.matrix(bfr)[i, majal] - as.matrix(dfr)[i, majal]), 4)

  out <- paste(
    paste(bs[i, ], collapse = " "),
    paste(ds[i, -c(1:2)], collapse = " "),
    afd,
    sep = " "
  )

  write.table(out, outfile, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

message("Finished fast AFD")
