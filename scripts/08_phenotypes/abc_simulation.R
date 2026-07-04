#!/usr/bin/env Rscript

# Step 08b: Simulate F1/F2 phenotypes and estimate SCAD allele frequency with ABC.
# Original section: T* Simulation combined with Approximate Bayesian Computation
#
# Run with:
#   Rscript scripts/08_phenotypes/abc_simulation.R

rm(list = ls())

library(abc)

setwd("C:/Users/bernerd/switchdrive/Institution/carlos/research/DSpheno_DB")

n_simulations <- 10000
nfo <- 20
nfam <- 5
fas <- c(6, 11, 13, 13, 17)
nf1 <- 60
nf2 <- 530

simdat <- NULL

for (j in seq_len(n_simulations)) {
  p <- runif(1, min = 0, max = 1)

  p1 <- rep("aa", nfo)
  p2 <- c(
    rep("AA", round(nfo * p^2)),
    rep("Aa", round(nfo * 2 * p * (1 - p))),
    rep("aa", round(nfo * (1 - p)^2))
  )
  p2 <- sample(p2, size = nfo, replace = TRUE)

  f1 <- NULL
  for (i in seq_len(nfam)) {
    parent1 <- sample(p1, 1)
    parent2 <- sample(p2, 1)

    gam1 <- c("a", "a")
    if (parent2 == "aa") gam2 <- c("a", "a")
    if (parent2 == "Aa") gam2 <- c("A", "a")
    if (parent2 == "AA") gam2 <- c("A", "A")

    for (k in seq_len(fas[i])) {
      g1 <- sample(gam1, 1)
      g2 <- sample(gam2, 1)
      f1 <- c(f1, paste0(sort(c(g1, g2), decreasing = TRUE), collapse = ""))
    }
  }

  f1pa <- matrix(sample(f1, 2 * nf2, replace = TRUE), ncol = 2, byrow = TRUE)
  f2 <- character(nf2)

  for (i in seq_len(nf2)) {
    parent1 <- f1pa[i, 1]
    parent2 <- f1pa[i, 2]

    gam1 <- if (parent1 == "aa") c("a", "a") else c("A", "a")
    gam2 <- if (parent2 == "aa") c("a", "a") else c("A", "a")

    g1 <- sample(gam1, 1)
    g2 <- sample(gam2, 1)
    f2[i] <- paste0(sort(c(g1, g2), decreasing = TRUE), collapse = "")
  }

  output <- c(
    round(p, 3),
    round(length(which(f1 == "aa")) / nf1, 3),
    round(length(which(f2 == "aa")) / nf2, 3)
  )
  simdat <- rbind(simdat, output)
}

simdat <- rbind(c("AF.SCAD", "F1.2DS", "F2.2DS"), simdat)
write.table(
  simdat,
  "DSsimusForABC_n10000_HWstoch_nfam.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

d <- read.table("DSsimusForABC_n10000_HWstoch_nfam.txt", h = TRUE)

obs <- c(0.5, 0.68)
pars <- data.frame(d[, 1])
names(pars) <- "AF.SCAD"
props <- d[, 2:3]

estim <- NULL
for (i in seq_len(500)) {
  eval <- abc(target = obs, param = pars, sumstat = props, tol = 0.01, method = "neuralnet")
  summary_eval <- round(summary(eval), 2)
  iestim <- c(summary_eval[3, 1], summary_eval[4, 1])
  estim <- rbind(estim, iestim)
}

estim <- data.frame(estim)
names(estim) <- c("medianEst", "meanEst")
write.table(estim, "ABCestim_stoch_nfam_tol0.01.txt", quote = FALSE, row.names = FALSE)

hist(estim$medianEst, breaks = 100, xlim = c(0, 1))
print(median(estim$medianEst))
