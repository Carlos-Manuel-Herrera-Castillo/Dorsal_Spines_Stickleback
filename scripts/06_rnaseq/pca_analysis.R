#!/usr/bin/env Rscript

# Step 06a: PCA helper functions for RNA-seq analyses.
# Original section: P* Principal Component Analysis (PCA) conditions
#
# Use from another R script with:
#   source("scripts/06_rnaseq/pca_analysis.R")

library(ggplot2)
library(ggrepel)
library(matrixStats)

plot_pca_pair <- function(object, pc_x = 1, pc_y = 2, intgroup = "condition",
                          ntop = 5000, returnData = FALSE) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup_df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup_df, 1, paste, collapse = " : "))
  } else {
    colData(object)[[intgroup]]
  }

  pc_x_name <- paste0("PC", pc_x)
  pc_y_name <- paste0("PC", pc_y)

  d <- data.frame(
    x = pca$x[, pc_x],
    y = pca$x[, pc_y],
    group = group,
    intgroup_df,
    name = colData(object)[, 1]
  )

  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(pc_x, pc_y)]
    return(d)
  }

  ggplot(d, aes(x = x, y = y, color = group, label = name)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3) +
    xlab(paste0(pc_x_name, ": ", round(percentVar[pc_x] * 100), "% variance")) +
    ylab(paste0(pc_y_name, ": ", round(percentVar[pc_y] * 100), "% variance")) +
    coord_fixed()
}

plotPCA.12 <- function(object, intgroup = "condition", ntop = 5000, returnData = FALSE) {
  plot_pca_pair(object, 1, 2, intgroup, ntop, returnData)
}

plotPCA.13 <- function(object, intgroup = "condition", ntop = 5000, returnData = FALSE) {
  plot_pca_pair(object, 1, 3, intgroup, ntop, returnData)
}

plotPCA.23 <- function(object, intgroup = "condition", ntop = 5000, returnData = FALSE) {
  plot_pca_pair(object, 2, 3, intgroup, ntop, returnData)
}

plotPCA.14 <- function(object, intgroup = "condition", ntop = 5000, returnData = FALSE) {
  plot_pca_pair(object, 1, 4, intgroup, ntop, returnData)
}

plotPCA.15 <- function(object, intgroup = "condition", ntop = 5000, returnData = FALSE) {
  plot_pca_pair(object, 1, 5, intgroup, ntop, returnData)
}
