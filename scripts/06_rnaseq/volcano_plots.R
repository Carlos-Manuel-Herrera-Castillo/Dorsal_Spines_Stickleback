#!/usr/bin/env Rscript

# Step 06b: Volcano plot and heatmap helpers for RNA-seq analyses.
# Original source: plotting parts of Q* Bulk RNAseq data analysis
#
# Use from another R script with:
#   source("scripts/06_rnaseq/volcano_plots.R")

library(EnhancedVolcano)
library(ggplot2)

plot_basic_volcano <- function(results_table, title, padj_cutoff, fc_cutoff,
                               label_size = 0) {
  EnhancedVolcano(
    results_table,
    title = title,
    lab = rownames(results_table),
    x = "log2FoldChange",
    y = "padj",
    labSize = label_size,
    pCutoff = padj_cutoff,
    FCcutoff = fc_cutoff
  )
}

add_peak_distance <- function(results_table, gene_ranges, midpoint_column = "midpt",
                              gene_column = "name", peak = 18159274) {
  midpoints <- vapply(
    rownames(results_table),
    function(gene) {
      match_index <- which(gene_ranges[[gene_column]] == gene)
      if (length(match_index) == 0) {
        return(NA_real_)
      }
      gene_ranges[[midpoint_column]][match_index[1]]
    },
    numeric(1)
  )

  results_table$midpoint <- midpoints
  results_table$peak_distance_scaled <- abs(results_table$midpoint - peak) /
    max(results_table$midpoint, na.rm = TRUE)

  results_table
}

plot_chrVI_distance_volcano <- function(results_table, title,
                                        color_breaks = 1000,
                                        xlim = c(-8, 8)) {
  color_fun <- colorRampPalette(c("#000000", "#d0d0d0", "#e1e1e1", "#ececec", "#f6f6f6"))
  colors <- color_fun(color_breaks)[
    cut(as.numeric(results_table$peak_distance_scaled), breaks = color_breaks)
  ]
  custom_colors <- setNames(colors, rownames(results_table))

  plot <- EnhancedVolcano(
    results_table,
    title = title,
    lab = NA,
    x = "log2FoldChange",
    y = "padj",
    xlim = xlim,
    colCustom = custom_colors,
    labSize = 1,
    pCutoff = 10,
    FCcutoff = 0,
    pointSize = 3,
    axisLabSize = 20,
    titleLabSize = 22
  )

  plot + guides(color = FALSE)
}

plot_gene_expression_by_age <- function(gene_name, expression_long) {
  ggplot(expression_long[expression_long$Gene == gene_name, ],
         aes(x = Age, y = value, colour = Origin)) +
    geom_boxplot() +
    theme_classic() +
    labs(title = gene_name) +
    ylab("Normalized Gene Expression (CPM)") +
    theme(plot.title = element_text(hjust = 0.5))
}
