#!/usr/bin/env Rscript

# Step 06c: Main bulk RNA-seq DESeq2 analysis.
# Original section: Q* Bulk RNAseq data analysis
#
# Run with:
#   Rscript scripts/06_rnaseq/deseq_analysis.R

rm(list = ls())

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(viridis)
library(reshape2)

setwd("my/folder/path")

source("scripts/06_rnaseq/pca_analysis.R")
source("scripts/06_rnaseq/volcano_plots.R")

gene_expression_matrix <- read.table("gene_expression_matrix.txt", header = TRUE)
Metadata <- read.csv("Metadata_for_gene_expression.csv")
Genes_in_ChrVI <- read.table("AnnotationChrVI_onlygenes.txt")
Genes_in_ChrUn <- read.table("AnnotationChrUnMappingToChrXIV_onlygenes.txt")

rownames(Metadata) <- Metadata$Sample

# Remove final rows containing summary statistics.
gene_expression_matrix <- head(gene_expression_matrix, -5)
rownames(gene_expression_matrix) <- gene_expression_matrix$Gene
gene_expression_matrix <- gene_expression_matrix[, -1]

# Keep genes expressed above 10 counts in more than 3 samples.
gene_expression_matrix <- gene_expression_matrix[rowSums(gene_expression_matrix > 10) > 3, ]

# Full model: check sample clustering by tissue, age, and origin.
deseq_matrix <- DESeqDataSetFromMatrix(
  countData = gene_expression_matrix,
  colData = Metadata,
  design = ~ Tissue + Days + Origin
)

dds <- DESeq(deseq_matrix)
vst <- varianceStabilizingTransformation(dds)

plotPCA.12(vst, intgroup = c("Tissue"))
plotPCA.12(vst, intgroup = c("Days"))
plotPCA.12(vst, intgroup = c("Origin"))

res_tissue <- results(dds, name = "Tissue_PS_vs_DS")
res_tissue <- na.omit(res_tissue)
res_tissue$log2FoldChange <- as.numeric(res_tissue$log2FoldChange)
res_tissue_sig <- res_tissue[res_tissue$padj < 0.05, ]

plot_basic_volcano(
  res_tissue,
  title = "Tissue: pelvic spine vs dorsal spine",
  padj_cutoff = 0.05,
  fc_cutoff = 1
)

# Dorsal spine only: compare origin while accounting for days.
Metadata_DS <- Metadata[Metadata$Tissue == "DS", ]
gene_expression_matrix_DS <- gene_expression_matrix[
  ,
  colnames(gene_expression_matrix) %in% Metadata_DS$Sample
]
gene_expression_matrix_DS <- gene_expression_matrix_DS[
  rowSums(gene_expression_matrix_DS > 10) > 3,
]

deseq_matrix_DS <- DESeqDataSetFromMatrix(
  countData = gene_expression_matrix_DS,
  colData = Metadata_DS,
  design = ~ Days + Origin
)

dds_DS <- DESeq(deseq_matrix_DS)
vst_DS <- varianceStabilizingTransformation(dds_DS)

plotPCA.12(vst_DS, intgroup = c("Origin"))
plotPCA.23(vst_DS, intgroup = c("Origin"))
plotPCA.12(vst_DS, intgroup = c("Days"))
plotPCA.23(vst_DS, intgroup = c("Days"))
plotPCA.14(vst_DS, intgroup = c("Days"))
plotPCA.15(vst_DS, intgroup = c("Days"))

res_tissue_DS <- results(dds_DS, name = "Origin_SCAD_vs_DUIN")
res_tissue_DS <- na.omit(res_tissue_DS)
res_tissue_DS$log2FoldChange <- as.numeric(res_tissue_DS$log2FoldChange)
res_tissue_DS_sig <- res_tissue_DS[
  res_tissue_DS$padj < 10e-2 & abs(res_tissue_DS$log2FoldChange) > 1,
]

plot_basic_volcano(
  res_tissue_DS,
  title = "Dorsal spine: present vs absent",
  padj_cutoff = 10e-2,
  fc_cutoff = 1
)

res_tissue_DS_chrVI <- res_tissue_DS[
  rownames(res_tissue_DS) %in% as.character(Genes_in_ChrVI$V1),
]

plot_basic_volcano(
  res_tissue_DS_chrVI,
  title = "Dorsal spine: present vs absent in Chr VI",
  padj_cutoff = 10e-10,
  fc_cutoff = 1.5
)

matching_genes_ChrVI <- Genes_in_ChrVI$V1 %in% rownames(res_tissue_DS_sig)
genes_in_both_ChrVI <- Genes_in_ChrVI$V1[matching_genes_ChrVI]
matching_row_names_ChrVI <- rownames(res_tissue_DS_sig)[
  rownames(res_tissue_DS_sig) %in% genes_in_both_ChrVI
]
print(matching_row_names_ChrVI)

# Heatmap for strongly differentiated genes.
annotation_columns <- Metadata_DS[, c(3:4)]
column_colors <- list(
  Origin = c(DUIN = "#f3766e", SCAD = "#1cbdc2"),
  Days = c(
    "10days" = "#fde725",
    "13days" = "#35b779",
    "16days" = "#31688e",
    "19days" = "#440154"
  )
)

heat_plot <- pheatmap(
  assay(vst_DS)[res_tissue_DS$padj < 1e-5 & abs(res_tissue_DS$log2FoldChange) > 1, ],
  col = inferno(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D",
  scale = "row",
  annotation_col = annotation_columns,
  annotation_colors = column_colors,
  annotation_names_col = TRUE,
  fontsize_row = 10,
  fontsize_col = 7,
  angle_col = 45,
  legend_breaks = c(-2, 0, 2),
  legend_labels = c("Low", "Medium", "High"),
  show_colnames = TRUE,
  show_rownames = FALSE,
  main = "Heatmap dorsal spine"
)

write.csv(res_tissue_DS, "res_tissue_DS.csv", row.names = TRUE)

# Optional: gene expression boxplot data for DS samples.
gene_expression_matrix_cpm <- t(t(gene_expression_matrix_DS) / colSums(gene_expression_matrix_DS) * 1000000)
melt_gene_exp <- reshape2::melt(gene_expression_matrix_cpm)
colnames(melt_gene_exp) <- c("Gene", "Sample", "value")
melt_gene_exp$Age <- with(Metadata, Days[match(melt_gene_exp$Sample, rownames(Metadata))])
melt_gene_exp$Origin <- with(Metadata, Origin[match(melt_gene_exp$Sample, rownames(Metadata))])

plot_gene_expression_by_age("lnpk", melt_gene_exp)
