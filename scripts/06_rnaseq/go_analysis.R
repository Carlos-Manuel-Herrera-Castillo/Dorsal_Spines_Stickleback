#!/usr/bin/env Rscript

# Step 06d: GO term analysis.
# Original section: R* GO term analysis
#
# Run with:
#   Rscript scripts/06_rnaseq/go_analysis.R

rm(list = ls())

library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)

setwd("/my/folder/path")

DS_genes <- read.csv("GenesDS.csv", row.names = 1)
GaToDr <- read.csv("ENSAMBL_IDs_Gac_to_Dre.csv")
NCBItoENSEMBL <- read.csv("NCBI_SymbolToENSEMBLid.csv")

universe_genes <- rownames(DS_genes)
genes_to_test <- rownames(DS_genes[abs(DS_genes$log2FoldChange) > 1, ])

# Convert differentially expressed stickleback gene symbols to zebrafish ENSEMBL IDs.
matching_genes <- intersect(NCBItoENSEMBL$Symbol, genes_to_test)
non_matching_genes <- setdiff(union(NCBItoENSEMBL$Symbol, genes_to_test), matching_genes)

cat("Number of matching genes:", length(matching_genes), "\n")
cat("Number of non-matching genes:", length(non_matching_genes), "\n")

filtered_NCBItoENSEMBL <- NCBItoENSEMBL[NCBItoENSEMBL$Symbol %in% matching_genes, ]
filtered_GaToDr <- GaToDr[
  GaToDr$Gene.stable.ID %in% filtered_NCBItoENSEMBL$Ensembl.GeneIDs,
]

filtered_GaToDr_veryconservative <- filtered_GaToDr[
  filtered_GaToDr$Zebrafish.homology.type == "ortholog_one2one" &
    filtered_GaToDr$Zebrafish.orthology.confidence..0.low..1.high. == 1,
]
genes_to_test <- unique(filtered_GaToDr_veryconservative$Zebrafish.gene.stable.ID)

# Convert all expressed genes to define the GO universe.
matching_genes_universe <- intersect(NCBItoENSEMBL$Symbol, universe_genes)
non_matching_genes_universe <- setdiff(
  union(NCBItoENSEMBL$Symbol, universe_genes),
  matching_genes_universe
)

cat("Number of matching universe genes:", length(matching_genes_universe), "\n")
cat("Number of non-matching universe genes:", length(non_matching_genes_universe), "\n")

filtered_NCBItoENSEMBL_universe <- NCBItoENSEMBL[
  NCBItoENSEMBL$Symbol %in% matching_genes_universe,
]
filtered_GaToDr_universe <- GaToDr[
  GaToDr$Gene.stable.ID %in% filtered_NCBItoENSEMBL_universe$Ensembl.GeneIDs,
]

filtered_GaToDr_universe_veryconservative <- filtered_GaToDr_universe[
  filtered_GaToDr_universe$Zebrafish.homology.type == "ortholog_one2one" &
    filtered_GaToDr_universe$Zebrafish.orthology.confidence..0.low..1.high. == 1,
]
universe_genes_veryconservative <- unique(
  filtered_GaToDr_universe_veryconservative$Zebrafish.gene.stable.ID
)

GO_results <- enrichGO(
  gene = genes_to_test,
  OrgDb = "org.Dr.eg.db",
  universe = universe_genes_veryconservative,
  keyType = "ENSEMBL",
  ont = "BP"
)

results <- as.data.frame(GO_results)
write.csv(results, "GO_results.csv", row.names = FALSE)

plot(barplot(GO_results, showCategory = 15))
