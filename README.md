# Project

The repeated evolution of reduced skeletal armor in threespine stickleback provides a powerful model for understanding the genetic basis of morphological change. We investigated the genetic and developmental mechanisms underlying the loss of the second dorsal spine in a freshwater stickleback population from North Uist, Scotland. Crosses between spineless freshwater and fully spined marine individuals confirmed a genetic basis for the trait, with inheritance patterns inconsistent with simple Mendelian expectations. A bulk segregant analysis of F3 hybrids revealed a strong genomic signal on chromosome VI, overlapping the hoxdb cluster. Developmental staining and bulk RNA sequencing showed delayed and incomplete cartilage formation and altered gene expression during critical stages of dorsal development in spineless fish, notably of the hoxdb genes. Our findings indicate that while the loss of dorsal spines has evolved repeatedly, it may involve distinct genetic mechanisms across populations, with chromosome VI playing a central role in this case.

This repository contains scripts for DNA sequencing, RNA-seq, annotation, and phenotypic analyses.

The code was split from one larger analysis file into smaller scripts so the workflow is easier to follow, run step by step, and share on GitHub.

## Project Structure

```text
project/
├── config/
│   ├── paths.yaml
│   └── params.yaml
├── scripts/
│   ├── 00_setup/
│   ├── 01_alignment/
│   ├── 02_bam_processing/
│   ├── 03_pileup/
│   ├── 04_freqsum/
│   ├── 05_af_differentiation/
│   ├── 06_rnaseq/
│   ├── 07_annotation/
│   └── 08_phenotypes/
└── README.md
```

## Before Running

Edit the paths in:

```text
config/paths.yaml
```

The scripts still contain the original paths from the analysis code. Some are placeholders such as:

```text
/my/folder/path
```

Replace those with the correct paths on your computer or cluster.

Analysis settings such as chromosomes, coverage thresholds, SLURM resources, and plotting parameters are collected in:

```text
config/params.yaml
```

## Data

You can find the project files in the European Nucleotide Archive under [Project: PRJEB91223](https://www.ebi.ac.uk/ena/browser/view/PRJEB91223). Large data files should usually not be uploaded to GitHub. This includes:

- FASTQ files
- SAM/BAM files
- pileup files
- large result tables
- genome indexes

Keep those files on your computer, server, or cluster storage, and use the paths in `config/paths.yaml` to point to them.

## Workflow

Run scripts from the main project folder.

### 00 Setup

Split large FASTQ files and trim RNA-seq adapters.

```bash
bash scripts/00_setup/split_fastq.sh
sbatch scripts/00_setup/trim_adapters.sh
```

### 01 Alignment

Run DNA alignment with Novoalign and RNA-seq alignment with STAR.

```bash
bash scripts/01_alignment/novoalign.sh
sbatch scripts/01_alignment/star_index.sh
sbatch scripts/01_alignment/star_align.sh
```

### 02 BAM Processing

Convert SAM files to BAM and merge BAM files.

```bash
bash scripts/02_bam_processing/sam_to_bam.sh
bash scripts/02_bam_processing/merge_bam.sh
```

### 03 Pileup

Generate pileup files from BAM files.

```bash
Rscript scripts/03_pileup/pileup.R
```

### 04 FreqSum

Convert pileup files to FreqSum format and split by chromosome.

```bash
Rscript scripts/04_freqsum/freqsum_parallel.R
Rscript scripts/04_freqsum/freqsum_fast.R
bash scripts/04_freqsum/split_by_chr.sh
```

### 05 Allele Frequency Differentiation

Calculate allele frequency differences between populations.

```bash
sbatch scripts/05_af_differentiation/run_afd_parallel.sh
Rscript scripts/05_af_differentiation/compute_afd.R VI 80 300 80 300
Rscript scripts/05_af_differentiation/afd_fast.R
```

### 06 RNA-seq

Run DESeq2, PCA, volcano plots, and GO analysis.

```bash
Rscript scripts/06_rnaseq/deseq_analysis.R
Rscript scripts/06_rnaseq/go_analysis.R
```

The files below are helper scripts used by the main RNA-seq analysis:

```text
scripts/06_rnaseq/pca_analysis.R
scripts/06_rnaseq/volcano_plots.R
```

### 07 Annotation

Extract genes from a selected chromosome in the GFF annotation file.

```bash
bash scripts/07_annotation/extract_genes.sh
```

### 08 Phenotypes

Run phenotypic analyses.

```bash
Rscript scripts/08_phenotypes/spine_location.R
Rscript scripts/08_phenotypes/abc_simulation.R
Rscript scripts/08_phenotypes/spine_length.R
Rscript scripts/08_phenotypes/development_sequence.R
```

## Software

The scripts use a mixture of command-line tools and R packages.

Command-line tools include:

- Novoalign
- STAR
- SAMtools
- Trimmomatic

R packages include:

- `data.table`
- `Rsamtools`
- `GenomicRanges`
- `DESeq2`
- `ggplot2`
- `ggrepel`
- `EnhancedVolcano`
- `pheatmap`
- `viridis`
- `reshape2`
- `clusterProfiler`
- `org.Dr.eg.db`
- `AnnotationDbi`
- `abc`
- `tidyverse`
- `patchwork`
- `Rfast`

## Notes

Some scripts are intended to run on a SLURM cluster with `sbatch`. Others can be run directly with `bash` or `Rscript`.

The current scripts are a cleaned first GitHub version. The next improvement would be to make every script read paths and parameters directly from `config/paths.yaml` and `config/params.yaml`.

## Citing this work

@article{herrera2025evolution,
  title={Evolution of threespine stickleback dorsal spines via hoxdb gene regulation},
  author={Herrera-Castillo, Carlos Manuel and Brechb{\"u}hl, Tanja and Fages, Antoine and Cameron MacColl, Andrew Donald and Dean, Laura L and Tschopp, Patrick and Berner, Daniel},
  journal={bioRxiv},
  pages={2025--08},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}
