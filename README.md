# Analysis of DepMap

This repostory contains codes to find out 'classes' of essential modalities of genes from [the DepMap project data of the Broad Intitute](https://depmap.org/portal/download/).

## 01_pre_processing
Download and preprocess associated files from DepMap project.

## 02_compute_perturbation_score
Impute and combine CRISPR and shRNA dependency scores. (Updated: now the dependency scores are defined using multple mixing ratios of CRISPR and shRNA)

## 03_find_essential_genes
Computed two parameters, 'efficacy' and 'selectivity', and the essential genes from the combined perturbation scores are identified for each mixing ratio.

## 04_adam
Lineage-dependent and common gene essentiality are computed using ADaM2 algorithm.

## 05_gsea
Pathways that are overrepresented among generally or selectively essential genes are identified using gene set enrichment analysis.

## 06_echodots
ECHODOTS, a tSNE and DBSCAN-based cluster analysis, is performed among essential genes. See [05_echodots](06_echodots) for the details.

## 07_prep_for_shinyDepMap
Files to be loaded in [shinyDepMap](https://labsyspharm.shinyapps.io/depmap) are computed.