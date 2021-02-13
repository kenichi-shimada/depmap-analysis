# Analysis of DepMap for generating shinyDepMap

This repostory contains codes to find out 'classes' of essential modalities of genes from [the DepMap project data of the Broad Intitute](https://depmap.org/portal/download/).
(Note: this repo is updated in Dec 2020. Now the dependency scores are calculated with various mixing ratios between CRISPR and shRNA. For the details, see an upcoming publicaiton.)

## 01 pre processing
Download and preprocess associated files from DepMap project.

## 02 compute perturbation score
Impute and combine CRISPR and shRNA dependency scores. (Updated: now the dependency scores are defined using multple mixing ratios of CRISPR and shRNA)

## 03 find essential genes
Computed two parameters, 'efficacy' and 'selectivity', and the essential genes from the combined perturbation scores are identified for each mixing ratio.

## 04 adam
Lineage-dependent and common gene essentiality are computed using ADaM2 algorithm.

## 05 gsea
Pathways that are overrepresented among generally or selectively essential genes are identified using gene set enrichment analysis.

## 06 echodots
ECHODOTS, a tSNE and DBSCAN-based cluster analysis, is performed among essential genes. See [06 echodots](06_echodots) for the details.

## 07 prep for shinyDepMap
Files to be loaded in [shinyDepMap](https://labsyspharm.shinyapps.io/depmap) are computed.
