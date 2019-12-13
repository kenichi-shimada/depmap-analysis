# Analysis of DepMap

This repostory contains codes to find out 'classes' of essential modalities of genes from [the DepMap project data of the Broad Intitute](https://depmap.org/portal/download/).

## 01_pre_processing
Associated files are downloaded from DepMap project and preprocessed.
- [preprocessing_19q3.r](01_pre_processing/preprocessing_19q3.r): preprocess the data for further analysis.

## 02_compute_perturbation_score
Perturbation score is calculated from the CRISPR and shRNA efficacy.

- [score_calculation.r](02_compute_perturbation_score/score_calculation.r)
- [generate_fig2_figs1.r](02_compute_perturbation_score/generate_fig2_figs1.r)

## 03_find_essential_genes
Essential genes from the combined perturbation scores are identified.

- [find_ess_ts_genes.r](03_find_essential_genes/find_ess_ts_genes.r)

## 04_gsea
Pathways that are overrepresented among generally or selectively essential genes are identified using gene set enrichment analysis.

- [msigdb_gsea_prep.r](04_gsea/msigdb_gsea_prep.r) : objects prep for GSEA
- [comp_gsea_eff_sel.r](04_gsea/comp_gsea_eff_sel.r)
- [comp_gsea_eff_sel.sh](04_gsea/comp_gsea_eff_sel.sh)
- [plot_stepfun.r](04_gsea/plot_stepfun.r) 
- [comp_gsea_eff_sel_assemble_1e7.r](04_gsea/comp_gsea_eff_sel_assemble_1e7.r)

## 05_echodots
ECHODOTS, a tSNE and DBSCAN-based cluster analysis, is performed among essential genes. See [05_echodots](05_echodots) for the details.
