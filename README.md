# Analysis of DepMap

This repostory contains codes to find out 'classes' of essential modalities of genes from [the DepMap project data of the Broad Intitute](https://depmap.org/portal/download/).

## 01_pre_processing
Associated files are downloaded from DepMap project and preprocessed.
- [prep.r](prep.r) : set up subdirectories under `data_dir`, `rda_dir`, `plot_dir`, `src_dir`.
- [preprocessing_19q3.r](01_pre_processing/preprocessing_19q3.r) : preprocess the data for further analysis.

## 02_compute_perturbation_score
Perturbation score is calculated from the CRISPR and shRNA efficacy.

- [score_calculation.r](02_compute_perturbation_score/score_calculation.r) : calculate perturbation scores.
- [fig2_figs1.r](02_compute_perturbation_score/fig2_figs1.r) : generates figs 2 and S1.

## 03_find_essential_genes
Essential genes from the combined perturbation scores are identified.

- [find_essential_genes_fig3A-D_S2.r](find_essential_genes_fig3A-D_S2.r) : compute essential genes.

## 04_gsea
Pathways that are overrepresented among generally or selectively essential genes are identified using gene set enrichment analysis.

- [gsea_prep.r](04_gsea/gsea_prep.r) : objects prep for GSEA.
- [gsea.r](04_gsea/gsea.r), [master_gsea.sh](04_gsea/master_gsea.sh) : Slurm batch script to compute GSEA.
- [assemble_gsea_figs3E_S3.r](04_gsea/assemble_gsea_figs3E_S3.r) : assemble results of GSEA but generate figs 3E and S3.

## 05_echodots
ECHODOTS, a tSNE and DBSCAN-based cluster analysis, is performed among essential genes. See [05_echodots](05_echodots) for the details.

## 06_prep_for_shinyDepMap
This is to prepare data to use for shinyDepMap. 

- [compile_loading_data_19q3.r](06_prep_for_shinyDepMap/compile_loading_data_19q3.r) : compile loading data for shinyDepMap.
- [upload_objects_on_s3_19q3.r](06_prep_for_shinyDepMap/upload_objects_on_s3_19q3.r) : upload objects to AWS S3.
- [graphs.r](06_prep_for_shinyDepMap/graphs.r) : generate figure 4H-J.

## Slurm batch job run on HPC
Following src files are uploaded to `src_dir` on HPC for Slurm batch job run. `sh` files are run by `$ sbatch sh_file_name`.

- [04_gsea/gsea.r](04_gsea/gsea.r), [04_gsea/master_gsea.sh](04_gsea/master_gsea.sh)
- [05_echodots/01_perform_tsne/rtsne.r](05_echodots/01_perform_tsne/rtsne.r), [05_echodots/01_perform_tsne/master_rtsne.sh](05_echodots/01_perform_tsne/master_rtsne.sh)
- [05_echodots/02_perform_dbscan/dbscan-clue.r](05_echodots/02_perform_dbscan/dbscan-clue.r), [05_echodots/02_perform_dbscan/master_dbscan-clue.sh](05_echodots/02_perform_dbscan/master_dbscan-clue.sh)


