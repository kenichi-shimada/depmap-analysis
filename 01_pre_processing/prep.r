## we have to set path to three directories (can be set in .Rprofile file):
# src_dir : directory containing src (particularly for Slurm scripts)
# rda_dir : directory containing R objects (*.rda, *.rds)
# plot_dir : directory containing plots


## make sub-directories under plot_subdirs
plot_subdirs <- c("01_genewise_relationship",
	"02_compare_efficacy_measure",
	"03_histogram_essential_genes",
	"04_efficacy_selectivity_relationship",
	"05_gsea",
	"06_echodots_prep",
	"07_cluster_eff_sel")

setwd(plot_dir)
for(subdir in plot_subdirs){
	dir.create(subdir)
}


## make sub-directories under rda_subdirs
rda_subdirs <- c("fgsea","rtsne","clue")

setwd(rda_dir)
for(subdir in rda_subdirs){
	dir.create(subdir)
}

## make sub-directories under src_subdirs
src_subdirs <- c("fgsea_logs","rtsne_logs","clue_logs","impute_logs")

setwd(src_dir)
for(subdir in src_subdirs){
	dir.create(subdir)
}

