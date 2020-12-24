## compute spearman correlation coefficients among essential genes
## run from 'compute_coef.slurm'

k <- as.numeric(commandArgs(TRUE))[1] # 1-6

c1 <- seq(0,1,.2)[k]

setwd(rda_dir)
x <- load(file="find_essential_genes.rda") #is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns

setwd(rda_dir);setwd(paste0("r",c1))
dep.scores <- readRDS("scores_15847g_423c.rds") 

gns <- rownames(dep.scores)

ess.genes <- gns[is.ess[[k]]$q1]

system.time(
	coef <- cor(t(dep.scores[ess.genes,]),method="spearman",use="pairwise.complete.obs")
)
saveRDS(coef,file="ess_coef.rds")

