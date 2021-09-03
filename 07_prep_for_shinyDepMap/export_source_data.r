setwd(rda_dir)
x <- load("depmap_initial_19q3_v3.rda")

setwd(rda_dir)
graphs <- readRDS(file="depmap_initial_19q3_v3_graphs.rda")
coefs <- readRDS(file="depmap_initial_19q3_v3_coefs.rda")
p.scores <- readRDS(file="depmap_initial_19q3_v3_p.scores.rda")

setwd(rda_dir);setwd("source_data")

## source data 0: cell line
setwd(rda_dir);setwd("source_data")
write.csv(cl.info[c(1:3,6:9)],"cell_lines.csv")

## source data 1: dependency scores
for(i in 1:6){
	write.csv(p.scores[[i]],file=paste0("sd1_dep_score_",i,".csv"))
}

library(org.Hs.eg.db)
syms <- unlist(mget(rownames(p.scores[[1]]),org.Hs.egSYMBOL))
write.csv(syms,"symbols.csv")
head(syms)
sapply(syms,length))

## source data 2: efficacy and selectivity
library(dplyr)

df.combined <-lapply(dfs,function(x){ # c("q1","q2.5","q5","q10","q25")
	df <- do.call(cbind,lapply(x[2:6],function(y)y[2:4]))
})

setwd(rda_dir);setwd("source_data")
for(i in 1:6){
	write.csv(df.combined[[i]],file=paste0("sd2_eff_sel_",i,".csv"))
}

## source data 3: overrepresented pathways (done)

## source data 4: lineage
setwd(rda_dir)
is.lin.dep <- readRDS(file="is.lin.dep_v3.rds")
is.comm <- readRDS(file="is.comm.rds")

lin_deps <- lapply(1:6,function(i){
	rs <- rowSums(is.lin.dep[[i]][-1])
	ic <- is.comm[,i]
	x <- is.lin.dep[[i]] %>% 
		mutate(n_dependent_lineage= rs) %>%
		mutate(is_commonly_essential = as.logical(ic))
	return(x)		
})

setwd(rda_dir);setwd("source_data")
for(i in 1:6){
	write.csv(lin_deps[[i]],file=paste0("sd4_lineage_dependency_",i,".csv"))
}

## source data 5: membership/probability
setwd(rda_dir);setwd("source_data")
for(i in 1:6){
	write.csv(mems[[i]],file=paste0("sd5_clusters_",i,".csv"))
}
