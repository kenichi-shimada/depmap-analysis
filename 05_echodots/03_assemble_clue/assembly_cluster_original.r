library(dbscan)
library(clue)
library(RColorBrewer)

setwd(rda_dir)
ess.genes <- readRDS(file="ess_2492_genes.rds")
mis.cl <- readRDS(file="cluster_missed.rds")

d <- "rtsne_original"
pre <- "rtsne_2492gns_20000_original_"

setwd(rda_dir);setwd("clue_original")
fns <- paste0("clue_original_",30:200,".rda")
# table(dir() %in% fns)

ds <- as.character(30:200)


mems <- lapply(seq(fns),function(idx){
	fn <- fns[idx]
	load(fn)
	mis <- mis.cl[idx]
	tab <- sort(table(mem),decreasing=T)
	uniq.mem <- seq(tab)
	uniq.mem <- uniq.mem[!uniq.mem %in% c(mis,names(tab[tab==1]))]
	new.lab <- seq(uniq.mem)
	names(new.lab) <- uniq.mem
	mem.1 <- as.numeric(new.lab[as.character(mem)])
	return(mem.1)
	# table(mem.1,useNA="ifany")
})
# sapply(mems,function(x)sum(is.na(x)))

max.ns <- lapply(fns,function(fn){
	load(fn)
	return(max.n)
})

liks <- lapply(1:171,function(ind){
	fn <- fns[ind]
	load(fn)
	mem <- mems[[ind]]
	has.na <- is.na(mem)
	lik[has.na] <- NA
	return(lik)
})

tabs <- sapply(mems,function(x){
	x1 <- table(x)
	missing <- length(x)-sum(x1)
	return(c(x1[1:2],sum(x1[-(1:2)]),missing))
})

ratios <- sapply(mems,function(x){
	x1 <- table(x)
	return(x1[1]/x1[2])
})

# colnames(tabs) <- ds
names(mems) <- names(max.ns) <- names(liks) <- names(tabs) <- names(ratios) <- ds

setwd(rda_dir)
save(mems,max.ns,liks,tabs,ratios,ds,file="ecodots_original.rda")	

## missing (genes missed out from dbscan)
setwd(rda_dir);setwd("rtsne_original")
m0 <- sapply(30:200,function(ind){
	cat("*")
	m <- sapply(fns.1,function(fn){
		rt <- readRDS(fn)
		mem <- clue.mem(x=rt,ind=ind)
		return(mem==0)
	})
	m <- rowSums(m)
	return(m)
})

##
mis.cl <- sapply(1:171,function(i){
	tmp <- unique(mems[[i]][m0[,i]>=100])
	if(length(tmp)==0)tmp <- 0	
	return(tmp)
})
mis.cl <- unlist(mis.cl)

setwd(rda_dir)
saveRDS(mis.cl,file="cluster_missed.rds")
