library(dbscan)
library(clue)
library(RColorBrewer)


## missing (genes missed out from dbscan)
setwd(rda_dir);setwd("rtsne")
fns.1 <- paste0("rtsne_",1:200,".rds")

clue.mem <- function(x=rt,ind=ind){
	knn.th <- diff(range(x))/ind
	cl <- dbscan(x,eps=knn.th,minPts=2)
	mem <- cl$cl
	return(mem)
}

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

## clusters representing missing values
setwd(rda_dir)
ess.genes <- readRDS(file="ess_2492_genes.rds")

setwd(rda_dir);setwd("clue")

mis.cl <- sapply(1:171,function(i){
	fn <- paste0("d_",(30:200)[i],".rda")
	load(fn)
	tmp <- unique(mem[m0[,i]>=100])
	if(length(tmp)==0)tmp <- 0	
	return(tmp)
})

##
fns <- paste0("d_",30:200,".rda")
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
})

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

ds <- as.character(30:200)
names(mems) <- names(max.ns) <- names(liks) <- names(tabs) <- names(ratios) <- ds

setwd(rda_dir)
save(mems,max.ns,liks,tabs,ratios,ds,file="ecodots.rda")	

