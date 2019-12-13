library(dbscan)
library(clue)

ind <- as.numeric(commandArgs(TRUE)[1]) ## i in 30:200
fns <- paste0("rtsne_2492gns_20000_original_",1:200,".rds")

if(0)table(fns %in% dir()) ## made sure 200 files are all available	

##
setwd(rda_dir);setwd("rtsne_original")

rts <- lapply(fns,function(fn){
	rt <- readRDS(fn)
	return(rt)
})

clue.mem <- function(x=rt,ind=ind){
	knn.th <- diff(range(x))/ind
	cl <- dbscan(x,eps=knn.th,minPts=2)
	mem <- cl$cl
	return(mem)
}

mems <- lapply(fns,function(fn){
	rt <- readRDS(fn)
	mem <- clue.mem(x=rt,ind=ind)
	return(mem)
})

max.n <- sapply(mems,max) ## 133 - 213, mode(185)

##
set.seed(12345)
cl.pt <- lapply(mems,as.cl_hard_partition)
system.time(cons <- cl_consensus(cl.pt,method="SE",control=list(maxiter=1000))) ## soft/euclidean

## 
ori.pred.idx.all <- apply(cons$.Data,1,which.max) 
ori.pred.max <- apply(cons$.Data,1,max)
tab.pred.idx.all <- table(ori.pred.idx.all)
s <- sort(tab.pred.idx.all,decreasing=T)
p.mems <- t(apply(cons$.Data[,as.numeric(names(s))],1,function(x)x/sum(x)))
lik <- apply(p.mems,1,function(x)max(x))

new.pred.idx.all <- apply(p.mems,1,which.max)
names(new.pred.idx.all) <- rownames(rt)

mem <- as.numeric(factor(new.pred.idx.all,levels=seq(max(new.pred.idx.all,na.rm=T))))
names(mem) <- names(new.pred.idx.all)

setwd(rda_dir);setwd("clue_original")
save(max.n,mem,lik,file=paste0("d_",ind,".rda"))