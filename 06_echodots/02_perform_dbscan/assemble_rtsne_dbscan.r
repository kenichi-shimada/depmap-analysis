library(dbscan)
library(clue)

j <- as.numeric(commandArgs(TRUE)[1]) ## i in 1:1026

ps <- expand.grid(d=30:200,k=1:6)

d <- ps$d[j] ## d in 30:200
k <- ps$k[j] ## k in 1:6
c1 <- seq(0,1,.2)[k]

fns <- paste0("rtsne_ori_",1:200,".rds")

##
setwd(rda_dir);setwd(paste0("r",c1));setwd("rtsne_multi_thres_1")
if(0)table(fns %in% dir()) ## made sure 200 files are all available	


rts <- lapply(fns,function(fn){
	rt <- readRDS(fn)
	return(rt)
})

dbscan.mem <- function(x=rt,d=d){
	knn.th <- diff(range(x))/d
	cl <- dbscan(x,eps=knn.th,minPts=2)
	mem <- cl$cl
	return(mem)
}

mems <- lapply(fns,function(fn){
	rt <- readRDS(fn)
	mem <- dbscan.mem(x=rt,d=d)
	return(mem)
})

##
set.seed(12345)
cl.pt <- lapply(mems,as.cl_hard_partition)
cons <- cl_consensus(cl.pt,method="SE",control=list(maxiter=1000)) ## soft/euclidean

## 
# ori.pred.idx.all <- apply(cons$.Data,1,which.max) 
# ori.pred.max <- apply(cons$.Data,1,max)
# tab.pred.idx.all <- table(ori.pred.idx.all)
# s <- sort(tab.pred.idx.all,decreasing=T)
# p.mems <- t(apply(cons$.Data[,as.numeric(names(s))],1,function(x)x/sum(x)))
# lik <- apply(p.mems,1,function(x)max(x))

# new.pred.idx.all <- apply(p.mems,1,which.max)
# names(new.pred.idx.all) <- rownames(rt)

# mem <- as.numeric(factor(new.pred.idx.all,levels=seq(max(new.pred.idx.all,na.rm=T))))
# names(mem) <- names(new.pred.idx.all)

setwd(rda_dir);setwd(paste0("r",c1));setwd("clue_multi_thres")

# save(mem,lik,file=paste0("d_",d,".rda"))
saveRDS(cons,file=paste0("cons_",d,".rds"))