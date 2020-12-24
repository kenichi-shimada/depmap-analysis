library(dbscan)
library(clue)
library(RColorBrewer)

## missing (genes missed out from dbscan)
dbscan.mem <- function(x=rt,ind=ind){
	knn.th <- diff(range(x))/ind
	cl <- dbscan(x,eps=knn.th,minPts=2)
	mem <- cl$cl
	return(mem)
}

setwd(rda_dir)
x <- load(file="find_essential_genes.rda") 
#is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns

for(k in c(1:6)){
	cat(paste("k =",k,"\n"))
	## miss-classified
	c1 <- seq(0,1,.2)[k]

	# setwd(rda_dir);setwd(paste0("r",c1))
	# x <- load("ess.genes.rda")

	setwd(rda_dir);setwd(paste0("r",c1));
	setwd("clue_multi_thres")
	ncls <- sort(as.numeric(sub("cons_([0-9]+)\\.rds","\\1",dir(pattern="cons"))))

	setwd(rda_dir);setwd(paste0("r",c1));
	setwd("rtsne_multi_thres_1")
	fns <- paste0("rtsne_ori_",1:200,".rds")
	fns <- fns[fns %in% dir()]

	not.clustered <- sapply(ncls,function(ind){
		cat("*")
		m <- sapply(fns,function(fn){
			rt <- readRDS(fn)
			mem <- dbscan.mem(x=rt,ind=ind)
			return(mem==0)
		})
		m <- rowSums(m)
		n <- 200/2
		idx <- (m >= n)
		return(idx)
	})

	##
	if(0){
		setwd(rda_dir);setwd(paste0("r",c1))
		ess.genes.all <- readRDS(file="ess_genes_6thres.rds")
		if(0){
			ess.genes <- unique(unlist(ess.genes.all[1:5]))
		}else{
			ess.genes <- unique(unlist(ess.genes.all))
		}
	}
	ess.genes <- names(which(is.ess[[k]]$q1))

	# mis.cl <- readRDS(file="cluster_missed.rds")
	ds <- as.character(ncls)

	setwd(rda_dir);setwd(paste0("r",c1));
	setwd("clue_multi_thres")
	fns <- paste0("cons_",ds,".rds")
	table(fns %in% dir())

	mls <- parallel::mclapply(seq(fns),function(idx){
		cat("*")
		fn <- fns[idx]
		cons <- readRDS(fn)

		missed <- not.clustered[,idx]

		ori.pred.idx.all <- apply(cons$.Data,1,which.max) 
		ori.pred.idx.all[missed] <- NA

		tab.pred.idx.all <- table(ori.pred.idx.all)
		s <- sort(tab.pred.idx.all,decreasing=T)
		p.mems <- t(apply(cons$.Data[,as.numeric(names(s))],1,function(x)x/sum(x)))

		new.pred.idx.all <- apply(p.mems,1,which.max)
		lik <- apply(p.mems,1,max)

		names(lik) <- names(new.pred.idx.all) <- rownames(rt)
		lik[missed] <- new.pred.idx.all[missed] <- NA

		mem <- as.numeric(factor(new.pred.idx.all,levels=seq(max(new.pred.idx.all,na.rm=T))))
		names(mem) <- names(new.pred.idx.all)
		more.missed <- mem %in% as.numeric(names(which(table(mem)==1)))

		lik[more.missed] <- mem[more.missed] <- NA

		max.n <- max(mem,na.rm=T)

		# mis <- mis.cl[idx]
		tab <- sort(table(mem),decreasing=T)

		new.lab <- names(tab)

		mem.1 <- as.numeric(factor(mem,levels=new.lab))
		# return(all(diff(table(mem.1))<=0)) ## ALL TRUE
		return(list(mem=mem.1,lik=lik))
	},mc.cores=8)

	mems <- lapply(mls,function(x)x$mem)
	liks <- lapply(mls,function(x)x$lik)

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
	names(mems) <- names(liks) <- names(tabs) <- names(ratios) <- ds

	setwd(rda_dir);setwd(paste0("r",c1));
	save(mems,liks,tabs,ratios,ds,file="ecodots_ori_v3.rda")
}
