library(dbscan)
library(clue)
library(RColorBrewer)
library(gplots)

##
setwd(rda_dir)
x <- load(file="find_essential_genes.rda") 
ess.genes <- list()
for(k in c(1:6)){
	c1 <- seq(0,1,.2)[k]
	ess.genes[[k]] <- names(which(is.ess[[k]]$q1))
}
ess.gns <- ess.genes
n.ess.gns <- sapply(ess.gns,length)

n1 <- list()
df <- list()
for(k in c(1:6)){
	c1 <- seq(0,1,.2)[k]

	setwd(rda_dir);setwd(paste0("r",c1));
	load(file="ecodots_ori_v3.rda")	 # mems,liks,ds,ratios,optims,types,

	attr(tabs, "names")<- NULL
	## n1, n2, rest, missing

	idx <- min(which(ratios < 2)) # idx : ratio < 2

	## n1: # of genes in largest cluster and its index
	n.clust <-sapply(mems,function(x){
		x1 <- table(x)
		sum(x1>1)
	})
	names(n.clust) <- ds

	sig.i <- sapply(c(300,500,700),function(x){
		which.min(abs(n.clust-x))
	}) # sig.i : 300, 500, 700

	n1[[k]]	<- tabs[1,c(idx,sig.i)]
	names(n1[[k]]) <- ds[c(idx,sig.i)]

	##
	r.n1.all <- tabs[1,]/n.ess.gns[k]
	rel.n.clust <- n.clust/n.ess.gns[k] ## 

	df[[k]] <- data.frame(r.n1.n2=ratios,r.n1.all=r.n1.all,r.nclust.ngenes=rel.n.clust)
	## r.n1.n2: ratio n1/n2
	## r.n1.all: n1/all genes
	## r.nclust.ngenes: ratio # clusters / # points
}

names(n1) <- names(df) <- paste0("r",seq(0,1,.2))

ds <- 30:200

##
interval <- 0.0625
n.thres <- 0.15 # for r.nclust.ngenes

thres <- sapply(df,function(x){
	idx <- which(x$r.nclust.ngenes > n.thres)
	max.i <- idx[which.max(x$r.n1.n2[idx])]
	c(x$r.nclust.ngenes[max.i],x$r.n1.n2[max.i])
})
rownames(thres) <- c("r.nclust.ngenes","r.n1.n2")
max.r.n1.n2 <- max(thres[2,])

##
if(0){
	sig.idx <- lapply(df,function(x){ ## index/number of genes
		i1 <- min(which(x$r.n1.n2 < max.r.n1.n2))
		i2 <- min(which(x$r.n1.all < n.thres))
		# idx <- i1
		idx <- max(i1,i2)

		rind <- x$r.nclust.ngenes[idx] + interval * (0:2)
		rs <- sapply(rind,function(r1){
			which.min(abs(x$r.nclust.ngenes-r1))
		})
		return(rs)
	})

	setwd(rda_dir)
	saveRDS(sig.idx,file="sig.idx_6thres_v3.rds")
}else{
	setwd(rda_dir)
	sig.idx <- readRDS(file="sig.idx_6thres_v3.rds")
}


n.clust <- sapply(1:6,function(k){
	c1 <- seq(0,1,.2)[k]
	setwd(rda_dir);setwd(paste0("r",c1));
	x <- load(file="ecodots_ori_v3.rda")	 # mems,liks,ds,ratios,optims,types,

	nclust.1 <-sapply(mems,function(x){
		x1 <- table(x)
		sum(x1>1)
	})
	names(nclust.1) <- ds
	return(nclust.1)
})
colnames(n.clust) <- 1:6

ratio <- rev(c("100:0 (CRISPR)","80:20","60:40","40:60","20:80","0:100 (shRNA)"))
n.clust.3 <- t(sapply(1:6,function(i){
	n.clust[sig.idx[[i]],i]
}))
rownames(n.clust.3) <- ratio
colnames(n.clust.3) <- c("L","M","S")

setwd(plot_dir);setwd("07_echodots_prep")
write.csv(n.clust.3,file="nclust_3.csv")

##
x.rng <- range(sapply(df,function(x)x$r.nclust.ngenes))
y.rng <- range(sapply(df,function(x)x$r.n1.all))
plot(df[[k]]$r.nclust.ngenes,df[[k]]$r.n1.all,type="n",
	xlab="r.nclust.ngenes",ylab="# genes in the largest clust",
	main="# of genes in the largest cluster",
	xlim=x.rng,ylim=y.rng)
for(k in 1:6){
	lines(df[[k]]$r.nclust.ngenes,df[[k]]$r.n1.all,type="l",col=k)
	abline(v=df[[k]]$r.nclust.ngenes[sig.idx[[k]]],las=2,col=k)
	points(df[[k]]$r.nclust.ngenes[sig.idx[[k]]],
		df[[k]]$r.n1.all[sig.idx[[k]]],pch=20,col=k)
}
legend("topright",paste0("k=",1:6),col=1:6,lty=1)
# abline(h=2.2,lty=2)
# abline(h=0.1,lty=2)
abline(h=n.thres,lty=2)

##
setwd(plot_dir);setwd("07_echodots_prep")
par(mfrow=c(1,1))
k=1
plot(df[[k]]$r.nclust.ngenes,df[[k]]$r.n1.n2,type="n",ylim=c(0,10),
	xlab="r.nclust.ngenes",ylab="ratio (1st/2nd clusters)",
	main="ratio (N1/N2)")
for(k in 1:6){
	lines(df[[k]]$r.nclust.ngenes,df[[k]]$r.n1.n2,type="l",ylim=c(0,10),col=k)
	abline(v=df[[k]]$r.nclust.ngenes[sig.idx[[k]]],las=2,col=k)
	points(df[[k]]$r.nclust.ngenes[sig.idx[[k]]],
		df[[k]]$r.n1.n2[sig.idx[[k]]],pch=20,col=k)
	# points(thres[1,k],thres[2,k],pch=20,col=k)
}
legend("topright",paste0("k=",1:6),col=1:6,lty=1)
abline(h=max.r.n1.n2,lty=2)


## Fig. 4B
setwd(plot_dir);setwd("07_echodots_prep")

pdf("echodots_barplot-ratio_v3.pdf",width=9,height=6)
par(mfrow=c(2,3))

for(k in 1:6){
	ratios <- df[[k]]$r.n1.n2
	sig.i <- sig.idx[[k]]

	## Fig. 4C
	par(mar=c(4,4,2,2))

	cols <- rep("black",length(ds))

	plot(as.numeric(ds),ratios,col=cols,ylab="ratio",type="n",
		main=paste0("k=",k),cex.main=2,xlab="d",ylim=c(0,10))
		# main="ratio of 1st/2nd clusters",cex.main=2,xlab="d",ylim=c(0,10))
	lines(as.numeric(ds),ratios)
	abline(v=as.numeric(ds[sig.i]),lty=2)
	points(as.numeric(ds[sig.i]),ratios[sig.i],cex=1,pch=20,col=2)
	y <- par()$usr[3:4] %*% c(1,3)/4
	text(ds[sig.i],y,paste0("d = ",ds[sig.i]," (",round(ratios[sig.i],2),")"),srt=90,pos=2)
	# abline(h=c(1,1.5,2),lty=2,col="grey70")
}
dev.off()

##
setwd(plot_dir);setwd("07_echodots_prep")
pdf("echodots_barplot-ratio_cmbd_v3.pdf",width=5,height=5)
par(mar=c(5,5,3,3))

cols <- rev(brewer.pal(6,"Spectral"))
plot(as.numeric(ds),ratios,col=cols,ylab="ratio (N1/N2)",type="n",
	main="",cex.main=2,xlab="d",ylim=c(0,10))
	# main="ratio of 1st/2nd clusters",cex.main=2,xlab="d",ylim=c(0,10))

sapply(1:6,function(k){
	ratios <- df[[k]]$r.n1.n2
	lines(as.numeric(ds),ratios,col=cols[k])
})
sapply(1:6,function(k){
	ratios <- df[[k]]$r.n1.n2
	sig.i <- sig.idx[[k]]
	points(as.numeric(ds[sig.i]),ratios[sig.i],cex=1.5,pch=20,col=cols[k])
})

abline(h=max(thres[2,]),lty=2)
ratio <- rev(c("100:0 (CRISPR)","80:20","60:40","40:60","20:80","0:100 (shRNA)"))
legend("topright",ratio,title="CRISPR:shRNA",col=cols,lty=1,pch=20,pt.cex=1.5)

dev.off()

## Fig. 4C
setwd(plot_dir);setwd("07_echodots_prep")
pdf("echodots_barplot-tab_v3.pdf",width=9,height=6)
par(mfrow=c(2,3))
for(k in c(1:6)){
	c1 <- seq(0,1,.2)[k]
	setwd(rda_dir);setwd(paste0("r",c1));
	x <- load(file="ecodots_ori_v3.rda")	 # mems,liks,ds,ratios,optims,types,

	sig.i <- sig.idx[[k]]

	par(mar=c(4,4,2,2))

	tabs.1 <- tabs.0 <- as.matrix(tabs)
	optim <- c(sig.i)
	tabs.0[,optim] <- 0
	tabs.1[,-optim]  <- 0
	tab <- rbind(tabs.0,tabs.1)
	tab <- tab[,ds %in% as.character(30:200)]

	ds1 <- ds
	ds1[as.numeric(ds) %% 20 != 0] <- NA
	colnames(tab) <- ds1

	cols <- c(grey.colors(4),rev(brewer.pal(4,"Reds")))
	x <- barplot(tab,beside=F,las=2,col=cols,border=NA,
		xlab="d",ylab="# essential genes",ylim=c(0,length(ess.gns[[k]])),cex.main=2,space=0,
		main=paste0("k=",k))
	# abline(v=x[idx],lty=2)
	# abline(h=diag(ns)[k]*.1,col=4)
}
dev.off()


## Fig. 4D
setwd(plot_dir);setwd("07_echodots_prep")
pdf("echodots_plot-nclust_no-lik-filter_v3.pdf",width=9,height=6)
par(mfrow=c(2,3))
par(mar=c(4,4,2,2))
for(k in c(1:6)){
	c1 <- seq(0,1,.2)[k]
	setwd(rda_dir);setwd(paste0("r",c1));
	x <- load(file="ecodots_ori_v3.rda")	 # mems,liks,ds,ratios,optims,types,

	sig.i <- sig.idx[[k]]

	plot(ds,n.clust[,k],type="l",
		main=paste0("k=",k),cex.main=2,
		xlab="d",ylab="# clusters",lwd=2,col="steelblue")
	abline(v=ds[sig.i],lty=2)
	y <- par()$usr[4]#par()$usr[3:4] %*% c(1,9)/10
	text(ds[sig.i],y,paste0("d = ",ds[sig.i]," (",round(n.clust[sig.i,k],1),") "),
		srt=90,pos=2,cex=1)
}
dev.off()

##
pdf("echodots_plot-nclust_no-lik-filter_cmbd_v3.pdf",width=9,height=6)
par(mar=c(5,5,3,3))

cols <- brewer.pal(6,"Spectral")
ds <- 30:200
plot(range(ds),range(n.clust),type="n",
	main="",cex.main=2,
	xlab="d",ylab="# clusters",lwd=2,col="steelblue")
sapply(1:6,function(k){
	lines(ds,n.clust[,k],xlab="d",ylab="# clusters",lwd=2,col=cols[k])
})
sapply(1:6,function(k){
	sig.i <- sig.idx[[k]]
	points(ds[sig.i],n.clust[sig.i,k],cex=1.5,col=cols[k],pch=20)
})
dev.off()

## Fig. 4E
setwd(plot_dir);setwd("07_echodots_prep")
pdf("echodots_med_cl_size_v3.pdf",width=9,height=6)
par(mfrow=c(2,3))
par(mar=c(4,4,2,2))
for(k in c(1:6)){
	c1 <- seq(0,1,.2)[k]
	setwd(rda_dir);setwd(paste0("r",c1));
	x <- load(file="ecodots_ori_v3.rda")	 # mems,liks,ds,ratios,optims,types,
	med.cl.size <- sapply(mems,function(x){
		mean(table(x))
	})

	sig.i <- sig.idx[[k]]
	plot(ds,med.cl.size,type="l",
		main=paste0("k=",k),cex.main=2,
		xlab="d",ylab="mean cluster size",lwd=2,col="steelblue")
	abline(v=ds[sig.i],lty=2)
	y <- par()$usr[3:4] %*% c(1,3)/4
	text(ds[sig.i],y,paste0("d = ",ds[sig.i]," (",round(med.cl.size[sig.i],1),")"),srt=90,pos=2,cex=1)
}
dev.off()

