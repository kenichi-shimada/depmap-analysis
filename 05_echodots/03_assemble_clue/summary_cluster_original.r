library(dbscan)
library(clue)
library(RColorBrewer)

setwd(rda_dir)
load(file="ecodots_original.rda")	 # mems,liks,ds,ratios,optims,types,

##
setwd(plot_dir);setwd("ecodots_prep")

##
ds <- 30:200
r <- ratios
le2 <- range(as.numeric(names(which(r <= 2))))
ind.le2 <- which(ds %in% as.character(le2))

mems.1 <- lapply(ds,function(d){
	d <- as.character(d)
	return(mems[[d]])
})
names(mems.1) <- ds

nclust.1 <-sapply(mems.1,function(x){
	x1 <- table(x)
	sum(x1>1)
})
names(nclust.1) <- 30:200
nc <- names(nclust.1)
nc[!(30:200) %in% as.character(seq(30,200,10))] <- NA
names(nclust.1) <- nc

sig.i <- sapply(c(200,400,600),function(x){
	which.min(abs(nclust.1-x))
})

## Fig. 4B
pdf("barplot-tab.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
tabs.1 <- tabs.0 <- tabs
optim <- sig.i
tabs.0[,optim] <- 0
tabs.1[,-optim]  <- 0
tab <- rbind(tabs.0,tabs.1)
tab <- tab[,ds %in% as.character(30:200)]

ds[as.numeric(ds) %% 20 != 0] <- NA
colnames(tab) <- ds

cols <- c(grey.colors(4),rev(brewer.pal(4,"Reds")))
x <- barplot(tab,beside=F,las=2,col=cols,border=NA,main="cluster size",
	xlab="d",ylab="# essential genes",ylim=c(0,2500),cex.main=2,space=0)
dev.off()

## Fig. 4C
pdf("barplot-ratio.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
ratio <- ratios
optim <- sig.i
idx <- min(which(ratio<2))
cols <- rep("black",length(ds))

plot(as.numeric(ds),ratio,col=cols,ylab="ratio",type="n",
	main="ratio of 1st/2nd clusters",cex.main=2,xlab="d")
lines(as.numeric(ds),ratio)
abline(v=as.numeric(ds[sig.i]),lty=2)
abline(h=max(ratio[idx:171]),lty=2,col="grey70")
points(as.numeric(ds[sig.i]),ratio[sig.i],cex=1.5,pch=20)
text(ds[sig.i],10,paste0("d = ",ds[sig.i]),srt=90,pos=2)
dev.off()

## Fig. 4D
pdf("plot-nclust_no-lik-filter.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
x <- barplot(nclust.1,beside=F,
	las=2,border=NA,main="# clusters",cex.main=2,
	xlab="d",ylab="# clusters",space=0)
abline(v=x[sig.i],lty=2)
text(x[sig.i],700,paste0("d = ",ds[sig.i]),srt=90,pos=2)
dev.off()
 
## Fig. 4E
pdf("med_cl_size.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
med.cl.size <- sapply(mems.1,function(x){
	median(table(x))
})
ds <- 30:200
plot(ds,med.cl.size,type="l",
	main="median cluster size",cex.main=2,
	xlab="d",ylab="# gene/cluster")
abline(v=ds[sig.i],lty=2)
dev.off()


