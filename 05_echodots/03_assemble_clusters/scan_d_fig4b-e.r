library(dbscan)
library(clue)
library(RColorBrewer)

setwd(rda_dir)
load(file="ecodots.rda")	 # mems,liks,ds,ratios,optims,types,


##
le2 <- range(as.numeric(names(which(ratios <= 2))))
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

setwd(rda_dir)
three.ds <- ds[sig.i]
saveRDS(three.ds,file="three_ds.rds")

##
setwd(plot_dir);setwd("06_echodots_prep")

## Fig. 4B
pdf("fig4B_barplot-tab.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
tabs.1 <- tabs.0 <- tabs
optim <- sig.i
tabs.0[,optim] <- 0
tabs.1[,-optim]  <- 0
tab <- rbind(tabs.0,tabs.1)
tab <- tab[,ds %in% as.character(30:200)]

ds.1 <- ds
ds.1[as.numeric(ds.1) %% 20 != 0] <- NA
colnames(tab) <- ds.1

cols <- c(grey.colors(4),rev(brewer.pal(4,"Reds")))
x <- barplot(tab,beside=F,las=2,col=cols,border=NA,main="cluster size",
	xlab="d",ylab="# essential genes",ylim=c(0,2500),cex.main=2,space=0)
dev.off()

## Fig. 4C
pdf("fig4C_barplot-ratio.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
optim <- sig.i
idx <- min(which(ratios <2))
cols <- rep("black",length(ds))

plot(as.numeric(ds),ratios,col=cols,ylab="ratio",type="n",
	main="ratio of 1st/2nd clusters",cex.main=2,xlab="d")
lines(as.numeric(ds),ratios)
abline(v=as.numeric(ds[sig.i]),lty=2)
abline(h=max(ratios[idx:171]),lty=2,col="grey70")
points(as.numeric(ds[sig.i]),ratios[sig.i],cex=1.5,pch=20)
text(ds[sig.i],10,paste0("d = ",ds[sig.i]),srt=90,pos=2)
dev.off()

## Fig. 4D
pdf("fig4D_plot-nclust_no-lik-filter.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
x <- barplot(nclust.1,beside=F,
	las=2,border=NA,main="# clusters",cex.main=2,
	xlab="d",ylab="# clusters",space=0)
abline(v=x[sig.i],lty=2)
text(x[sig.i],700,paste0("d = ",ds[sig.i]),srt=90,pos=2)
dev.off()
 
## Fig. 4E
pdf("fig4E_med_cl_size.pdf",width=5,height=5)
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


