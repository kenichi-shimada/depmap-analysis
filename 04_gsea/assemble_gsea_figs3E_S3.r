library(RColorBrewer)
library(gplots)

setwd(rda_dir)
load(file="gsea_data.rda") # df,msigdb.1

## in O2
setwd(rda_dir);setwd("fgsea")
fns <- paste0("fgsea_1e7_",1:2,".rds")
# fns %in% dir()

slps  <- list()
for(i in 1:2){
	fn <- fns[[i]]
	res <- readRDS(fn)
	slps[[i]] <- sign(res$ES) * -log10(res$pval)
}

slps <- data.frame(do.call(cbind,slps))
names(slps) <- c("eff","sel")

## Fig. 3E
setwd(plot_dir);setwd("05_gsea")
pdf("fig3E_smoothscatter_pathway-pvals.pdf",width=5,height=5)
par(mar=c(4,4,1,1))
smoothScatter(slps$eff,slps$sel,pch=NA,
	xlab="Signed log-P (efficacy)",
	ylab="Signed log-P (selectivity) ")
points(slps[[1]],slps[[2]],pch=20,cex=.3)
abline(h=0,v=0)

polygon(-7+c(1,1,-1,-1,1)*.1,7+c(-1,1,1,-1,-1)*.1,border=2,lwd=2)
polygon(-7+c(1,1,-1,-1,1)*.1,0+c(-1,1,1,-1,-1)*.25,border=2,lwd=2)
polygon(-0+c(1,1,-1,-1,1)*.25,7+c(-1,1,1,-1,-1)*.1,border=2,lwd=2)

text(c(-7,-7,0),c(7.2,-0.2,7.2),
	c("selective and essential","essential","selective"),
	col=2,pos=4)
dev.off()


##
p1 <- names(msigdb.1)[which(slps$eff < -6.99999 & slps$sel > 6.99999)]
p2 <- names(msigdb.1)[which(slps$eff < -6.99999 & abs(slps$sel) < .25)] # ef specific
p3 <- names(msigdb.1)[which(abs(slps$eff) < .25 & slps$sel > 6.99999)] # sel specific

sideCols <- brewer.pal(9,"YlGnBu")
cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)

## Fig. S3B
setwd(plot_dir);setwd("05_gsea")

n1 <- sapply(msigdb.1[p1],length)
n2 <- sapply(msigdb.1[p2],length)
n3 <- sapply(msigdb.1[p3],length)

min.n <- min(c(n1,n2,n3))
rng <- diff(range(c(n1,n2,n3)))

sc1 <- sideCols[round((n1-min.n)/rng*8)+1]
sc2 <- sideCols[round((n2-min.n)/rng*8)+1]
sc3 <- sideCols[round((n3-min.n)/rng*8)+1]

## Fig. S3
## selective & essnetial
ois1 <- sapply(msigdb.1[p1],function(x){
	sapply(msigdb.1[p1],function(y){
		n1 <- length(x)
		n2 <- length(y)
		n3 <- length(intersect(x,y))
		oi <- n3/min(n1,n2)
		return(oi)
	})
})
rownames(ois1) <- colnames(ois1) <- sub("^(kegg|go|reactome)_","",tolower(sub("^[^\\.]+\\.","",rownames(ois1))))

pdf("figS3_heatmap_overlap_1.pdf",width=5,height=5) 
heatmap.2(ois1,trace="none",margins=c(12,12),cexRow=.4,cexCol=.4,
	col=cols,RowSideColors=sc1,ColSideColors=sc1,main="essential & selective")
dev.off()

## commonly essential
ois2 <- sapply(msigdb.1[p2],function(x){
	sapply(msigdb.1[p2],function(y){
		n1 <- length(x)
		n2 <- length(y)
		n3 <- length(intersect(x,y))
		oi <- n3/min(n1,n2)
		return(oi)
	})
})

rownames(ois2) <- colnames(ois2) <- sub("^(kegg|go|reactome)_","",tolower(sub("^[^\\.]+\\.","",rownames(ois2))))

pdf("figS3_heatmap_overlap_2.pdf",width=5,height=5) 
heatmap.2(ois2,trace="none",margins=c(12,12),cexRow=.5,cexCol=.5,
	col=cols,RowSideColors=sc2,ColSideColors=sc2,main="essential")
dev.off()

## selective; not so essential
ois3 <- sapply(msigdb.1[p3],function(x){
	sapply(msigdb.1[p3],function(y){
		n1 <- length(x)
		n2 <- length(y)
		n3 <- length(intersect(x,y))
		oi <- n3/min(n1,n2)
		return(oi)
	})
})
rownames(ois3) <- colnames(ois3) <- sub("^(kegg|go|reactome)_","",tolower(sub("^[^\\.]+\\.","",rownames(ois3))))

pdf("figS3_heatmap_overlap_3.pdf",width=5,height=5)
heatmap.2(ois3,trace="none",margins=c(12,12),cexRow=.5,cexCol=.5,
	col=cols,RowSideColors=sc3,ColSideColors=sc3,main="selective")
dev.off()

## color index
pdf("index.pdf",width=3,height=3)
par(mar=c(5,5,4,1))
image(array(1:9,c(9,1)),col=sideCols,axes=F)
box()
xs <- par()$usr[1:2]
n.rng <- range(n1,n2,n3)
axis(1,
	at=seq(xs[1],xs[2],length=10),las=2,
	labels=round(seq(n.rng[1],n.rng[2],length=10)))
dev.off()

