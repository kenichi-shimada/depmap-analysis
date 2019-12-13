library(RColorBrewer)
library(gplots)

setwd(rda_dir)
load(file="gsea_data.rda") # df,msigdb.1

## in O2
setwd(rda_dir);setwd("fgsea")
fns <- paste0("fgsea_1e7_",1:2,".rds")
fns %in% dir()

slps  <- list()
for(i in 1:2){
	fn <- fns[[i]]
	res <- readRDS(fn)
	slps[[i]] <- sign(res$ES) * -log10(res$pval)
}

slps <- data.frame(do.call(cbind,slps))
names(slps) <- c("eff","sel")

## Fig. 3E
setwd(plot_dir);setwd("gsea")
pdf("smoothscatter_pathway-pvals.pdf",width=5,height=5)
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
setwd(plot_dir);setwd("gsea")

n1 <- sapply(msigdb.1[p1],length)
n2 <- sapply(msigdb.1[p2],length)
n3 <- sapply(msigdb.1[p3],length)

min.n <- min(c(n1,n2,n3))
rng <- diff(range(c(n1,n2,n3)))

sc1 <- sideCols[round((n1-min.n)/rng*8)+1]
sc2 <- sideCols[round((n2-min.n)/rng*8)+1]
sc3 <- sideCols[round((n3-min.n)/rng*8)+1]

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

# selective & essnetial
pdf("heatmap_overlap_1.pdf",width=5,height=5) 
heatmap.2(ois1,trace="none",margins=c(12,12),cexRow=.4,cexCol=.4,
	col=cols,RowSideColors=sc1,ColSideColors=sc1,main="essential & selective")
dev.off()

##
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

## commonly essential
pdf("heatmap_overlap_2.pdf",width=5,height=5) 
heatmap.2(ois2,trace="none",margins=c(12,12),cexRow=.5,cexCol=.5,
	col=cols,RowSideColors=sc2,ColSideColors=sc2,main="essential")
dev.off()

## selective; not so essential
pdf("heatmap_overlap_3.pdf",width=5,height=5)

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

## cell cycle genes in "selective & index" and "selective"
cc.1 <- grep("cell_",p1,value=T,ignore.case=T)[2]
cc.2 <- c(grep("cell_",p2,value=T,ignore.case=T),grep("m_g1",p2,value=T,ignore.case=T))

g1 <- unique(unlist(msigdb.1[cc.1]))
g2 <- unique(unlist(msigdb.1[cc.2]))
g12 <- intersect(g1,g2)

sym1 <- unlist(mget(g1,org.Hs.egSYMBOL))
sym2 <- unlist(mget(g2,org.Hs.egSYMBOL))
sym12 <- unlist(mget(g12,org.Hs.egSYMBOL))
sym2[!sym2 %in% sym1]
gn <- rownames(df);gn <- gn[!gn %in% c(g1,g2)]


if(0){
	## Fig. S4
	setwd(src_dir)
	setwd(plot_dir);setwd("gsea")
	pdf("cell_cycle_genes.pdf",width=5,height=5)
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(4,1), c(1,4), TRUE)

	cols <- c("green3","red")#brewer.pal(5,"Set1")[c(2,1)]
	par(mar = c(4,4,1,1))
	smoothScatter(df,pch=20,cex=NA,
		xlab="Efficacy (2.5 percentile)",
		ylab="Selectivity")
	points(df[gn,],pch=20,col=1,cex=.3)
	points(df[g1,],pch=20,col=cols[1],cex=.8)
	points(df[g2,],pch=20,col=cols[2],cex=.8)

	legend(-2.7,15,c("Cell cycle (ess)","Cell cycle (ess & sel)","the others"),
		col=c("red","green3","black"),pch=20,pt.cex=c(.8,.5,.3))

	par(mar = c(0,4,1,1))
	plot.ecdf(df[,1], col=1,verticals=T,pch=NA,main="",las=1,ylab="F(x)",xlim=range(df[[1]]),axes=F)
	axis(2,las=1)
	box()
	plot.ecdf(df[g1,1], col="green3",pch=NA,verticals=T,add=T)
	plot.ecdf(df[g2,1], col="red",pch=NA,verticals=T,add=T)

	par(mar = c(4,0,1,1))
	plot.stepfun(df[,2], col=1,verticals=T,pch=NA,horiz=T,main="",las=2,ylab="F(x)",xlim=range(df[[2]]),axes=F)
	axis(1,las=2)
	box()
	abline(v=0:1,col="grey80",lty=2)
	plot.stepfun(df[g1,2], col="green3",pch=NA,verticals=T,add=T,horiz=T,las=2)
	plot.stepfun(df[g2,2], col="red",pch=NA,verticals=T,add=T,horiz=T,las=2)

	dev.off()
}

