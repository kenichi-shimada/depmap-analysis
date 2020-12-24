library(ggplot2)
library(KernSmooth) ## bkde2D
library(RColorBrewer)
library(org.Hs.eg.db)
library(sp)
library(dplyr)


## load data
setwd(rda_dir)
load(file="eff_dep_19q3.rda") ## eff.1, dep.1, both 15847 x 423 cells
ed <- data.frame(e=as.vector(eff.1),d=as.vector(dep.1))
ed1 <- as.data.frame(readRDS("non-missing-ed.rds"))

##
ed.lims <- c(-3,2)
e <- ed1$e
d <- ed1$d
idx <- which(e > ed.lims[1] & e < ed.lims[2] & d > ed.lims[1] & d < ed.lims[2])
e <- e[idx]
d <- d[idx]

xhist <- hist(e,breaks=500,col="pink",xlim=c(ed.lims[1],ed.lims[2]),plot=F)
yhist <- hist(d,breaks=500,col="pink",xlim=c(ed.lims[1],ed.lims[2]),plot=F)
top <- max(xhist$counts,yhist$counts)

## CRISPR and shRNA thresholds
setwd(src_dir);source("functions.r")
e.th <- find.threshold(e,
	xlab="dependency score",plot=TRUE,
	main="",
	side.fit="upper",side.thres="lower")

d.th <- find.threshold(d,
	xlab="dependency score",plot=TRUE,
	main="",
	side.fit="upper",side.thres="lower")

## 
all.a <- sum(ed$e < e.th & ed$d > d.th,na.rm=T)
all.b <- sum(ed$d < d.th & ed$e > e.th,na.rm=T)

gns <- rownames(eff.1)
pvals <- sapply(gns,function(g){
	in.a <- table(eff.1[g,] < e.th & dep.1[g,] > d.th)[c("TRUE","FALSE")]
	in.b <- table(eff.1[g,] > e.th & dep.1[g,] < d.th)[c("TRUE","FALSE")]
	in.a[is.na(in.a)] <- 0
	in.b[is.na(in.b)] <- 0

	not.a <- all.a - in.a
	not.b <- all.b - in.b

	mat.a <- rbind(in.a,not.a)
	mat.b <- rbind(in.b,not.b)
	p.a <- fisher.test(mat.a,alternative="greater")$p.value
	p.b <- fisher.test(mat.b,alternative="greater")$p.value	
	return(c(p.a,p.b))
})
nlps <- round(-log10(t(pvals)),3) ## four pvals are 0

table(crispr=nlps[,1]!=0,shrna=nlps[,2]!=0)[2:1,2:1]
head(sort(nlps[,2],decreasing=T),20)



## Fig. 1B
setwd(plot_dir);setwd("01_genewise_relationship")
pdf("efficacy_selectivity.pdf",width=5,height=5)

nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

## 1
par(mar = c(4,4,1,1))
gene.ed.plot(ed.lims=c(-3,2))

uniq.cols <- brewer.pal(6,"Spectral")[c(1,6)]
text(-2,1,"A\n(essential\nby CRISPR)",col=uniq.cols[1],cex=1)
text(1,-2,"B\n(essential\nby shRNA)",col=uniq.cols[2],cex=1)

abline(v=e.th,h=d.th,col="grey50")

## 2
par(mar = c(0,4,1,1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,col=c("grey","grey50")[(xhist$breaks < e.th)+1],border=NA,lwd=.25)


## 3
par(mar = c(4,0,1,1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0,col=c("grey","grey50")[(yhist$breaks < d.th)+1],border=NA,lwd=.25,horiz=TRUE)
dev.off()


## Fig. 1C
setwd(plot_dir);setwd("01_genewise_relationship")
pdf("gene_enrichment.pdf",width=5,height=5)
is.non0 <- nlps[,1]!=0 | nlps[,2]!=0

cols <- c(1,uniq.cols[c(2,1)])[(nlps[,1]>3)*2+(nlps[,2]>3)*1+1] # "#532E8C","#3683B5"

par(mar=c(5,5,3,3))
plot(nlps[is.non0,],pch=20,
xlab="P-value (Enriched in points in Area A)",
ylab="P-value (Enriched in points in Area B)",
axes=F,col=cols[is.non0])
points(0,0,pch=20)
cr <- mget(names(head(sort(nlps[,1],decreasing=T),7)),org.Hs.egSYMBOL)
sh <- mget(names(head(sort(nlps[,2],decreasing=T),7)),org.Hs.egSYMBOL)

text(seq(80,125,length=7),10,unlist(rev(cr)),pos=4,srt=90)
text(1,nlps[head(order(nlps[,2],decreasing=T),7),2],unlist(sh),pos=4)
box()

labelsX=c("1",parse(text=paste0("10^-",seq(20,120,20))))
axis(1,at=seq(0,120,20),labels=labelsX,cex.axis=.8)
axis(2,at=seq(0,120,20),labels=labelsX,cex.axis=.8,las=2)
abline(v=c(124.45,127.38))
abline(h=3,v=3)
dev.off()


is.crispr <- nlps[,1] > 3 ## 958
is.shrna <- nlps[,2] > 3 ## 23

table(crispr=is.crispr,shrna=is.shrna)

ess.cr <- names(sort(nlps[is.crispr,1],decreasing=T))
ess.sh <- names(sort(nlps[is.shrna,2],decreasing=T))

head(sort(nlps[,2],decreasing=T),7)



## Fisher's test - pathway enrichment
setwd(rda_dir)
load(file="msigdb.v7.0.rda") # msigdb

## find overlap
msigdb.ol <- lapply(msigdb,function(x)x[x %in% gns])

## 
n <- sapply(msigdb.ol,length)
n1 <- n >= 15 & n <= 200

msigdb.1 <- msigdb.ol[n1] ## 6551
ns <- sapply(msigdb.1,length)
i <- order(ns,decreasing=T)
msigdb.1 <- msigdb.1[i]

table(sapply(names(msigdb.1),function(x)sub("_.+","",x)))

##
used.genes <- unique(unlist(msigdb.1)) ## 14831
used.crispr <- ess.cr[ess.cr %in% used.genes] # 920/923
used.shrna <- ess.sh[ess.sh %in% used.genes] # 20/23

pvals.path <- sapply(msigdb.1,function(x){
	mat.a <- table(path=used.genes %in% x,crispr=used.genes %in% used.crispr)[2:1,2:1]
	mat.b <- table(path=used.genes %in% x,crispr=used.genes %in% used.shrna)[2:1,2:1]

	p.a <- fisher.test(mat.a,alternative="greater")$p.value
	p.b <- fisher.test(mat.b,alternative="greater")$p.value	

	return(c(crispr=p.a,shrna=p.b))
})

nlps.path <- -apply(pvals.path,1,log10) ## 5976

##
ns <- sapply(msigdb.1,length)

nns <- ns/max(ns)+1
uniq.cols <- colorRampPalette(rev(brewer.pal(9,"YlOrRd")[-(1:3)]))(20)
# cols <- uniq.cols[round(((nns-min(nns))/diff(range(nns)))*19+1)]

db <- names(table(sapply(names(msigdb.1),function(x)sub("_.+","",x))))
dbs <- match(sapply(names(msigdb.1),function(x)sub("_.+","",x)),db)


##
setwd(plot_dir);setwd("01_genewise_relationship")

##
is.b <- nlps.path[,1]>10 & nlps.path[,2] > 2
is.c <- nlps.path[,1]>20 & nlps.path[,2] < 2
is.s <- nlps.path[,1]<10 & nlps.path[,2] > 3

sig.discrepant.paths <- data.frame(path=rownames(nlps.path),nlps.path) %>%
	filter(is.b|is.c|is.s) %>%
	arrange(desc(shrna),(crispr))

ind <- match(sig.discrepant.paths$path,names(msigdb.1))

uniq.cols <- brewer.pal(6,"Spectral")

cols <- c("grey60",uniq.cols[c(1,6,5)])[is.c + is.s*2 + is.b*3 + 1]


## Fig. 1D
pdf("pathway_enrichment_discrepant_genes.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
is.non0 <- nlps.path[,1]>0.1  | nlps.path[,2]>0.1 
is.used <- dbs %in% c(2,4,6) & is.non0
# cols <- densCols(nlps.path[is.used,],colramp = colorRampPalette(blues9[-(1:3)]))
plot(nlps.path[is.used,],pch=20,col=cols[is.used],
	xlab="P-value (Essential pathway, claimed by CRISPR)",
	ylab="P-value (Essential pathway, claimed by shRNA)",
	axes=F)#,xlim=0:1,ylm=0:1)

text(nlps.path[ind,],pch=20,pos=4)

points(0,0,pch=20,col="grey60")
box()

# plot(nlps.path[is.go,],pch=20,cex=nns[is.go],col=cols[is.go])

labelsX=c("1",parse(text=paste0("10^-",seq(10,40,10))))
labelsY=c("1",parse(text=paste0("10^-",1:4)))
axis(1,at=seq(0,40,10),labels=labelsX,cex.axis=.8,las=1)
axis(2,at=0:4,labels=labelsY,cex.axis=.8,las=1)

##
dev.off()


## Fig. 1D - table
write.csv(sig.discrepant.paths,file="discrepancy-pathwas.csv")

##
setwd(rda_dir)
save(ess.cr,ess.sh,file="discrepant_essential_genes.rda")


##
