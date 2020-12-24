library(openxlsx)

##
setwd(rda_dir)
load(file="lineage_cells.rda") # lin,lin2cell,n.lin2cell
is.comm <- readRDS(file="is.comm.rds")

depmap.gns <- rownames(is.comm)

## Behan data
setwd(data_dir);setwd("41586_2019_1103_MOESM4_ESM")
wb <- loadWorkbook("priorityscores_cancers.xlsx")
sheets <- names(wb)

ess.genes <- list()
for(lin in sheets){
	tmp <- as.character(readWorkbook(wb,lin)$entrez_id)
	tmp <- tmp[!is.na(tmp) & tmp %in% depmap.gns]
	ess.genes[[lin]] <- tmp
}

nxy <- sapply(ess.genes,function(x){
	sapply(ess.genes,function(y){
		n <- length(intersect(x,y))
		return(n)
	})
})

rxy <- sapply(ess.genes,function(x){
	sapply(ess.genes,function(y){
		r <- length(intersect(x,y))/min(length(x),length(y))
		return(r)
	})
})

library(gplots)
a <- heatmap.2(rxy,trace="none")

n.lin <- sapply(lin2cell,length)
names(n.lin)

head(n.lin)

lin.names <- c(
	"bone",
	"breast",
	"central_nervous_system",
	"colorectal",
	"esophagus",
	"gastric",
	"leukemia",
	"lung",
	"ovary",
	"pancreas")

behan.names <- c(
	"bone",
	"Breast Carcinoma",
	"Central Nervous System",
	"Colorectal Carcinoma",
	"Esophagus",
	"Gastric Carcinoma",
	"Haematopoietic and Lyphoid",
	"Lung Adenocarcinoma",
	"Ovarian Carcinoma",
	"Pancreatic Carcinoma")

ess.genes.shared <- ess.genes[behan.names]


nle <- n.lin.ess[[4]][,lin.names]

is.lin.dep <- nle >= crossoverpoints[,4][lin.names] & is.comm[,4]==0

depmap.lin.genes <- apply(is.lin.dep,2,function(x){
	names(which(x))
})

behan.lin.genes <- ess.genes[behan.names]
names(behan.lin.genes) <- lin.names

n.depmap <- sapply(depmap.lin.genes,length)
n.behan <- sapply(behan.lin.genes,length)

n.both <- sapply(lin.names,function(x){
	length(intersect(depmap.lin.genes[[x]],behan.lin.genes[[x]]))
})

both.scores <- lapply(lin.names,function(x){
	gns <- intersect(depmap.lin.genes[[x]],behan.lin.genes[[x]])
	tbl <- ess.scores[[x]]
	tbl[tbl$entrez_id %in% gns,]
})
names(both.scores) <- lin.names

##
ess.scores <- list()
for(lin in sheets){
	tmp <- readWorkbook(wb,lin)
	tmp <- tmp[!is.na(tmp$entrez_id),]
	ess.scores[[lin]] <- tmp
}

ess.scores <- ess.scores[behan.names]
names(ess.scores) <- lin.names


##
n.all <- nrow(is.comm)
sum.df <- cbind(both=n.both,
	depmap=n.depmap-n.both,
	behan=n.behan-n.both,
	rest=n.all-n.depmap-n.behan-n.both)


## ecdf, behan tissue-specific genes
if(0){
	rng <- range(unlist(sapply(ess.scores[-1],function(x)x$FINAL.SCORE)))
	plot.ecdf(ess.scores[[1]]$FINAL,pch=NA,verticals=T,xlim=rng)

	for(lin in lin.names[-1]){
		plot.ecdf(ess.scores[[lin]]$FINAL,pch=NA,verticals=T,add=T)
	}

	for(lin in lin.names){
		plot.ecdf(both.scores[[lin]]$FINAL,pch=NA,verticals=T,add=T,col=2)
	}
}

##
depmap.comm <- names(which(is.comm[,4]==1))
behan.comm <-ess.genes[["pan-cancer"]]
sum(behan.comm %in% depmap.gns)
both.comm <- intersect(depmap.comm,behan.comm)
df <-data.frame(depmap=length(depmap.comm),
	behan=length(behan.comm),
	both=length(both.comm))

dim(depmap.comm)

table(names(tab) %in% behan.comm)
table(behan.comm %in% names(tab))

##
setwd(rda_dir)
load(file="eff-sel_6thres.rda") # coefs,ef.sels
ef.sel <- ef.sels$r0.6$q1

tab.ess.sc <- rep(0,n.all)
names(tab.ess.sc) <- rownames(is.comm)
tab <- table(unlist(ess.genes))
tab.ess.sc[names(tab)] <- tab
o <- order(tab.ess.sc)

cols <- c("grey70",colorRampPalette(rev(brewer.pal(11,"Spectral")))(15))

length(ess.genes)
##
setwd(plot_dir);setwd("06_adam")
pdf("n_behan_p1.pdf",width=5,height=5)
par(mar = c(5,5,3,3))
plot(ef.sel[o,],pch=20,type="n",
	xlab="Efficacy",ylab="Selectivity",
	main="")#,cex.lab=3,cex.axis=3)
abline(h=0,v=0,lty=1,col="grey80")
abline(v=sc.ths$r0.6[["0.001"]],col="grey20",lty=2)
points(ef.sel[o,],col=cols[tab.ess.sc[o]+1],pch=20,cex=.5)
dev.off()

##
pdf("n_behan_p2.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
image(1:15,1,matrix(1:15,ncol=1),col=cols[-1],axes=F,
	xlab="",ylab="")
box()
axis(1,at=c(1,5,10,15),las=1)
dev.off()

