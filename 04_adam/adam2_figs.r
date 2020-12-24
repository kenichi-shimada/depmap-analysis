library(dplyr)
library(impute)
library(RColorBrewer)
library(ggplot2)
library(org.Hs.eg.db)

setwd(rda_dir)
load(file="eff-sel_6thres.rda") # coefs,ef.sels
load(file="find_essential_genes.rda") # is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns
cl.info <- readRDS("celllines_19q3.rds") 
is.comm <- readRDS(file="is.comm.rds")
load(file="crossoverpoints_bagel_ge10.rda") # crossoverpoints,lin.crossoverpoints
load(file="lineage_cells.rda") # lin,lin2cell,n.lin2cell
load(file="find_essential_genes.rda")

k=4
c1 <- paste0("r",seq(0,1,.2)[k])
setwd(rda_dir);setwd(c1)
sc <- readRDS("scores_15847g_423c.rds")

used.cells <- unlist(lin2cell) # 387
used.sc <- sc[,used.cells] ## 


## 
# n.lin <- sapply(lin2cell,length)
# o <- order(n.lin,decreasing=T)
# par(mar=c(10,5,4,1))
# barplot(n.lin[o],las=2) # 17 tissues

## impute essentiality matrix
is.ess.imp <- parallel::mclapply(1:6,function(k){
  c1 <- seq(0,1,.2)[k]
  setwd(rda_dir);setwd(paste0("r",c1))
  scores <- readRDS(file=paste0("scores_15847g_423c.rds"))
  imp.scores <- impute.knn(scores)$data

  ess.mat <- array(as.numeric(imp.scores < sc.ths[[k]][1]),dim(imp.scores))
  rownames(ess.mat) <- rownames(scores)
  colnames(ess.mat) <- colnames(scores)

  return(ess.mat)
},mc.cores=6)

cl.info <- cl.info %>% filter(DepMap_ID %in% colnames(is.ess.imp[[1]]))

## number of essential genes per lineage
n.lin.ess <- parallel::mclapply(1:6,function(k){
  sapply(lin,function(ti){
    cells <- lin2cell[[ti]]
    this.ess <- is.ess.imp[[k]][,cells,drop=F]
    n.ess <- rowSums(this.ess)
    return(n.ess)
  })
},mc.cores=6)

## # dependent cells (not lineage)
n.dep.cells <- sapply(n.lin.ess,rowSums)

##
is.lin.dep <- lapply(1:6,function(k){
	nle <- n.lin.ess[[k]]
	is.lin.dep <- nle >= crossoverpoints[,k][col(nle)]
	num.lin.dep <- array(as.numeric(is.lin.dep),dim(is.lin.dep))
	sym.all <- unlist(mget(rownames(is.lin.dep),org.Hs.egSYMBOL))
	rownames(num.lin.dep) <- rownames(is.lin.dep)
	colnames(num.lin.dep) <- colnames(is.lin.dep)
	num.lin.dep <- data.frame(symbol=sym.all,num.lin.dep)
	return(num.lin.dep)
})

n.lin.dep <- sapply(1:6,function(k){
	nle <- n.lin.ess[[k]]
	is.lin.dep <- nle >= crossoverpoints[,k][col(nle)]
	return(rowSums(is.lin.dep))
})

setwd(rda_dir)
saveRDS(is.lin.dep,file="is.lin.dep_v3.rds")
saveRDS(n.dep.cells,file="n.dep.cells_v3.rds")
saveRDS(n.lin.dep,file="n.lin.dep_v3.rds")

write.csv(is.lin.dep[[4]],file="is.lin.dep_v3.csv")
##
ef.sel <- ef.sels$r0.6$q1
eff <- ef.sel$eff
sel <- ef.sel$sel



##
setwd(src_dir);source("functions.r")
k=4
ndc <- n.dep.cells[,k]
df.1 <- data.frame(cells=ndc,lins=factor(n.lin.dep[,k]))
sel.th <- find.threshold(ef.sel$sel,side.fit="lower",side.thres="upper")
zero.sel <-  df.1$cells > 100 & df.1$lins==0 & ef.sel$sel > sel.th

all.gns <- rownames(ef.sel)
zero.sym <- unlist(mget(all.gns[zero.sel],org.Hs.egSYMBOL))


this.sc <- used.sc[,used.cells]
cl.sub <- cl.info[c("DepMap_ID","lineage")] %>% 
	filter(DepMap_ID %in% used.cells) %>%
	mutate(lineage = factor(lineage,levels=sort(unique(lineage),decreasing=T)))
cl.sub <- cl.sub[match(cl.sub$DepMap_ID,used.cells),]

## one lineage

one.sel <-  which(df.1$lins==1 & ef.sel$sel > sel.th)
one.gns <- all.gns[one.sel]

nle <- n.lin.ess[[k]]
is.lin.dep <- nle >= crossoverpoints[,k][col(nle)]
nle <- as.matrix(apply(is.lin.dep[one.sel,],1,function(x)colnames(is.lin.dep)[which(x)]))

med.scores.1 <- sapply(one.gns,function(gn){
	this.lins <- nle[gn,]
	this.cells <- unlist(lin2cell[this.lins])
	other.cells <- unlist(lin2cell[!names(lin2cell) %in% this.lins])
	med.dif <- median(used.sc[gn,this.cells])-median(used.sc[gn,other.cells])
})

diff.sc <- list()
diff.sc[[1]] <- list(sc=med.scores.1,lin=nle)

## k-lineage specific
for(k in 2:14){
	cat("*")
	this.sel <-  which(df.1$lins==k & ef.sel$sel > sel.th)
	this.gns <- all.gns[this.sel]

	nle <- t(apply(is.lin.dep[this.sel,],1,function(x)colnames(is.lin.dep)[which(x)]))
	colnames(nle) <- paste0("lin",seq(ncol(nle)))

	med.scores <- sapply(this.gns,function(gn){
		this.lins <- nle[gn,]
		this.cells <- unlist(lin2cell[this.lins])
		other.cells <- unlist(lin2cell[!names(lin2cell) %in% this.lins])
		med.dif <- median(used.sc[gn,this.cells],na.rm=T)-median(used.sc[gn,other.cells],na.rm=T)
	})
	names(med.scores) <- this.gns

	diff.sc[[k]] <- list(sc=med.scores,lin=nle)
}

##
boxplot(lapply(diff.sc,function(x)x$sc))
abline(h=-0.3)
sig.nle <- lapply(diff.sc,function(x){
	idx <- which(x$sc < -0.3)
	x$lin[idx,,drop=F]
})

selected.gns <- c(leukemia="4602",
	colorectal="1499",
	mm="3662",
	pancreas="3845",
	skin="6663",
	lm="865",
	ckos="10979",
	all="5693",
	all.2="5682",
	zero.sel="1019",
	zero.0="3670",
	zero.1="22931")

zero.sel <-  df.1$cells > 50 & df.1$lins==0 & ef.sel$sel > sel.th
zero.sel.gns <- all.gns[which(zero.sel)]
o <- order(ef.sel$ef[which(zero.sel)],decreasing=F)
zero.sel.gns <- zero.sel.gns[o]
head(unlist(mget(zero.sel.gns,org.Hs.egSYMBOL)))

zero.0 <-  df.1$cells > 50 & df.1$lins==0 & ef.sel$sel < sel.th
zero.0.gns <- all.gns[which(zero.0)]
o <- order(ef.sel$ef[which(zero.0)],decreasing=F)
zero.0.gns <- zero.0.gns[o]
head(unlist(mget(zero.0.gns,org.Hs.egSYMBOL)))

##
selected.syms <- unlist(mget(selected.gns,org.Hs.egSYMBOL))
nle.sel <- apply(is.lin.dep[selected.gns,],1,sum)
n.lin.sum <- apply(n.lin.ess[[k]][selected.gns,],1,sum)

syms <- c("ISL1","CDK4","CTNNB1","MYB","IRF4","KRAS","SOX10","FERMT2","PSMB5","PSMA1")
o1 <- match(syms,selected.syms)
setwd(plot_dir);setwd("06_adam")
write.csv(
	t(cbind(selected.syms,
		round(ef.sel[selected.gns,],3),
		nle.sel,
		n.lin.sum)[o1,]),
	"selected.gns.stat.csv")




sc.rng <- range(sc[selected.gns,used.cells])
mget("KIF11",org.Hs.egSYMBOL2EG)
##
setwd(plot_dir);setwd("06_adam")
for(gn in selected.gns){
	gn <- "3832"
	sym1 <- unlist(mget(gn,org.Hs.egSYMBOL))
	pdf(paste0("violin_",gn,".pdf"),width=3,height=5)
	pvio <- data.frame(ef=used.sc[gn,used.cells],DepMap_ID=used.cells,
		stringsAsFactors=F) %>%
		left_join(cl.sub,by="DepMap_ID") %>%
		# mutate(lineage=factor(lineage,levels=cell.levs)) %>%
		ggplot(aes(x=lineage,y=ef)) +
		geom_hline(yintercept=c(0,sc.ths$r0.6[1]), col="grey70",size=0.2) +
		geom_violin(scale="area",aes(fill=lineage,color=lineage)) +
		geom_jitter(height = 0, width = 0.2, size=.5) +
		theme(plot.margin = unit(c(0,0.2,0,0), "lines")) +
		# coord_cartesian(xlim = c(10, 14), clip = 'off') +
		ylim(sc.rng) + 
		coord_flip() +
		theme_bw() +
		theme(legend.position="none",plot.margin=unit(c(0,.2,0,0),"lines"),
		  axis.text=element_text(size=8),axis.title=element_text(size=10),
		  plot.title=element_text(size=15, hjust=0.5),
  		  panel.grid.minor = element_line(size = 0.2),
		  panel.grid.major = element_line(size = 0.2),
		  axis.text.x = element_text(angle = 0)) +
		ggtitle(paste0(sym1))
	print(pvio)
	dev.off()
}

##
length(ndc)

##
setwd(plot_dir);setwd("06_adam")
pdf("ecdf_n_dep-cells.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
cols <- rev(brewer.pal(6,"Spectral"))
plot.ecdf(n.dep.cells[,1],col=cols[1],verticals=T,pch=NA,
	xlab="# dependent cell lines",
	ylab="Empirical CDF")
for(i in 2:6){
	plot.ecdf(n.dep.cells[,i],col=cols[i],verticals=T,add=T,pch=NA)	
}
ratio <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))
legend(200,0.6,ratio,lty=1,col=cols)
dev.off()

##
setwd(plot_dir);setwd("06_adam")
pdf("n_dep-cells_p1.pdf",width=5,height=5)
par(mar = c(5,5,3,3))
min.n <- round(423*0.01)
cols <- c(rep("grey70",min.n),
	colorRampPalette(brewer.pal(11,"Spectral"))(max(ndc)-min.n+1))
o <- order(ndc)
plot(ef.sel[o,],pch=20,col=cols[ndc[o]+1],type="n",
	xlab="Efficacy",ylab="Selectivity",
	main="")#,cex.lab=3,cex.axis=3)
abline(h=0,v=0,lty=1,col="grey80")
abline(v=sc.ths$r0.6[["0.001"]],col="grey20",lty=2)
points(ef.sel[o,],pch=20,col=cols[ndc[o]+1],cex=.5)#,cex=1.5)

dev.off()

##
pdf("n_dep-cells_p2.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
image(min.n:max(ndc),1,matrix(min.n:max(ndc),ncol=1),col=cols[-seq(min.n)],axes=F,
	xlab="",ylab="")
box()
axis(1,at=c(min.n+1,seq(75,300,75),387),las=1)
dev.off()


##
pdf("n_dep-lin_barplot_p1.pdf",width=5,height=5)
# barplot(table(n.lin.dep))

par(mar = c(5,5,3,3))
cols <- c("grey70",colorRampPalette(brewer.pal(11,"Spectral"))(max(n.lin.dep[,k])))
o <- order(n.lin.dep[,k])
plot(ef.sel[o,],pch=20,col=cols[n.lin.dep[,k][o]+1],type="n",
	xlab="Efficacy",ylab="Selectivity",
	main="")#,cex.lab=3,cex.axis=3)
abline(h=0,v=0,lty=1,col="grey80")
abline(v=sc.ths$r0.6[["0.001"]],col="grey20",lty=2)
points(ef.sel[o,],pch=20,col=cols[n.lin.dep[,k][o]+1],cex=.5)#,cex=1.5)

dev.off()

##
pdf("n_dep-lin_barplot_p2.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
image(1:max(n.lin.dep[,k]),1,matrix(1:max(n.lin.dep[,k]),ncol=1),col=cols[-1],axes=F,
	xlab="",ylab="")
box()
axis(1,at=c(1,5,10,15,17))
dev.off()

##
library(dplyr)
library(ggplot2)
setwd(plot_dir);setwd("06_adam")
for(k in 1:6){
	ndc <- n.dep.cells[,k]
	df.1 <- data.frame(cells=ndc,lins=factor(n.lin.dep[,k]))
	# boxplot(df.1$cells ~ df.1$lins)

	pdf(paste0("violin_lin_cells_",k,".pdf"),width=5,height=5)
	nb.cols <- nlevels(df.1$lins)
	mycolors <- c("grey20",colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols-1))
	# Create a ggplot with 18 colors 
	# Use scale_fill_manual
	p <- df.1 %>% 
		ggplot(aes(x=lins,y=cells)) + 
		geom_violin(scale="width",aes(fill=lins),size=.5) +
		scale_fill_manual(values=mycolors) +
		geom_vline(aes(xintercept=lin.crossoverpoints[k]+.5),linetype='dashed')+
		ggtitle(paste0("k = ",k)) + 
		theme_bw() + 
	  	theme(legend.position="none",
	  		panel.grid.minor = element_line(size = 0.2),
			panel.grid.major = element_line(size = 0.2),
			plot.margin = margin(20,20,12,12))+
		xlab("# depending lineages") + ylab ("# depending cells")
	print(p)
	dev.off()
}

##
pdf("n_lin_barplot.pdf",width=5,height=5)
tabs <- list()
for(i in 1:6){
	nle <- n.lin.ess[[i]]
	is.lin.dep <- nle >= crossoverpoints[,i][col(nle)]
	n.lin.dep <- rowSums(is.lin.dep)
	tab <- table(n.lin.dep)
	tab <- tab[-1]
	tabs[[i]] <- tab
	cols <- colorRampPalette(brewer.pal(11,"Spectral"))(max(n.lin.dep))
	barplot(tab,col=cols,ylim=c(0,1000),axes=F,
		ylab="# essential genes")
	axis(2,at=c(0,500,1000))
}
dev.off()

##
tab <- table(n.lin.dep)
tab <- tab[-1]
dev.off()

##
rbind(com=colSums(is.comm),sel=sapply(tabs,sum)-colSums(is.comm))
par(mar=c(5,5,3,3))
ef.sel <- ef.sels[[k]]$q1
com <- is.comm[,k]==1

plot(ef.sel,main=k,pch=20,type="n",
	xlim=range(ef.sel$eff),ylim=range(ef.sel$sel))
cols <- densCols(ef.sel$eff,ef.sel$sel,
	colramp = colorRampPalette(blues9[-(1:3)]))
points(ef.sel[!com,],col=cols[!com],pch=20,cex=.5)
points(ef.sel[com,],pch=20,col=1,cex=.5)
abline(h=0,v=0,lty=2)

dim(n.lin.ess[[1]])

##
library(org.Hs.eg.db)
k <- 4
sel.th <- find.threshold(ef.sel$sel,side.fit="lower",side.thres="upper")
zero.sel <-  df.1$cells > 50 & df.1$lins==0 & ef.sel$sel > sel.th

plot(ef.sel,main=k,pch=20,type="n",
	xlim=range(ef.sel$eff),ylim=range(ef.sel$sel))
cols <- densCols(ef.sel$eff,ef.sel$sel,
	colramp = colorRampPalette(blues9[-(1:3)]))
points(ef.sel[!zero.sel,],col=cols[!zero.sel],pch=20,cex=.5)
points(ef.sel[zero.sel,],pch=20,col=1,cex=.5)
abline(h=0,v=0,lty=2)
abline(h=sel.th,lty=2)

##

all.gns <- rownames(ef.sel)
unlist(mget(all.gns[zero.sel],org.Hs.egSYMBOL))

sc.ths
