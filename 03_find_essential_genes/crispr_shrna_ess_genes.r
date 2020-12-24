setwd(rda_dir)
load(file="discrepant_essential_genes.rda") # ess.cr,ess.sh
x <- load(file="find_essential_genes.rda") #is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns

mats.ol <- sapply(1:6,function(k){
	c1 <- seq(0,1,.2)[k]

	setwd(rda_dir);setwd(paste0("r",c1))
	dep.scores <- readRDS("scores_15847g_423c.rds") 

	gns <- rownames(dep.scores)
	ess.genes <- gns[is.ess[[k]]$q1]

	mat.cr <- table(is.ess=gns %in% ess.genes,crispr=gns %in% ess.cr)[2:1,2:1]
	mat.sh <- table(is.ess=gns %in% ess.genes,shrna=gns %in% ess.sh)[2:1,2:1]	

	return(c(crispr=mat.cr[1],shrna=mat.sh[1]))

	nlp.cr <- -log10(fisher.test(mat.cr,alternative="greater")$p.value)
	nlp.sh <- -log10(fisher.test(mat.sh,alternative="greater")$p.value)

	return(c(crispr=nlp.cr,shrna=nlp.sh))
})

##
library(RColorBrewer)
ratio <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))

cols <- brewer.pal(6,"Spectral")
mr.ol <- mats.ol/c(length(ess.cr),length(ess.sh))
colnames(mr.ol) <- ratio

setwd(plot_dir);setwd("02_compare_efficacy_measure")
pdf(paste0("overlap_ratio.pdf"),width=4.5,height=3)
par(mar=c(5,5,3,8))
barplot(mr.ol,beside=T,col=cols[c(1,6)],
	ylab="Relative overlap",las=2,space=c(.1,.8),
	border=NA)
abline(h=seq(0,1,.2),col="grey80")
par(xpd=T)
legend(par()$usr[2],.8,c("CRISPR","shRNA"),fill=cols[c(1,6)],
	border=NA,cex=.8,title="Essential by")
par(xpd=F)
dev.off()

if(0){
	df <- data.frame(ratio=as.vector(mr.ol),
		method=factor(rep(c("CRISPR","shRNA"),times=6),levels=c("CRISPR","shRNA")),
		mix=factor(rep(ratio,each=2),levels=ratio))

	ggplot(data=df, aes(x=mix, y=ratio, fill=method)) +
	  geom_bar(stat="identity", position=position_dodge())+
	  scale_fill_manual(values=cols[c(1,6)])+
	  theme_minimal()
}

mat.n <- sapply(1:6,function(i){
	sapply(1:6,function(j){
		n.i <- is.ess[[i]]$q1
		n.j <- is.ess[[j]]$q1
		ol <- sum(n.i & n.j)
		return(ol)
	})
})

mat.ratio <- sapply(1:6,function(i){
	sapply(1:6,function(j){
		n.i <- is.ess[[i]]$q1
		n.j <- is.ess[[j]]$q1
		ol <- sum(n.i & n.j)
		r <- ol/min(sum(n.i),sum(n.j))
		return(r)		
	})
})

rownames(mat.ratio) <- colnames(mat.ratio) <- ratio
# library(gplots)
# heatmap.2(mat.ratio,Colv=F,Rowv=F,trace="none",dendrogram="none")

mat.ratio[upper.tri(mat.ratio)] <- NA
mat.n[upper.tri(mat.n)] <- NA

##
setwd(plot_dir);setwd("02_compare_efficacy_measure")
pdf("overlap_index.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
# cols <- colorRampPalette(brewer.pal(9,"YlOrRd")[1:7])(20)
image(1:6,1:6,t(mat.ratio[6:1,]),zlim=c(0.5,1),
	col=heat.colors(20),axes=F,
	xlab="",ylab="")
par(tcl=0)
axis(1,at=1:6,labels=ratio,las=2)
axis(2,at=6:1,labels=ratio,las=1)

for(x in 1:6){
	for(y in 1:6){
		if(y <= x){
			text(y,7-x,mat.n[x,y])
		}
		# text(y,7-x,round(mat.ratio[x,y],3))
	}
}
for(i in 1:6){
	j <- 7-i
	lines(c(-1,-1,1,1,-1)*.5+i,c(-1,1,1,-1,-1)*.5+j)
	lines(rep(0.5,2),c(0.5,6.5))
	lines(c(0.5,6.5),rep(0.5,2))
}

dev.off()

##
pdf("color_key.pdf",width=5,height=3)
par(mar=c(5,5,3,3))
image(1:20,1,matrix(1:20,ncol=1),col=heat.colors(20),axes=F,
	xlab="",ylab="")
box()
axis(1,at=seq(1,20,length=6),labels=seq(0.5,1,.1),las=1)
dev.off()



