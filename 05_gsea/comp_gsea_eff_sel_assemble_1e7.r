library(RColorBrewer)
library(gplots)

setwd(rda_dir)
x <- load(file="gsea_data_6thres.rda") # df,msigdb.1

i1 <- names(dfs)
i2 <- names(dfs[[1]])
ids <- expand.grid(idx1=i1,idx2=i2,ef.sel=1:2)

## in O2
setwd(rda_dir);setwd("fgsea")
fns <- paste0("fgsea_1e7_new_",1:72,".rds")
fns %in% dir()

slps  <- list()
for(i in seq(nrow(ids))){
	fn <- fns[[i]]
	res <- readRDS(fn)
	slps[[i]] <- sign(res$ES) * -log10(res$pval)
}

slps.df <- list()
for(i in i1){
	slps.df[[i]] <- list()
	for(j in i2){
		idx <- which(ids$idx1==i & ids$idx2==j)
		o <- ids$ef.sel[idx]
		tmp <- data.frame(do.call(cbind,slps[idx[o]]))
		names(tmp) <- c("eff","sel")
		slps.df[[i]][[j]] <- tmp
	}
}

eff.rng <- range(sapply(slps.df,function(x)sapply(x,function(y)y$eff)))
sel.rng <- range(sapply(slps.df,function(x)sapply(x,function(y)y$sel)))

if(0){
	par(mfrow=c(2,3))
	for(i in i1){
		smoothScatter(slps.df[[i]]$q1,xlim=eff.rng,ylim=sel.rng,main=i)	
		abline(h=0,v=0,lty=2)
	}
}

eff.sel.paths <- lapply(1:6,function(i){
	d <- slps.df[[i]]$q1
	names(msigdb.1)[which(d$eff < -6.9 & d$sel > 6.9)]
})

sapply(eff.sel.paths,function(x){
	sapply(eff.sel.paths,function(y){
		length(intersect(x,y))
	})
})

## Fig. 3E
for(i in i1){

	slp.df <- slps.df[[i]]$q1

	p1 <- which(slp.df$eff < -6.999 & slp.df$sel > 6.999)
	p2 <- which(slp.df$eff < -6.999 & abs(slp.df$sel) < .25) # ef specific
	p3 <- which(abs(slp.df$eff) < .25 & slp.df$sel > 6.999) # sel specific

	##
	setwd(plot_dir);setwd("05_gsea")
	pdf(paste0("smoothscatter_pathway-pvals_",i,".pdf"),width=5,height=5)

	par(mar=c(5,5,3,3))
	plot(slp.df$eff,slp.df$sel,
		pch=NA,main="",type="n",
		xlab="Signed log-P (efficacy)",
		ylab="Signed log-P (selectivity)")
	cols <- densCols(slp.df,colramp = colorRampPalette(blues9[-(1:3)]))
	points(slp.df,col=cols,pch=20,cex=.5)
	abline(h=0,v=0)

	# }
	polygon(-7+c(1,1,-1,-1,1)*.1,7+c(-1,1,1,-1,-1)*.1,border=2,lwd=2)
	polygon(-7+c(1,1,-1,-1,1)*.1,0+c(-1,1,1,-1,-1)*.25,border=2,lwd=2)
	polygon(-0+c(1,1,-1,-1,1)*.25,7+c(-1,1,1,-1,-1)*.1,border=2,lwd=2)

	text(c(-7,-7,0),c(6.5,-0.2,6.5),
		paste(c("selective and essential","essential","selective"),
			"\n(",c(length(p1),length(p2),length(p3)),")",sep=""),
		col=2,pos=4)
	dev.off()

	##
	setwd(plot_dir);setwd("05_gsea")
	write.csv(cbind(names(msigdb.1)[c(p1,p2,p3)],slp.df[c(p1,p2,p3),]),
		file=paste0("significant_pathways_",i,".csv"))

	##
	slp.df <- slps.df[[i]]$q1

	p1 <- which(slp.df$eff < -6.999 & slp.df$sel > 6.999)
	p2 <- which(slp.df$eff < -6.999 & abs(slp.df$sel) < .25) # ef specific
	p3 <- which(abs(slp.df$eff) < .25 & slp.df$sel > 6.999) # sel specific

	##
	sideCols <- brewer.pal(9,"YlGnBu")
	cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))(10)

	## Fig. S3B
	setwd(plot_dir);setwd("05_gsea")

	n1 <- unlist(sapply(msigdb.1[p1],length))
	n2 <- unlist(sapply(msigdb.1[p2],length))
	n3 <- unlist(sapply(msigdb.1[p3],length))

	min.n <- min(c(n1,n2,n3))
	rng <- diff(range(c(n1,n2,n3)))

	sc1 <- sideCols[round((n1-min.n)/rng*8)+1]
	sc2 <- sideCols[round((n2-min.n)/rng*8)+1]
	sc3 <- sideCols[round((n3-min.n)/rng*8)+1]

	if(!is.null(n1)){
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
		pdf(paste0("heatmap_overlap_1_",i,"_v3.pdf"),width=5,height=5) 
		heatmap.2(ois1,trace="none",margins=c(12,12),cexRow=.4,cexCol=.4,
			col=cols,RowSideColors=sc1,ColSideColors=sc1,main="essential & selective")
		dev.off()
	}
	##
	if(!is.null(n2)){
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
		pdf(paste0("heatmap_overlap_2_",i,"_v3.pdf"),width=5,height=5) 
		heatmap.2(ois2,trace="none",margins=c(12,12),cexRow=.3,cexCol=.3,
			col=cols,RowSideColors=sc2,ColSideColors=sc2,main="essential")
		dev.off()
	}

	## selective; not so essential
	if(!is.null(n3)){
		pdf(paste0("heatmap_overlap_3_",i,"_v3.pdf"),width=5,height=5)

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
	}
}
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

##
library(openxlsx)
setwd(plot_dir);setwd("05_gsea")
sig.paths <- lapply(1:6,function(i){
	tmp <- read.xlsx("significant_pathways.xlsx",i)	
	tapply(tmp$Pathway,tmp$Gene.set,function(x)x)
})

n.sig.paths <- sapply(sig.paths,function(sg){
	st <- rep(0,3)
	names(st) <- c("Essential","Selective","Selective & Essential")
	tmp <- sapply(sg,function(x){
		length(x)
	})
	st[names(tmp)] <- tmp
	return(st)
})
colnames(n.sig.paths) <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))

cols <- brewer.pal(3,"Pastel2")[c(3,1,2)]

setwd(plot_dir);setwd("08_comparison")
pdf("sig-path.pdf",width=5,height=4)
par(mar=c(5,5,3,3))
barplot(n.sig.paths,beside=T,col=cols,las=2,
	ylab="# overrepresented pathways")
legend("topright",
	c("Essential","Selective","Selective & Essential"),
	fill=cols,cex=.8)
dev.off()

eff.paths <- lapply(sig.paths,function(sp)sp$Essential)
sel.paths <- lapply(sig.paths,function(sp)sp$Selective)
es.paths <- lapply(sig.paths,function(sp)sp[["Selective & Essential"]])

eff.al <- unique(unlist(eff.paths))
sel.al <- unique(unlist(sel.paths))
es.al <- unique(unlist(es.paths))

##
eff.1 <- sapply(eff.paths,function(x) as.numeric(eff.al %in% x))
is.ge3 <- rowSums(eff.1)>=3
eff.2 <- apply(eff.1,1,paste,collapse="")
pat <- names(table(eff.2))
ok.pat <- pat[!grepl("10+1",pat)]
is.ok <- eff.2 %in% ok.pat
tab.pat <- table(eff.2[is.ge3 & is.ok])
pat2path <- lapply(names(tab.pat),function(tp){
	eff.al[eff.2==tp]
})
names(pat2path) <- names(tab.pat)

##
sel.1 <- sapply(sel.paths,function(x) as.numeric(sel.al %in% x))
limma::vennDiagram(sel.1[,4:6])
is.eq3 <- rowSums(sel.1)==3
sel.al[is.eq3]

es.1 <- sapply(es.paths,function(x) as.numeric(es.al %in% x))
limma::vennDiagram(es.1[,4:6])
is.eq3 <- rowSums(es.1)==3
es.al[is.eq3]

setwd(plot_dir);setwd("08_comparison")
df.eff <- data.frame(eff.al,eff.1)
df.sel <- data.frame(sel.al,sel.1)
df.es <- data.frame(es.al,es.1)
ratio <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))
names(df.es) <- names(df.sel) <- names(df.eff) <- c("Pathway",ratio)
write.csv(df.eff,file="eff_overlap.csv")
write.csv(df.sel,file="sel_overlap.csv")
write.csv(df.es,file="es_overlap.csv")


