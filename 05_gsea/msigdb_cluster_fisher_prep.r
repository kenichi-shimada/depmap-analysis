## this file is to prepare two objects to compute GSEA in the cloud, namely
# - df: efficacy-selectivity scores
# - msigdb.1: list object contains MSigDB pathway data

setwd(rda_dir)
load(file="msigdb.v7.0.rda") # msigdb
load(file="eff-sel_6thres.rda") # coefs,ef.sels

##
all.genes <- unique(unlist(msigdb)) ## 19842

##
exp.genes <- rownames(ef.sels[[1]][[1]])
msigdb.ol <- lapply(msigdb,function(x)x[x %in% exp.genes])

##
n <- sapply(msigdb.ol,length)
n1 <- n >= 15 & n <= 500

msigdb.1 <- msigdb.ol[n1] ## 6551
table(sapply(names(msigdb.1),function(x)sub("_.+","",x)))

n.msigdb.1 <- sapply(msigdb.1,length)
##
used.genes <- unique(unlist(msigdb.1)) ## 14831
dfs <- lapply(ef.sels,function(x)
	lapply(x,function(y)y[exp.genes %in% used.genes,])) ## 14831 x 2


library(dplyr)

mems.ge15 <- lapply(1:6,function(i){
	setwd(rda_dir)
	x <- load(paste0("depmap_initial_19q3_local_run_v3_",i,".rda"))
	mems.1 <- mems %>% filter(eid %in% used.genes)
	mem3eid <- tapply(mems.1$eid,mems.1$mem3,function(x)x)
	is.ge15 <- sapply(mem3eid,length) >= 15

	# table(is.ge15)[c("TRUE","FALSE")]
	mem3eid.ge15 <- mem3eid[is.ge15]
	return(mem3eid.ge15)	
})

mems.all <- lapply(1:6,function(i){
	setwd(rda_dir)
	x <- load(paste0("depmap_initial_19q3_local_run_v3_",i,".rda"))
	return(mems)
})
eids.all <- sapply(mems.all,function(x)x$eid)
table(table(unlist(eids.all))) ## 2008 genes common
comm.eids <- names(which(table(unlist(eids.all))==6)) ## 2008 genes common
mems.mat <- sapply(mems.all,function(m){
	m$mem3[match(comm.eids,m$eid)]
})
head(mems.mat)

library(clue)

mems.partition <- lapply(1:6,function(k){
	as.cl_partition(mems.mat[,k])
})
ratio <- rev(c("100:0\n(CRISPR)","80:20","60:40","40:60","20:80","0:100\n(shRNA)"))
names(mems.partition) <- ratio

setwd(plot_dir);setwd("08_comparison")
pdf("dend_clusters.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
hc <- hclust(cl_dissimilarity(mems.partition,method="comemberships"),method="complete")
plot(hc,hang=0.1,main="Dendrogram of 6 cluster sets\n(2008 common genes)")
dev.off()

pre <- rep(paste0("G",1:6),sapply(mems.ge15,length))
mems.ge15 <- do.call(c,mems.ge15)
names(mems.ge15) <- paste0(pre,names(mems.ge15))

n.mems <- sapply(1:6,function(i){
	setwd(rda_dir)
	x <- load(paste0("depmap_initial_19q3_local_run_v3_",i,".rda"))
	mems.1 <- mems %>% filter(eid %in% used.genes)
	mem3eid <- tapply(mems.1$eid,mems.1$mem3,function(x)x)
	mem3eid.1 <- tapply(mems.1$eid,mems.1$mem3,function(x)x)
	is.ge15 <- sapply(mem3eid.1,length) >= 15

	# table(is.ge15)[c("TRUE","FALSE")]
	mem3eid.ge15 <- mem3eid[is.ge15]
	return(sapply(mem3eid.ge15,length))
})

sum(sapply(n.mems,function(x)sum(x>=15)))

##
length(mems.ge15) # 274 groups

##
if(0){
	nlps <- parallel::mclapply(mems.ge15,function(m){
		cat("*")
		sapply(msigdb.1,function(path){
			mat <- table(is.mem=factor(used.genes %in% m,levels=c("TRUE","FALSE")),
				is.path=factor(used.genes %in% path,levels=c("TRUE","FALSE")))
			nlp <- -log10(fisher.test(mat,alternative="greater")$p.value)
		})
	},mc.cores=10)

	nlps.all <- do.call(rbind,nlps)

	setwd(rda_dir)
	saveRDS(nlps.all,"nlps_all_cluster_fisher.rds")
}else{
	setwd(rda_dir)
	nlps.all <- readRDS("nlps_all_cluster_fisher.rds")
}

##
idx.nlps <- apply(nlps.all,1,function(nlp){
	which.max(nlp)
})
max.nlps <- apply(nlps.all,1,function(nlp){
	max(nlp)
})

##
max.paths <- colnames(nlps.all)[idx.nlps]
n.paths <- n.msigdb.1[max.paths]
max.paths <- tolower(gsub("_"," ",sub("[^_]+_","",max.paths)))
gps <- as.numeric(sub("G(.).+","\\1",rownames(nlps.all)))
gps <- rev(c("100:0 (CRISPR)","80:20","60:40","40:60","20:80","0:100 (shRNA)"))[as.numeric(gps)]

cls <- sub("..(.+)","\\1",rownames(nlps.all))

##
df <- data.frame(mixture=gps,Cluster=cls,
	n=unlist(n.mems),
	Pathway=max.paths,#n=n.paths,
	nlp=max.nlps,
	stringsAsFactors=F)
rownames(df) <- c()
names(df) <- c("Score","Cluster","Num.Genes","Pathway","-logP")

setwd(plot_dir);setwd("08_comparison")
write.csv(df,"cluster_pathway_enrichment_v3.csv")

paths <- unique(df$path)
gr.path <- tapply(df$path,df$group,function(x)x)
sapply(gr.path,length)

sapply(gr.path,function(x)x[1:21])

mat.gr.path <- sapply(gr.path,function(x)match(paths,x))

gr.cl <- tapply(df$clu,df$group,function(x)x)
table(mat.gr.path)
head(mat.gr.path)
g6.path <- paths[which(rowSums(mat.gr.path)==6)]
gr.path <- tapply(df$path,df$group,function(x)x[x %in% g6.path])


sig.g4 <- df[idx,] %>% filter(group=="G4" & nlp >= 30)

g4.cl <- sig.g4$cluster
g4.path <- sig.g4$path

tapply(max.paths[idx],gps[idx],function(x)x)

##
setwd(rda_dir)
x <- load("depmap_initial_19q3_local_run_v3_4.rda")

eids.4 <- mems %>% 
	filter(mem3 %in% g4.cl) %>%
	mutate(mem3 = factor(mem3)) %>%
	select(eid,mem3)

mem2eid <- tapply(eids.4$eid,eids.4$mem3,identity)

eids.ess <- mems$eid
cols <- rep("grey80",length(eids.ess))
pchs <- rep(20,length(eids.ess))
cexs <- rep(.5,length(eids.ess))
names(cols) <- eids.ess
nl <- length(mem2eid) # 9

set.seed(1)
# uniq.cols <- sample(colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral")[-(4:5)])(9))
# uniq.cols <- sample(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(nl))
uniq.cols <- RColorBrewer::brewer.pal(10,"Set3")[-9]

for(l in seq(mem2eid)){
	cols[mem2eid[[l]]] <- uniq.cols[l]
	cexs[mem2eid[[l]]] <- 1
}
xys <- t(sapply(mem2eid,function(xs){
	idx <- mems$eid %in% xs
	med.x <- median(mems$x[idx])
	med.y <- median(mems$y[idx])
	c(med.x,med.y)
}))


##
setwd(plot_dir);setwd("07_echodots_prep")
pdf("rtsne_map.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
plot(mems$x,mems$y,col=cols,pch=pchs,cex=cexs,
	xlab="t-SNE coordinate 1",ylab="t-SNE coordiante 2")
text(xys,names(mem2eid),cex=1.5)
dev.off()

write.csv(sig.g4,"sig_g4.csv")

##
