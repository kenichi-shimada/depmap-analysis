library(RColorBrewer)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(sp)

setwd(rda_dir)
ess.genes <- readRDS(file="ess_2492_genes.rds")
ess.syms <- unlist(mget(ess.genes,org.Hs.egSYMBOL))

x <- load("eff-sel.rda") # ef.sel, thres

mems <- readRDS(file="hierarchy_cluster.rds")

is.ess <- rownames(ef.sel) %in% ess.genes

eff <- ef.sel$eff[is.ess]
sel <- ef.sel$sel[is.ess]

mem <- mems$m3
mem <- factor(mem,levels=levels(mem)[1:203])

is.used <- mem!="NA"
sel.mem <- tapply(sel[is.used],mem[is.used],identity)
eff.mem <- tapply(eff[is.used],mem[is.used],identity)

sbx <- boxplot(sel.mem,pch=20,cex=.3,lwd=.25)
ebx <- boxplot(eff.mem,pch=20,cex=.3,lwd=.25)

## pick
n.mem <- table(mem[is.used])
st <- sbx$stats
et <- ebx$stats
n <- round(log2(n.mem))
rng.n <- tapply(n.mem,n,function(x){
	rng <- range(x)
	if(rng[1]==rng[2]){
		return(rng[1])
	}else{
		return(paste0(rng[1],"-",rng[2]))
	}
})
uniq.cols <- c("grey90",blues9[-c(1:2)])
cols <- uniq.cols[n]

## Fig. 4F
if(0){
	# Manually picked the clusters to show their labels
	plot(et[3,],st[3,],pch=21,cex=n/2,col=1,bg=cols,type="n",
		xlim=range(et[2:4,]),ylim=range(st[2:4,]),
		xlab="efficacy",ylab="selectivity")		
	points(et[3,],st[3,],pch=21,cex=(n+1)/4,col=1,bg=cols)
	xy <- locator(20,type="l",col=4)
	xy <- locator(8,type="l",col=4)

	idx <- which(point.in.polygon(et[3,],st[3,],xy$x,xy$y)==1)
}else{
	idx <- c(71,84,137,140,145,147,164,173,177,192,194,199,115,138,201)
}

##
setwd(plot_dir);setwd("07_cluster_eff_sel")
pdf("fig4F_scatter-sel-eff_norm_txt.pdf",width=5,height=5)
par(mar=c(4,4,1,1))
plot(et[3,],st[3,],pch=21,cex=n/2,col=1,bg=cols,type="n",
	xlim=range(et[3,]),ylim=range(st[3,]),
	xlab="efficacy",ylab="selectivity")		
abline(h=0,lty=2,col="grey70")
points(et[3,],st[3,],pch=21,cex=(n+1)/4,col=1,bg=cols)
text(et[3,idx],st[3,idx],colnames(et)[idx],pos=3,cex=1)
legend("topleft",rng.n,
	pch=21,pt.cex=unique(sort(n)+1)/4,pt.bg=uniq.cols,
	title="# genes")
dev.off()

## Fig. 4G, selected genes
sel.m3 <- paste0("L",idx)
colnames(et) <- colnames(st) <- paste0("L",seq(ncol(et)))

mems.sel <- mems %>% filter(m3 %in% sel.m3) %>%
	group_by(m1,m2,m3) %>% summarize(syms=paste(sym,collapse=",")) %>%
	mutate(eff=round(et[3,m3],3)) %>%
	mutate(sel=round(st[3,m3],3)) %>%
	arrange(desc(sel),m3,m2,m1)

write.csv(mems.sel[c(3:1,4:6)],"fig4G_sel_genes_112519.csv")
