library(RColorBrewer)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(sp)

##
setwd(rda_dir)
x <- load(file="find_essential_genes.rda") 
th <- sc.ths$r0.6[1]
ess.genes <- list()
for(k in c(1:6)){
	c1 <- seq(0,1,.2)[k]
	ess.genes[[k]] <- names(which(is.ess[[k]]$q1))
}

this.ess.genes <- ess.genes[[4]]
##
setwd(rda_dir)
x <- load("eff-sel_6thres.rda") # ef.sel, thres
ef.sel <- ef.sels$r0.6$q1

setwd(rda_dir);setwd("r0.6");
mems <- readRDS("ori_hierarchy_cluster_v3.rds")

is.ess <- rownames(ef.sel) %in% this.ess.genes

eff <- ef.sel$eff[is.ess]
sel <- ef.sel$sel[is.ess]


mem <- mems$m3
# mem[mem=="L339"] <- "NA"
# mem <- factor(mem,levels=levels(mem)[1:203])
lik <- mems$l3

df <- data.frame(selectivity=sel,likelihood=factor(n.lik/10))

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
rng.n <- rng.n[-1]

# uniq.cols <- c("grey90",colorRampPalette(blues9[-c(1:2)])(max(n))) # brewer.pal(9,"YlGnBu")
uniq.cols <- colorRampPalette(blues9[-c(1:2)])(max(n)) # brewer.pal(9,"YlGnBu")
cols <- uniq.cols[n]

## Fig. 4F
if(0){
	# Manually picked the clusters to show their labels
	plot(et[3,],st[3,],pch=21,cex=n/2,col=1,bg=cols,type="n",
		xlim=range(et[3,]),ylim=range(st[3,]),
		xlab="efficacy",ylab="selectivity")		
	points(et[3,],st[3,],pch=21,cex=(n+1)/4,col=1,bg=cols)
	xy <- locator(20,type="l",col=4)
	# xy <- locator(8,type="l",col=4)

	idx <- which(point.in.polygon(et[3,],st[3,],xy$x,xy$y)==1)
}else{
	idx <- c(19,68,91,119,136,142,161,168,185,210,220,243,262,280,281,297,305,309,311,316,318,337)
	idx <- c(68,119,161,168,185,210,243,262,281,297,305,309,311,316,337)

}
sel.m3 <- paste0("L",idx)
colnames(et) <- colnames(st) <- paste0("L",seq(ncol(et)))

setwd(plot_dir);setwd("07_echodots_prep")
pdf("scatter-sel-eff_norm_txt_v3_p1.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
plot(et[3,],st[3,],pch=21,cex=n/2,col=1,bg=cols,type="n",
	xlim=range(et[3,]),ylim=range(st[3,]),
	xlab="efficacy",ylab="selectivity")		
abline(h=0,lty=2,col="grey70")
points(et[3,],st[3,],pch=21,cex=(n+1)/4,col=1,bg=cols)
abline(v=th,lty=2)
dev.off()
##
pdf("scatter-sel-eff_norm_txt_v3_p2.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
plot(et[3,],st[3,],pch=21,cex=n/2,col=1,bg=cols,type="n",
	xlim=range(et[3,]),ylim=range(st[3,]),
	xlab="efficacy",ylab="selectivity",axes=F)	
text(et[3,idx],st[3,idx],colnames(et)[idx],pos=3,cex=1)
dev.off()
##
pdf("scatter-sel-eff_norm_txt_v3_p3.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
plot(et[3,],st[3,],pch=21,cex=n/2,col=1,bg=cols,type="n",
	xlim=range(et[3,]),ylim=range(st[3,]),
	xlab="efficacy",ylab="selectivity",axes=F)		
legend("topleft",rng.n,
	pch=21,pt.cex=unique(sort(n)+1)/4,pt.bg=uniq.cols,
	title="# genes/cluster",border=F)
dev.off()

## Fig. 4G, selected genes

mems.sel <- mems %>% 
	filter(m3 %in% sel.m3) %>%
	dplyr::select(sym,m1,m2,m3) %>%
	group_by(m1,m2,m3) %>% summarize(syms=paste(sym,collapse=",")) %>%
	mutate(eff=round(et[3,m3],3)) %>%
	mutate(sel=round(st[3,m3],3)) %>%
	arrange(desc(sel),m3,m2,m1)

setwd(plot_dir);setwd("07_echodots_prep")
write.csv(mems.sel[c(3:1,4:6)],"sel_genes_092920.csv")
