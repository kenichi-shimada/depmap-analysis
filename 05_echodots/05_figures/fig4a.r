library(dbscan)
library(RColorBrewer)
library(gplots)

##
source("my_hullplot.r") # or source it manually

setwd(rda_dir)
coef <- readRDS("eff_coef_2492_genes.rds")
cols <- rev(brewer.pal(11,"Spectral"))
load("eff-sel.rda") # ef.sel, thres

setwd(rda_dir);setwd("rtsne")
rt <- readRDS("rtsne_1.rds")


## Fig. 4A (heatmap)
setwd(plot_dir);setwd("06_echodots_prep")
qt99 <- quantile(coef,c(.005,.995))

coef.1 <- coef
coef.1[coef.1 > qt99[2]] <- qt99[2]
coef.1[coef.1 < qt99[1]] <- qt99[1]

png("fig4A_coef_heatmap_2.png",width=500,height=500)
heatmap.2(coef.1,trace="none",col=cols,
	labRow=rep("",nrow(coef.1)),labCol=rep("",nrow(coef.1)),
	key.title=NA,key.xlab="",key.ylab="",
	key.ytickfun=function() {
          return(list(
               at=c(-10000),
               labels=""
               ))
      },
	margins=c(5,5))
dev.off()

## Fig. 4A (grey rtsne map)
pdf("fig4A_plain-rt.pdf",width=4,height=4)
par(mar=c(3,3,2,2))
plot(rt,pch=20,col="grey50",main="t-SNE plot",xlab="",ylab="",cex.main=1.5,las=2,cex=.5)
dev.off()


## Fig. 4A (hullplots)
sig.ds <- c("53","87","131")
set.seed(1)
cols <- sample(rainbow(10))
for(d in sig.ds){
	pdf(paste0("fig4A_hullplots_d-",d,".pdf"),width=5,height=5)
	par(mar=c(0.5,1,1.5,1))

	clue.mem <- function(x=rt,i){
		knn.th <- diff(range(x))/i
		cl <- dbscan(x,eps=knn.th,minPts=2)
		mem <- cl$cl
		return(mem)
	}

	mem <- clue.mem(rt,as.numeric(d))
	tab <- table(mem)
	m2 <- names(which(tab==2))
	is.0 <- mem==0

	my.hullplot(rt[!is.0,],cl=mem[!is.0],pch=20,main=paste0("d = ",d),
		xlab="",ylab="",solid=F,axes=F,cex.main=2,col=cols)
	box()
	dev.off()
}
