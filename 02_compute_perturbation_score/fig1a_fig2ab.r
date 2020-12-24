library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
# library(dplyr)

## load data
setwd(rda_dir)
load(file="eff_dep_19q3.rda") ## eff.1, dep.1, both 15847 x 423 cells
ed1 <- as.data.frame(readRDS("non-missing-ed.rds"))

## Fig. 1B
setwd(src_dir);source("functions.r")
setwd(plot_dir);setwd("01_genewise_relationship")

syms <- c("CCND1","RAN","EIF5B","FOXD4")
gids <- unlist(mget(syms,org.Hs.egSYMBOL2EG))

##
e.lims <- quantile(as.vector(eff.1),c(0.0005,.9995),na.rm=T)
d.lims <- quantile(as.vector(dep.1),c(0.0005,.9995),na.rm=T)
e.lim1 <- range(eff.1[gids,],na.rm=T)
d.lim1 <- range(dep.1[gids,],na.rm=T)

ed.lims <- range(e.lims,d.lims,e.lim1,d.lim1)		

## Fig. 1A
setwd(plot_dir);setwd("01_genewise_relationship")
for(sym in syms){
	gid <- gids[sym]
	pdf(paste0(sym,".pdf"),width=5,height=5)
	par(mar=c(5,5,3,3))
	gene.ed.plot(gid=gid,ed.lims=ed.lims)
	Sys.sleep(1)
	dev.off()
}

## Fig. 2A
s <- svd(ed1,1,0)
f <- s$u %*% t(s$d) 

contrib <- s$d/sum(s$d) # c(0.659275, 0.340725)

setwd(plot_dir);setwd("01_genewise_relationship")
pdf(paste0("combine_scores.pdf"),width=5,height=5)
par(mar=c(5,5,3,3))
ratio <- c("100:0\n(CRISPR)","80:20","60:40","40:60","20:80","0:100\n(shRNA)")
gene.ed.plot(ed.lims=ed.lims)
m <- 3.5
cols <- brewer.pal(6,"Spectral")
for(i in 1:6){
	if(i==6){
		abline(v=0,col=cols[i],lwd=2)
		text(0,-m,ratio[i],srt=90)
	}else{
		tmp <- seq(0,1,.2)[i]
		abline(a=0,b=tmp/(1-tmp),col=cols[i],lwd=2)

		f <- sqrt(tmp^2+(1-tmp)^2)
		x <- -m*(1-tmp)/f
		y <- -m*(tmp)/f
		if(i==1){
			text(x,y,ratio[i],srt=atan(tmp/(1-tmp))/(pi/2)*90)
		}else{
			text(x,y,ratio[i],srt=atan(tmp/(1-tmp))/(pi/2)*90,pos=3)			
		}
	}
}

## PC1
abline(a=0,b=contrib[2]/contrib[1],lty=1,col="grey50")
tmp <- .35
m <- 4.5

f <- sqrt(tmp^2+(1-tmp)^2)
x <- -m*(1-tmp)/f
y <- -m*(tmp)/f
text(x,y,"PC1",srt=atan(tmp/(1-tmp))/(pi/2)*90,pos=3)

dev.off()


## Fig. 2B
syms <- c("CCND1","RAN","EIF5B","FOXD4")
gids <- unlist(mget(syms,org.Hs.egSYMBOL2EG))

setwd(rda_dir)
gid.scs <- lapply(1:6,function(k){
	c1 <- seq(0,1,.2)
	sc <- readRDS(file=paste0("r",c1[k],"/scores_15847g_423c.rds"))
	return(sc[gids,])
})

scs <- lapply(gids,function(gid){
	tmp.sc <- lapply(gid.scs,function(sc){
		sc[gid,]
	})
	names(tmp.sc) <- paste0("r",seq(0,1,.2))
	return(tmp.sc)
})

ns <- sapply(scs,function(x)sapply(x,length))
rs <- c("CRISPR","80:20","60:40","40:60","20:80","shRNA")
df <- data.frame(exp=unlist(scs),
	rs=factor(rep(rep(rs,ns[,1]),times=4),levels=rev(rs)),
	sym=factor(rep(syms,colSums(ns)),levels=syms))

setwd(plot_dir);setwd("01_genewise_relationship")
pdf("violin_4syms.pdf",width=7,height=3)
library(dplyr)
p <- df %>% ggplot(aes(fill=rs,x=sym,y=exp)) + 
	geom_violin(scale="width",width=.9) +
	scale_fill_brewer(palette="Spectral",direction=-1) +
	theme(legend.position="none") +
	theme_bw() + #scale_x_continuous(breaks = NULL) +
	xlab("") + ylab ("Dependency score (60:40)")
	 #geom_jitter(height=0,width=0.1,size=.1) +
print(p)
dev.off()

##　


