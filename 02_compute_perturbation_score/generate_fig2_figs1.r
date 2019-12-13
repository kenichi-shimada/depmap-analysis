library(ggplot2)
library(KernSmooth) ## bkde2D
library(RColorBrewer)
library(org.Hs.eg.db)
library(sp)
library(dplyr)


## load data
setwd(rda_dir)
x <- load(file="eff_dep_19q3.rda") ## eff.1, dep.1, both 15847 x 423 cells
sc <- readRDS(file="scores_15847g_423c_101119.rds")

##
e.lims <- quantile(as.vector(eff.1),c(0.0005,.9995),na.rm=T)
d.lims <- quantile(as.vector(dep.1),c(0.0005,.9995),na.rm=T)
ed.lims <- range(e.lims,d.lims)


## Fig. 2A: relationship between cripsr and shrna
gs <- 500
bws <- c(diff(e.lims),diff(d.lims))/50
est <- bkde2D(ed1,bandwidth=bws,gridsize=rep(gs,2))
ncols <- rep(rep(c("grey80","grey50"),c(4,1)),length.out=200)[-(1:4)]
levs <- c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1:4)

uniq.cols <- c("grey80","grey70","grey50","grey20")#brewer.pal(8,"YlOrRd")[c(2,4,6,8)]
cols <- uniq.cols[c(1,2,2,2,3,3,3,4,4,4,4)]

setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/genewise_relationship")
pdf("CRIPSR_shRNA_relationship.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
ed.lims <- c(-3,2)
smoothScatter(as.vector(eff.1),as.vector(dep.1),
	transformation=function(x)x^.5,colramp=colorRampPalette(c("white",blues9[3:6])),
	xlab="CRISPR efficacy",ylab="shRNA efficacy",asp=1,xlim=ed.lims,ylim=ed.lims,
	pch=NA)
contour(est$x1,est$x2,est$fhat,levels=levs,xlim=e.lims,ylim=d.lims,col=cols,add=T)
abline(a=0,b=1,lty=2)
abline(h=0,v=0,lty=2)

a <- list(x=c(-2,-0.25,-0.25,-2,-2),y=c(0.2,0.2,-0.2,-0.2,0.2))
b <- list(x=c(0.2,0.2,-0.2,-0.2,0.2),y=c(-2,-0.25,-0.25,-2,-2))
polygon(a,border="red")
polygon(b,border="red")
text(-2,0.27,"A",col=2,cex=1)
text(-0.27,-2,"B",col=2,cex=1)
dev.off()

gns <- rownames(eff.1)
all.a <- table(point.in.polygon(as.vector(eff.1),as.vector(dep.1),a$x,a$y)==1)[2:1]
all.b <- table(point.in.polygon(as.vector(eff.1),as.vector(dep.1),b$x,b$y)==1)[2:1]
pvals <- sapply(gns,function(g){
	in.a <- table(point.in.polygon(eff.1[g,],dep.1[g,],a$x,a$y)==1)[2:1]
	in.b <- table(point.in.polygon(eff.1[g,],dep.1[g,],b$x,b$y)==1)[2:1]
	a1 <- all.a - in.a
	b1 <- all.b - in.b
	mat.a <- rbind(in.a,a1)
	mat.b <- rbind(in.b,b1)
	mat.a[is.na(mat.a)] <- 0
	mat.b[is.na(mat.b)] <- 0
	p.a <- fisher.test(mat.a,alternative="greater")$p.value
	p.b <- fisher.test(mat.b,alternative="greater")$p.value	
	return(c(p.a,p.b))
})
nlps <- round(-log10(t(pvals)),3) ## four pvals are 0
nlps[nlps > 315] <- 315



## Fig. S1A
setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/genewise_relationship")
pdf("eCDF_boxes.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
plot.ecdf(nlps[,1],verticals=T,pch=NA,axes=F,col=NA,
	xlab="P-value (Fisher's test)",ylab="eCDF(genes)",main="")
box()
axis(2)
labelsX=c("1",parse(text=paste0("10^-",seq(50,300,50))))
axis(1,at=seq(0,300,50),labels=labelsX,cex.axis=.8)
plot.ecdf(nlps[,1],verticals=T,axes=F,col=1,add=T,pch=NA)
plot.ecdf(nlps[,2],verticals=T,col=2,add=T,pch=NA)
legend(150,0.4,c("box A","box B"),lty=1,col=1:2)
dev.off()

## Fig. S1B
pdf("CRISPR_shRNA_relationship_polygon.pdf",width=5,height=5)
smoothScatter(as.vector(eff.1),as.vector(dep.1),transformation=function(x)x^.5,colramp=colorRampPalette(c("white",blues9[3:6])),
	xlab="CRISPR efficacy",ylab="shRNA efficacy",asp=1,xlim=ed.lims,ylim=ed.lims,pch=NA)
contour(est$x1,est$x2,est$fhat,levels=levs,xlim=e.lims,ylim=d.lims,col=cols,add=T)
abline(a=0,b=1,lty=2)
abline(h=0,v=0,lty=2)
dev.off()



## Fig. 2A
setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/genewise_relationship")

gene.ed.plot <- function(gid="2879",sym="GPX4",ed.lims=ed.lims){
	uniq.cols <- c("grey80","grey70","grey50","grey20")#brewer.pal(8,"YlOrRd")[c(2,4,6,8)]
	cols <- uniq.cols[c(1,2,2,2,3,3,3,4,4,4,4)]
	par(mar=c(5,5,3,3))
	smoothScatter(as.vector(eff.1),as.vector(dep.1),
		transformation=function(x)x^.5,
		colramp=colorRampPalette(c("white",blues9[3:6])),
		xlab="",ylab="",asp=1,
		main=sym,xlim=ed.lims,ylim=ed.lims,pch=NA)
	mtext("CRISPR efficacy",1,3)
	mtext("shRNA efficacy",2,3)
	contour(est$x1,est$x2,est$fhat,levels=levs,add=T,col=cols)
	abline(h=0,v=0,lty=2)
	abline(a=0,b=s$d[2]/s$d[1],col=2,lty=1)
	abline(a=0,b=1,lty=2)
	points(eff.1[gid,],dep.1[gid,],pch=20,col=1,cex=.5)
}

sym <- c("CCND1","RAN","SNAPC4","ACTG2")
gids <- c("595","5901","6621","72")

pdf("ccnd1.pdf",width=5,height=5)
gene.ed.plot(gid="595",sym="CCND1") # mget("CCND1",org.Hs.egSYMBOL2EG)
dev.off()

pdf("ran.pdf",width=5,height=5)
gene.ed.plot(gid="5901",sym="RAN") # mget("RAN",org.Hs.egSYMBOL2EG)
dev.off()

pdf("snapc4.pdf",width=5,height=5)
gene.ed.plot(gid="6621",sym="SNAPC4") # mget("GPX4",org.Hs.egSYMBOL2EG)
dev.off()

pdf("actg2.pdf",width=5,height=5)
gene.ed.plot(gid="72",sym="ACTG2") # mget("GPX4",org.Hs.egSYMBOL2EG)
dev.off()



## Fig. 2F
pdf("utp14a.pdf",width=5,height=5)
gene.ed.plot(gid="10813",sym="UTP14A",ed.lims=ed.lims)
dev.off()



## Fig. 2B
sym <- c("CCND1","RAN","SNAPC4","ACTG2")
gids <- c("595","5901","6621","72")

sc1 <- t(sc[gids,])
df <- data.frame(exp=as.vector(sc1),sym=factor(rep(sym,each=423),levels=sym))

setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/genewise_relationship")
pdf("violin_score_4genes.pdf",width=4,height=3)
p <- ggplot(df,aes(sym,exp)) + 
	geom_violin(aes(fill=sym)) + #scale_fill_brewer(palette="Spectral",direction=-1) +
	theme_bw() + geom_jitter(height=0,width=0.1,size=.1) +
	theme(legend.position="none")
p
dev.off()

## below here is comparisong between CRISPR, shRNA and score (perturbation score).
## corresponding to Fig. 2 and Fig. S1.

mid.e <- apply(eff.1,1,median,na.rm=T)
mid.d <- apply(dep.1,1,median,na.rm=T)
mid.sc <- apply(sc,1,median,na.rm=T)

##
lm1  <- lm(mid.sc ~ mid.e)
pred.sc <- predict(lm1)
idx <- (mid.sc - pred.sc < -0.1)+1

dns.gns <- unlist(mget(names(which(idx==2)),org.Hs.egSYMBOL))
all.gns <- unlist(mget(rownames(sc),org.Hs.egSYMBOL))
dns.prefix <- sub("^(...).+","\\1",dns.gns)
all.prefix <- sub("^(...).+","\\1",all.gns)

## Three letter synonym is too short to distinguish some gene family.
## Made some of the gene names longer.

dns.prefix[dns.gns=="COP1"] <- "COP1"
dns.prefix[grep("^COPS",dns.gns)] <- "COPS"
dns.prefix[grep("^POLR1",dns.gns)] <- "POLR1"
dns.prefix[grep("^POLR2",dns.gns)] <- "POLR2"
dns.prefix[grep("^POLR3",dns.gns)] <- "POLR3"
dns.prefix[grep("^PRPF",dns.gns)] <- "PRPF"
dns.prefix[grep("^PRPH",dns.gns)] <- "PRPH"
dns.prefix[grep("^PRPS",dns.gns)] <- "PRPS"
dns.prefix[grep("^SNRK",dns.gns)] <- "SNRK"

all.prefix[all.gns=="COP1"] <- "COP1"
all.prefix[grep("^COPS",all.gns)] <- "COPS"
all.prefix[grep("^POLR1",all.gns)] <- "POLR1"
all.prefix[grep("^POLR2",all.gns)] <- "POLR2"
all.prefix[grep("^POLR3",all.gns)] <- "POLR3"
all.prefix[grep("^PRPF",all.gns)] <- "PRPF"
all.prefix[grep("^PRPH",all.gns)] <- "PRPH"
all.prefix[grep("^PRPS",all.gns)] <- "PRPS"
all.prefix[grep("^SNRK",all.gns)] <- "SNRK"

##
tab.dns <- sort(table(dns.prefix),decreasing=T)
tab.all <- sort(table(all.prefix),decreasing=T)

##
dn.prefix <- names(tab.dns[tab.dns >= 5])
pref2gns <- tapply(all.gns,all.prefix,identity)[dn.prefix]
pref2nams <- tapply(rownames(sc),all.prefix,function(x){
	unlist(mget(x,org.Hs.egGENENAME))
})[dn.prefix]

nlps <- sapply(dn.prefix,function(x){
	has.pref <- all.prefix %in% x
	is.dn <- all.gns %in% dns.gns
	mat <- table(prefix=has.pref,neg=is.dn)[2:1,2:1]
	nlp <- -log10(fisher.test(mat,alternative="greater")$p)
	return(c(mat[1,],nlp))
})

df.dns <- data.frame(sym=dn.prefix,hit=nlps[1,],non.hit=nlps[2,],
	nlp=round(nlps[3,],1),stringsAsFactors=F)
df.dns.1 <- df.dns %>% filter(nlp >5) %>% arrange(desc(nlp))
nm <- df.dns.1$sym

pref2gns <- tapply(all.gns,all.prefix,identity)[nm]
pref2nams <- tapply(rownames(sc),all.prefix,function(x){
	unlist(mget(x,org.Hs.egGENENAME))
})[nm]

##
desc <- c(RPS="small ribosomal protein",
	PSM="proteasome",
	RPL="large ribosomal protein",
	EIF="eukaryotic translation initiation",
	PRPF="pre-mRNA processing factor",
	SNR="small nuclear ribonucleoprotein",
	COPS="COP9 signalosome",
	SF3="splicing factor 3a/b",
	POLR2="RNA polymerase II",
	NUP="nucleoprotein",
	CCT="chaperonin containing TCP1",
	DDX="DEAD-box helicase")
df.dns.1 <- df.dns.1 %>% mutate(desc=desc)

if(0){
	## Fig. 2D
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")
	write.csv(df.dns.1,file="shrna-sensitive.csv")

	## Fig. 2C
	lm1  <- lm(mid.sc ~ mid.e)
	pred.sc <- predict(lm1)
	idx <- (mid.sc - pred.sc < -0.1)+1
	dc <- densCols(mid.e,mid.sc,
		colramp = colorRampPalette(rev(grey.colors(9))[-(1:3)]))
	dc[idx==2] <- 2

	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	png("scatter_crispr_cmpbd.png",width=500,height=500)
	par(mar=c(5,5,3,3))
	plot(mid.e,mid.sc,pch=20,cex=.7,col=dc,
		# xlab="median CRISPR score",
		# ylab="median Combined score")
		xlab="",ylab="",cex.axis=1.5)
		abline(h=0,v=0,col=1,lty=2)
		abline(coef(lm1),col=1)
		abline(coef(lm1)+c(-0.1,0),col=1,lty=2)
	dev.off()

	## Fig. S1C
	lm2  <- lm(mid.sc ~ mid.d)
	pred.sc <- predict(lm2)
	idx <- (mid.sc - pred.sc < -0.1)+1
	dc <- densCols(mid.e,mid.sc,
		colramp = colorRampPalette(rev(grey.colors(9))[-(1:3)]))
	dc[idx==2] <- 2

	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	png("scatter_shrna_cmpbd.png",width=500,height=500)
	par(mar=c(5,5,3,3))
	plot(mid.d,mid.sc,pch=20,cex=.7,col=dc,
		xlab="",ylab="",cex.axis=1.5)
	abline(h=0,v=0,col=1,lty=2)
	abline(coef(lm2),col=1)
	abline(coef(lm2)+c(-0.1,0),col=1,lty=2)
	dev.off()
}

## correlation - median value
mids <- cbind(d=mid.d,e=mid.e,sc=mid.sc)
cors <- cor(mids,method="spearman")
if(0){
	## Fig. S1D
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	write.csv(cors,file="corelation_mid_scores.csv")
}

if(0){
	## Fig. S1C
	lm1  <- lm(mid.sc ~ mid.e)
	pred.sc <- predict(lm1)
	idx <- (mid.sc - pred.sc < -0.1)+1
	dc <- densCols(mid.e,mid.d,
		colramp = colorRampPalette(rev(grey.colors(9))[-(1:3)]))
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	png("scatter_crispr_shrna.png",width=500,height=500)
	par(mar=c(5,5,3,3))
	plot(mid.e[idx==1],mid.d[idx==1],col=dc[idx==1],asp=1,cex=.7,pch=20,
		# xlab="median CRISPR score",
		# ylab="median Combined score")
		xlab="",ylab="",cex.axis=1.5)
		abline(h=0,v=0,col=1,lty=2)
	points(mid.e[idx==2],mid.d[idx==2],col=2,asp=1,cex=.7,pch=20)
		abline(h=0,v=0,col=1,lty=2)
	dev.off()

}

## assessing the differences in correlation
setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")
ess <- readRDS("ess_2492_genes.rds")
ess.syms <- unlist(mget(ess,org.Hs.egSYMBOL))

if(0){
	e1 <- cor(t(eff.1[ess,]),method="spearman",use="pairwise.complete.obs")
	d1 <- cor(t(dep.1[ess,]),method="spearman",use="pairwise.complete.obs")
	sc1 <- cor(t(sc[ess,]),method="spearman",use="pairwise.complete.obs")
	diag(e1) <- diag(d1) <- diag(sc1) <- 0
	rownames(e1) <- rownames(d1) <- rownames(sc1) <- ess.syms
	colnames(e1) <- colnames(d1) <- colnames(sc1) <- ess.syms

	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")
	save(e1,d1,sc1,file="correlation_3ways.rda")
}else{
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")
	load(file="correlation_3ways.rda") # e1,d1,sc1
}

eff.e <- apply(eff.1,1,function(x)quantile(x,.025,na.rm=T))
eff.d <- apply(dep.1,1,function(x)quantile(x,.025,na.rm=T))
eff.sc <- apply(sc,1,function(x)quantile(x,.025,na.rm=T))

sel.e <- apply(eff.1,1,function(x)diff(quantile(x,c(.025,.975),na.rm=T)))
sel.d <- apply(dep.1,1,function(x)diff(quantile(x,c(.025,.975),na.rm=T)))
sel.sc <- apply(sc,1,function(x)diff(quantile(x,c(.025,.975),na.rm=T)))

##
cor.i <- lower.tri(e1)
cor.df <- data.frame(CRISPR=e1[cor.i],shRNA=d1[cor.i],combined=sc1[cor.i])
ind <- sample(nrow(cor.df),10^5)
rng <- range(unlist(cor.df))

if(0){
	## Fig. S1E
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	pdf("coef_comparison.pdf",width=7,height=7)
	par(oma=c(2,2,0,0))
	par(mfrow=c(3,3))
	for(i in 1:3){
		for(j in 1:3){
			par(mar=c(1,1,1,1))
			if(i==j){
				plot(cor.df[c(j,i)],type="n",xlim=rng,ylim=rng,axes=F)
				box()
				text(mean(par()$usr[1:2]),mean(par()$usr[3:4]),c("CRISPR","shRNA","combined")[i],cex=3)	
			}else if(i < j){
				plot.new()
			}else{
				smoothScatter(cor.df[c(j,i)],asp=1,
					colramp = colorRampPalette(c("white", blues9[-c(1:3)])),axes=F,
					xlim=rng,ylim=rng)
				box()
				abline(h=0,v=0)
				if(i==3 & j != 3)axis(1)
				if(j==1 & i != 1)axis(2)
			}
		}
	}
	dev.off()
}

th <- .4
sig.i <- which(cor.df$c - cor.df$C > th)
dc[sig.i] <- 2

if(0){
	## Fig. 2E (plot)
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	pdf("coef_crispr_cmbd.pdf",width=5,height=5)
	par(mar=c(5,5,3,3))
	smoothScatter(cor.df[c("CRISPR","combined")],asp=1,
		colramp = colorRampPalette(c("white", blues9[-c(1:4)])),axes=F,
		xlim=rng,ylim=rng,pch=20,cex=.3,
		main="Pairwise correlation coefficients",
		xlab="CRISPR",
		ylab="Perturbation - combined")
	points(cor.df[sig.i,c("CRISPR","combined")],col=2,pch=20,cex=.5)
	box()
	abline(h=0,v=0)
	axis(1)
	axis(2)
	abline(a=0,b=1,lty=2)
	abline(a=th,b=1,lty=2)
	dev.off()
}

##
di.all <- data.frame(row=unlist(sapply(2:2492,function(x)(x:2492))),col=rep(1:2491,2491:1))
di <- di.all[sig.i,]
sig.syms <- cbind(ess.syms[di$r],ess.syms[di$c])
sort.sig.syms <- data.frame(sort(table(as.vector(sig.syms)),decreasing=T))

if(0){
	## Fig. 2E (list)
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	write.csv(sort.sig.syms,file="sorted.sig.syms.csv")
}

##
utps <- grep("^UTP14",ess.syms)
utps.id <- which(di.all$row%in% utps|di.all$col %in% utps)

eutp <- e1["UTP14A",ess.syms]
sutp <- sc1["UTP14A",ess.syms]

setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")
clust <- readRDS("hierarchy_cluster.rds")
this.mem <- clust$m1[clust$sym=="UTP14A"]
this.sym <- clust$sym[clust$m1==this.mem]
this.sym <- this.sym[this.sym!="UTP14A"]

if(0){
	## Fig. 2G
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/compare_eff_sel_score")
	pdf("coef_crispr_cmbd_utp14a.pdf",width=5,height=5)
	par(mar=c(5,5,3,3))
	dc <- densCols(eutp,sutp,
		colramp = colorRampPalette(blues9[-(1:3)]))
	plot(eutp,sutp,pch=20,col=dc,asp=1,cex=.7,
		xlab="correlation coefficients (CRISPR)",
		ylab="correlation coefficients (combined)")
	abline(h=0,v=0)
	points(eutp[this.sym],sutp[this.sym],pch=20,col=2)
	text(eutp[this.sym],sutp[this.sym],this.sym,pch=20,cex=.7,col=2,pos=1)
	dev.off()
}
