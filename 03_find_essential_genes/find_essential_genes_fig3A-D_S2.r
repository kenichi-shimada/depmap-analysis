library(RColorBrewer)
library(org.Hs.eg.db)
library(robustbase)
library(limma)
library(sp)

##
setwd(rda_dir)
ef <- readRDS(file="scores_15847g_423c_101119.rds")

## potency, selectivity

## note that the scores of the 10th, 20th and 42nd sensitive or resistant 
## cell lines are taken as 2.5th, 5th, and 10th percentile (sensitive) or 
## 97.5th, 95th, and 90th percentile (resistant).

q2.5 <- apply(ef,1,function(x)sort(x)[10])
q5 <- apply(ef,1,function(x)sort(x)[20])
q10 <- apply(ef,1,function(x)sort(x)[42])

q97.5 <- apply(ef,1,function(x)sort(x,decreasing=T)[10])
q95 <- apply(ef,1,function(x)sort(x,decreasing=T)[20])
q90 <- apply(ef,1,function(x)sort(x,decreasing=T)[42])

r95 <- q97.5-q2.5
r90 <- q95-q5
r80 <- q90-q10

qths <- list(q2.5,q5,q10) # for essential genes
qths2 <- list(q97.5,q95,q90) # for growth-suppressing genes

xrng <- range(unlist(qths)) + c(-1,1) * 0.01515 # range of histogram
xrng2 <- range(unlist(qths2)) + c(-1,1) * 0.01515 # range of histogram


## Identification of essential genes with three different criteria
## namely 2.5th, 5th, 10th percentiles

setwd(plot_dir);setwd("03_histogram_essential_genes")

lo.ngns <- list()
lo.idxs <- list()
lo.ths <- list()

for(i in 1:3){
	## fig.
	qth <- qths[[i]]
	thres <- c(2.5,5,10)
	nseq <- 100
	seq.br <- seq(xrng[1],xrng[2],length=nseq)

	x <- hist(qth,main=paste(thres[i]," percentile"),
		col="lightblue",border=NA,
		xlab="efficacy score",breaks=seq.br)

	nl <- length(x$breaks)
	br <- (x$breaks[-1]+x$breaks[-nl])/2
	ct <- x$counts
	names(ct) <- br
	mod <- as.numeric(names(which.max(ct)))

	pos <- (qth-mod)
	pos <- pos[pos > 0] ## 6275
	sd <- sqrt(sum(pos^2)/(length(pos)-1))

	ths <- mod * 2 - qnorm(1-c(1e-2,1e-3,1e-5),mean=mod,sd=sd)

	names(ths) <- sig.ps <- c("1e-2","1e-3","1e-5")
	lo.ths[[i]] <- ths
	lo.ngns[[i]] <- n.genes <- sapply(ths,function(n)sum(qth < n))

	lo.idxs[[i]] <- idx <- cbind(th2=qth < ths[1],
		th3=qth < ths[2],
		th5=qth < ths[3])

	ybr <- dnorm(br,mean=mod,sd=sd)
	f <- max(ct)/max(ybr)

	## Fig. 3A, S2
	filenames <- c("Fig3A_10th-lowest.pdf","FigS2A_20th-lowest.pdf","FigS2A_42th-lowest.pdf")
	pdf(filenames[i],width=5,height=5)
	par(mar=c(5,5,3,3))

	x <- hist(qth,main=paste(thres[i]," percentile"),
		col="lightblue",border=NA,
		xlab="efficacy score",breaks=seq.br)
	lines(seq.br,f*dnorm(seq.br,mean=mod,sd=sd),col=4,lwd=2)
	abline(v=ths,lty=2,col="grey70")
	text(ths,as.vector(par()$usr[3:4]%*%c(.25,.75)),
		paste("p<",sig.ps," (",n.genes,"genes)"),srt=90,pos=3,cex=.5)
	dev.off()

}

names(lo.ths) <- names(lo.idxs) <- names(lo.ngns) <- paste0("q",thres)
lo.idxs <- do.call(cbind,lo.idxs)
lo.idx.sum <- apply(lo.idxs,1,sum)
lo.ths <- round(do.call(cbind,lo.ths),3)
n.ess.gns <- sapply(lo.ngns,function(x){names(x) <- sig.ps;return(x)})

is.ess <- lo.idxs[,2] # q2.5, 1e-3
ess.genes <- names(which(is.ess))
ess.ef <- ef[is.ess,]

coef <- cor(t(ess.ef),method="spearman",use="pairwise.complete.obs")

##
setwd(rda_dir)
saveRDS(ess.genes,file="ess_2492_genes.rds")
saveRDS(is.ess,file="is.ess_2492_genes.rds")
save(lo.idxs,lo.ths,lo.idx.sum,n.ess.gns,file="eff-threshold.rda")
saveRDS(coef,file="eff_coef_2492_genes.rds")


## Identification of growth-suppressing genes
hi.ngns <- list()
hi.idxs <- list()
hi.ths <- list()
for(i in 1:3){
	qth <- qths2[[i]]
	thres <- c(97.5,95,90)
	nseq <- 100
	seq.br <- seq(xrng2[1],xrng2[2],length=nseq)

	x <- hist(qth,main=paste(thres[i]," percentile"),
		col="lightblue",border=NA,
		xlab="efficacy score",breaks=seq.br)

	nl <- length(x$breaks)
	br <- (x$breaks[-1]+x$breaks[-nl])/2
	ct <- x$counts
	names(ct) <- br
	mod <- as.numeric(names(which.max(ct)))
	upper <- 0.5
	abline(v=mod,col=2)
	pos <- (qth-mod)
	pos <- pos[pos > 0 & qth < upper] ## 6275
	sd <- sqrt(sum(pos^2)/(length(pos)-1))

	ths <- qnorm(1-c(1e-2,1e-3,1e-5),mean=mod,sd=sd)

	## what threhold to use? 
	names(ths) <- sig.ps <- c("1e-2","1e-3","1e-5")
	hi.ths[[i]] <- ths
	hi.ngns[[i]] <- n.genes <- sapply(ths,function(n)sum(qth > n))

	hi.idxs[[i]] <- idx <- cbind(th2=qth > ths[1],
		th3=qth > ths[2],
		th5=qth > ths[3])
}

names(hi.ths) <- names(hi.idxs) <- names(hi.ngns) <- paste0("q",thres)
hi.idxs <- do.call(cbind,hi.idxs)
hi.ths <- round(do.call(cbind,hi.ths),3)
hi.idx.sum <- apply(hi.idxs,1,sum)
n.gs.gns <- sapply(hi.ngns,function(x){names(x) <- sig.ps;return(x)})

is.gs <- hi.idxs[,2] ## q97.5, 1e-3
gs.genes <- names(which(is.gs))

##
setwd(rda_dir)
saveRDS(gs.genes,file="gs_65_genes.rds")
saveRDS(is.gs,file="is.gs_65_genes.rds")
save(hi.idxs,hi.ths,hi.idx.sum,n.gs.gns,file="gs-threshold.rda")


## Non-essential genes, ie genes that are neither essential nor growth-suppressing
lo.tab <- lo.idxs[,c(2,5,8)]
colnames(lo.tab) <- c("10th","20th","42th")
sum.lo.tab <- apply(lo.tab,1,sum)

hi.tab <- hi.idxs[,c(2,5,8)]
colnames(hi.tab) <- c("10th","20th","42th")
sum.hi.tab <- apply(hi.tab,1,sum)

is.non <- lo.idx.sum==0 & hi.idx.sum==0 ## 12283

setwd(rda_dir)
saveRDS(is.non,file="non-essential-genes.rds")


## fig. S2C
setwd(plot_dir);setwd("03_histogram_essential_genes")
pdf("figS2C_venn_ess_genes.pdf",width=5,height=5)
vennDiagram(lo.tab,main="Essential genes")
dev.off()


## compute 95% confidence interval  
df <- data.frame(q2.5=q2.5,q97.5=q97.5)[is.non,]
input <- data.frame(q2.5=q2.5)

##
fit <- lm(q97.5 ~ q2.5,data=df)
pred <- data.frame(predict(fit,input,interval="confidence"))
resid.non <- q97.5[is.non]-pred[is.non,1]

nseq <- 100
rng <- range(resid.non)
seq.br <- seq(rng[1],rng[2],length=nseq)

x <- hist(resid.non,breaks=nseq,col="lightpink",border=NA)
nl <- length(x$breaks)
br <- (x$breaks[-1]+x$breaks[-nl])/2
ct <- x$counts
names(ct) <- br
mod <- as.numeric(names(which.max(ct))) # mod
abline(v=mod,col=2)
neg <- (resid.non-mod)
neg <- neg[neg < 0] ## 6275
sd <- sqrt(sum(neg^2)/(length(neg)-1))
# abline(v=mod+c(1,-1)*sd,col=3)

ybr <- dnorm(br,mean=mod,sd=sd)
f <- max(ct)/max(ybr)

lines(seq.br,f*dnorm(seq.br,mean=mod,sd=sd),col=4,lwd=2)
d.ci95 <- qnorm(.975) * sd # 0.0632583



## fig. 3B
setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")

df <- data.frame(cbind(q2.5=q2.5,q97.5=q97.5)[is.non,])
fit <- lmrob(q97.5 ~ q2.5,data=df)

xs <- data.frame(q2.5=seq(-2.664009,-0.2337219,length=4))
ys <- predict(fit,xs)
inds <- c()

for(i in 1:4){
	x <- xs[[1]][i]
	y <- ys[i]
	dist.sq <- (q2.5-x)^2+(q97.5-y)^2
	ind <- which.min(dist.sq)
	inds[i] <- names(q2.5)[ind]
}

tmp.gns <- unlist(mget(inds,org.Hs.egSYMBOL))

cols <- brewer.pal(9,"YlOrRd")[9:6]

pdf("fig3B_smoothscatter_q2.5-q97.5.pdf",width=5,height=5)
par(mar = c(4,4,1,1))
smoothScatter(q2.5,q97.5,
	xlab="Perturbation score (the 10th most sensitive cell line)",
	ylab="Perturbation score (the 10th least sensitive cell line)",
	pch=20,cex=.5)
abline(a=0,b=1,lty=2)
abline(h=0,v=0,lty=2)

abline(coef(fit),lty=1,col=4)
abline(coef(fit)+c(d.ci95,0),lty=2,col=4)#"grey60")
abline(coef(fit)+c(-d.ci95,0),lty=2,col=4)#"grey60")
points(xs[[1]],ys,pch=20,col=cols)
text(xs[[1]],ys,tmp.gns,pos=4)
dev.off()

## fig. 3C
tmp <- ef[inds,]
rng.1 <- range(tmp)+c(-1,1)*.1
pdf("fig3C_histogram_efficacy_non-selective_4genes.pdf",width=4,height=5)
par(mfrow=c(4,1))
breaks <- seq(rng.1[1],rng.1[2],length=100)
for(i in 4:1){
	x <- tmp[i,]
	par(mar=c(3,5,0,1))
	hist(x,breaks=breaks,xlim=rng.1,main="",xlab="",ylim=c(0,80),border=NA,col=cols[i])
	text(-3.5,60,tmp.gns[i],cex=2,pos=4)
	abline(v=quantile(x,c(.025,.975)),lty=2)
}
dev.off()


##
input <- data.frame(q2.5=q2.5)
pred <- data.frame(predict(fit,input,interval="confidence"))

if(1){
	within <- q97.5 > pred$lwr-d.ci95 & q97.5 < pred$upr+d.ci95 # 12514/15847
}else{	
	within <- q97.5 > pred$lwr & q97.5 < pred$upr # 339
}

# points(q2.5[within],q97.5[within],pch=20) 
within.genes <- names(which(within))

##
rngs <- q97.5-q2.5
sds <- apply(ef,1,sd,na.rm=T)

df <- data.frame(q2.5=q2.5[within.genes],rng=rngs[within.genes])
fit <- lmrob(rng ~ q2.5, data=df)
pred.rngs <- predict(fit,newdata=data.frame(q2.5))

sel <- rngs/pred.rngs - 1


## fig. S2D
setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")
pdf("figS2D_q2.5_r95_only-within.pdf",width=5,height=5)
par(mar = c(4,4,1,1))
fit <- lmrob(rng ~ q2.5, data=df)
smoothScatter(q2.5[within],rngs[within],pch=20,
	xlab="efficacy score (2.5 percentile)",
	ylab="range (2.5-97.5 percentiles)")
abline(fit,col=4)
text(-1,0.8,"y = -0.2231 x + 0.2560\n(Adj. R2: 0.8197)",cex=1)
points(q2.5[inds],rngs[inds],pch=20,col=cols,cex=1.5)
dev.off()

## fig. S2E
pdf("figS2E_r95_sd_only-within.pdf",width=5,height=5)
par(mar = c(4,4,1,1))
smoothScatter(rngs[within],sds[within],pch=20,
	xlab="range (2.5-97.5 percentiles)",
	ylab="standard deviation")
fit <- lmrob(sds[within] ~ rngs[within])
abline(fit,col=4)
text(0.65,0.1,"Adj. R2: 0.9693",cex=1)
points(rngs[inds],sds[inds],pch=20,col=cols,cex=1.5)
dev.off()


## compute selectivity
nseq <- 200
xrng2 <- range(unlist(sel)) + c(-1,1) * 0.01

seq.br <- seq(xrng2[1],xrng2[2],length=nseq)

x <- hist(sel,main="",
	col="lightblue",border=NA,
	xlab="efficacy score",breaks=seq.br)

nl <- length(x$breaks)
br <- (x$breaks[-1]+x$breaks[-nl])/2
ct <- x$counts
names(ct) <- br
mod <- as.numeric(names(which.max(ct)))
abline(v=mod,col=2)

neg <- (sel-mod)
neg <- neg[neg < 0]
sd <- sqrt(sum(neg^2)/(length(neg)-1))

ybr <- dnorm(br,mean=mod,sd=sd)
f <- max(ct)/max(ybr)

lines(seq.br,f*dnorm(seq.br,mean=mod,sd=sd),col=4,lwd=2)

## what threhold to use? 
sel.ths <- mod * 2 - qnorm(c(1e-2,1e-3,1e-5),mean=mod,sd=sd)
names(sel.ths) <- sig.ps <- c("1e-2","1e-3","1e-5")

sel <- sel/sel.ths[2]

## Fig. 3D
xhist <- hist(q2.5,breaks=300,col="pink")
yhist <- hist(sel,breaks=300,col="pink")
top <- max(xhist$counts,yhist$counts)

setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")
pdf("Fig3D_efficacy_selectivity.pdf",width=5,height=5)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

par(mar = c(4,4,1,1))
smoothScatter(q2.5,sel,
	xlab="Efficacy (2.5 percentile)",
	ylab="Selectivity",
	pch=20,cex=.5)
if(0){
	xys1 <- locator(15,type="l",lty=2)
	xys2 <- locator(15,type="l",lty=2)
	library(sp)
	i1 <- which(point.in.polygon(q2.5,sel,xys1$x,xys1$y)==1)
	i2 <- which(point.in.polygon(q2.5,sel,xys2$x,xys2$y)==1)	
	i12 <- c(i1,i2)
	i12 <- i12[which(q2.5[i12] < 0)]
	setwd(rda_dir)
	save(xys1,xys2,i1,i2,file="polygons_efsel.rda")
	saveRDS(i12,file="index_selected_on_efsel_plot.rds")
}else{
	setwd(rda_dir)
	i12 <- readRDS("index_selected_on_efsel_plot.rds")
}

syms12 <- unlist(mget(names(q2.5)[i12],org.Hs.egSYMBOL))
text(q2.5[i12]-.05,sel[i12]-.2,syms12,pos=3,cex=.5)

abline(v=0,h=0,lty=1)

v.th <- lo.ths["1e-3","q2.5"]
h.th <- 1

abline(v=v.th,lty=2)
abline(h=h.th,lty=2)

par(mar = c(0,4,1,1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,col=c("grey","grey50")[(xhist$breaks< v.th)+1],border=NA,lwd=.25)

par(mar = c(4,0,1,1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0,col=c("grey","grey50")[(yhist$breaks > h.th)+1],border=NA,lwd=.25,horiz=TRUE)

dev.off()


## save objects
ef.sel <- data.frame(eff=q2.5,sel=sel)
thres <- list(eff=lo.ths["1e-3","q2.5"],sel=1)

setwd(rda_dir)
save(ef.sel,thres,file="eff-sel.rda")

eids.all <- rownames(ef.sel)
syms.all <- unlist(mget(eids.all,revmap(org.Hs.egSYMBOL2EG)))
syms.all <- paste0(syms.all,"(",eids.all,")")
ef.sel <- round(ef.sel,4)
rownames(ef.sel) <- syms.all
colnames(ef.sel) <- c("efficacy","selectivity")

write.csv(ef.sel,file="efficacy_selectivity_all-genes.csv")	

##
sum.ess.sel <- table(is.essential=is.ess,is.selselective=sel > 1)[2:1,2:1]
setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")
write.csv(sum.ess.sel,file="number_essential-selective-genes.csv")


