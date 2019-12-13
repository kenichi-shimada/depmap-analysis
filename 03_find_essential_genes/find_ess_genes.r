library(RColorBrewer)
library(org.Hs.eg.db)
library(robustbase)
library(limma)
library(sp)

##
setwd(rda_dir)
ef <- readRDS(file="scores_15847g_423c_101119.rds")

## potency, selectivity
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

## Identification of essential and growth-suppressing genes
setwd(plot_dir);setwd("histo_qantiles")
if(0){
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

		if(0){
			## Fig. 3A, S2
			pdf(paste0("lower_",i,".pdf"),width=5,height=5)
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

	if(0){
		setwd(rda_dir)
		saveRDS(ess.genes,file="ess_2492_genes.rds")
		saveRDS(is.ess,file="is.ess_2492_genes.rds")
		save(lo.idxs,lo.ths,lo.idx.sum,n.ess.gns,file="eff-threshold.rda")
		saveRDS(coef,file="eff_coef_2492_genes.rds")
	}

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
	n.ts.gns <- sapply(hi.ngns,function(x){names(x) <- sig.ps;return(x)})

	is.ts <- hi.idxs[,2] ## q97.5, 1e-3
	ts.genes <- names(which(is.ts))

	if(0){
		setwd(rda_dir)
		saveRDS(ts.genes,file="ts_65_genes.rds")
		saveRDS(is.ts,file="is.ts_65_genes.rds")
		save(hi.idxs,hi.ths,hi.idx.sum,n.ts.gns,file="ts-threshold.rda")
	}
}else{
	setwd(rda_dir)
	ess.genes <- readRDS(file="ess_2492_genes.rds")
	is.ess <- readRDS(file="is.ess_2492_genes.rds")
	load(file="eff-threshold.rda") # lo.idxs,lo.ths,lo.idx.sum,n.ess.gns

	ts.genes <- readRDS(file="ts_65_genes.rds")
	is.ts <- readRDS(file="is.ts_65_genes.rds")
	load(file="ts-threshold.rda") # hi.idxs,hi.ths,hi.idx.sum,n.ts.gns,

	coef <- readRDS(file="eff_coef_2492_genes.rds")

}

if(0){
	## genes that are neither essential nor growth-suppressing
	lo.tab <- lo.idxs[,c(2,5,8)]
	colnames(lo.tab) <- c("q97.5","q95","q90")
	sum.lo.tab <- apply(lo.tab,1,sum)

	hi.tab <- hi.idxs[,c(2,5,8)]
	colnames(hi.tab) <- c("q97.5","q95","q90")
	sum.hi.tab <- apply(hi.tab,1,sum)

	is.non <- lo.idx.sum==0 & hi.idx.sum==0 ## 12283

	setwd(rda_dir)
	saveRDS(is.non,file="non-essential-genes.rds")
}else{
	setwd(rda_dir)
	is.non <- readRDS(file="non-essential-genes.rds")
}

if(0){
	## fig. S2C
	setwd(plot_dir);setwd("histo_qantiles")
	pdf("venn_ess_genes.pdf",width=5,height=5)
	vennDiagram(lo.tab,main="Essential genes")
	dev.off()
}

## q10 vs r80 - non significant
i <- 1
nq <- c(2.5,5,10)[i]
nq2 <- c(97.5,95,90)[i]
nr <- c(95,90,80)[i]
qth <- get(paste0("q",nq))
qth2 <- get(paste0("q",nq2))
rth <- get(paste0("r",nr))

## compute ci95 for qth2 (independent of qth)
df <- data.frame(qth=qth,qth2=qth2)[is.non,]
qth1 <- data.frame(qth=qth)

if(1){
	fit <- lm(qth2 ~ qth,data=df)
	pred <- data.frame(predict(fit,qth1,interval="confidence"))
	resid.non <- qth2[is.non]-pred[is.non,1]

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
}

## fig. S3a-b
setwd(plot_dir);setwd("ess_scatter")

df <- data.frame(cbind(qth=qth,qth2=qth2)[is.non,])
fit <- lmrob(qth2 ~ qth,data=df)

xs <- data.frame(qth=seq(-2.664009,-0.2337219,length=4))
ys <- predict(fit,xs)
inds <- c()

for(i in 1:4){
	x <- xs[[1]][i]
	y <- ys[i]
	dist.sq <- (qth-x)^2+(qth2-y)^2
	ind <- which.min(dist.sq)
	inds[i] <- names(qth)[ind]
}

tmp.gns <- unlist(mget(inds,org.Hs.egSYMBOL))

cols <- brewer.pal(9,"YlOrRd")[9:6]

if(0){
	## fig. 3B
	setwd(plot_dir);setwd("ess_scatter")
	pdf("smoothscatter_q2.5-q97.5.pdf",width=5,height=5)
	par(mar = c(4,4,1,1))
	smoothScatter(qth,qth2,xlab="efficacy score (2.5 percentile)",ylab=" efficacy score (97.5 percentile)",
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
	pdf("histogram_efficacy_non-selective_4genes.pdf",width=4,height=5)
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
}

##
qth <- q2.5
qth2 <- q97.5

qth1 <- data.frame(qth=qth)
pred <- data.frame(predict(fit,qth1,interval="confidence"))

if(1){
	within <- qth2 > pred$lwr-d.ci95 & qth2 < pred$upr+d.ci95 # 12514/15847
}else{	
	within <- qth2 > pred$lwr & qth2 < pred$upr # 339
}

# points(qth[within],qth2[within],pch=20) 
within.genes <- names(which(within))

##
rngs <- q97.5-q2.5
sds <- apply(ef,1,sd,na.rm=T)

df <- data.frame(qth=qth[within.genes],rng=rngs[within.genes])
fit <- lmrob(rng ~ qth, data=df)
pred.rngs <- predict(fit,newdata=data.frame(qth))

sel <- rngs/pred.rngs - 1

##
if(0){
	par(mar = c(4,4,1,1))
	smoothScatter(qth,qth2,xlab=nq,ylab=nq2)
	abline(a=0,b=1,lty=2)
	abline(h=0,v=0,lty=2)
	points(qth[within.genes],qth2[within.genes],pch=20,cex=.5)
}

##
if(0){
	## fig. s1
	setwd(plot_dir);setwd("ess_scatter")
	pdf("q2.5_r95_only-within.pdf",width=5,height=5)
	par(mar = c(4,4,1,1))
	fit <- lmrob(rng ~ qth, data=df)
	smoothScatter(qth[within],rngs[within],pch=20,
		xlab="efficacy score (2.5 percentile)",
		ylab="range (2.5-97.5 percentiles)")
	abline(fit,col=4)
	text(-1,0.8,"y = -0.2231 x + 0.2560\n(Adj. R2: 0.8197)",cex=1)
	points(qth[inds],rngs[inds],pch=20,col=cols,cex=1.5)
	dev.off()

	## fig. s3
	pdf("r95_sd_only-within.pdf",width=5,height=5)
	par(mar = c(4,4,1,1))
	smoothScatter(rngs[within],sds[within],pch=20,
		xlab="range (2.5-97.5 percentiles)",
		ylab="standard deviation")
	fit <- lmrob(sds[within] ~ rngs[within])
	abline(fit,col=4)
	text(0.65,0.1,"Adj. R2: 0.9693",cex=1)
	points(rngs[inds],sds[inds],pch=20,col=cols,cex=1.5)
	dev.off()
}

##
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
# ths <- c(th.2,th.3,th.5)
sel.ths <- mod * 2 - qnorm(c(1e-2,1e-3,1e-5),mean=mod,sd=sd)
names(sel.ths) <- sig.ps <- c("1e-2","1e-3","1e-5")

sel <- sel/sel.ths[2]

##
if(0){
	xhist <- hist(qth,breaks=300,col="pink")
	yhist <- hist(sel,breaks=300,col="pink")
	top <- max(xhist$counts,yhist$counts)

	setwd(plot_dir);setwd("ess_scatter")
	pdf("efficacy_selectivity.pdf",width=5,height=5)
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

	par(mar = c(4,4,1,1))
	smoothScatter(qth,sel,
		xlab="Efficacy (2.5 percentile)",
		ylab="Selectivity",
		pch=20,cex=.5)
	if(0){
		xys1 <- locator(15,type="l",lty=2)
		xys2 <- locator(15,type="l",lty=2)
		library(sp)
		i1 <- which(point.in.polygon(qth,sel,xys1$x,xys1$y)==1)
		i2 <- which(point.in.polygon(qth,sel,xys2$x,xys2$y)==1)	
		i12 <- c(i1,i2)
		i12 <- i12[which(qth[i12] < 0)]
		setwd(rda_dir)
		save(xys1,xys2,i1,i2,file="polygons_efsel.rda")
		saveRDS(i12,file="index_selected_on_efsel_plot.rds")
	}else{
		setwd(rda_dir)
		i12 <- readRDS("index_selected_on_efsel_plot.rds")
	}

	syms12 <- unlist(mget(names(qth)[i1],org.Hs.egSYMBOL))
	text(qth[i12]-.05,sel[i12]-.2,syms1,pos=3,cex=.5)

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
}

f(0){
	ef.sel <- data.frame(eff=qth,sel=sel)
	thres <- list(eff=lo.ths["1e-3","q2.5"],sel=1)

	setwd(rda_dir)
	save(ef.sel,thres,file="eff-sel.rda")

	sum.ess.sel <- table(eff=is.ess,sel=sel > 1)[2:1,2:1]
	setwd(plot_dir);setwd("ess_scatter")
	write.csv(sum.ess.sel,file="n_ess-sel.csv")
}

