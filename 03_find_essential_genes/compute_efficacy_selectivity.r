library(L1pack)
library(org.Hs.eg.db)

## ef.sel
if(0){
	fits <- parallel::mclapply(1:6,function(k){
		c1 <- seq(0,1,.2)[k]

		##
		setwd(rda_dir);setwd(paste0("r",c1))
		ef <- readRDS(file=paste0("scores_15847g_423c.rds"))
		##
		setwd(plot_dir);setwd("ess_scatter")
		pdf(paste0("qth_qth2_smooth_",k,".pdf"),width=8,height=6)

		fit <- ef.sel <- list()

		qtl <- c(0.5,1,2.5,5,10,25)
		qtl2 <- 100-qtl

		par(mfrow=c(2,3))
		for(i in 1:6){
			cat("*")
			qth <- apply(ef,1,function(x)quantile(x,qtl[i]/100,na.rm=T))
			qth2 <- apply(ef,1,function(x)quantile(x,qtl2[i]/100,na.rm=T))
			rth <- qth2-qth

			lf1 <- l1fit(qth,qth2)
			smoothScatter(qth,qth2,pch=20,
				main=paste0("mix: ",c1,", th: ",qtl[i],"%"))
			abline(lf1,col=2)
			abline(h=0,v=0)

			fit[[i]]  <- l1fit(qth,rth)
			# smoothScatter(fit[[i]]$fitted.values,rth,pch=20)
			# abline(a=0,b=1)
			# abline(h=0,v=0)

			sel <- rth/fit[[i]]$fitted.values-1

			ef.sel[[i]] <- data.frame(eff=qth,sel=sel)
		}
		dev.off()

		coefs <- lapply(fit,function(x)x$coefficients)
		names(coefs) <- names(ef.sel) <- paste0("q",qtl)
		return(list(coef=coefs,ef.sel=ef.sel))
	},mc.cores=6)

	coefs <- lapply(fits,function(x)x$coef)
	ef.sels <- lapply(fits,function(x)x$ef.sel)

	names(coefs) <- names(ef.sels) <- paste0("r",seq(0,1,.2))

	setwd(rda_dir)
	save(coefs,ef.sels,file="eff-sel_6thres.rda")
}else{
	setwd(rda_dir)
	load("eff-sel_6thres.rda")
}

## choose ef.sel when score consists of CRISPR:shRNA = 60:40
ef.sel <- ef.sels$r0.6$q1

##
setwd(rda_dir)
sc <- readRDS(file=paste0("r0.6/scores_15847g_423c.rds")) # score matrix
load(file="find_essential_genes.rda") ## essential genes

##
qth <- apply(sc,1,function(x)quantile(x,0.01,na.rm=T))
qth2 <- apply(sc,1,function(x)quantile(x,0.99,na.rm=T))

fit <- l1fit(qth,qth2)

min.qth <- min(qth)
x <- hist(qth,breaks=100)
mid.qth <- x$mids[which.max(x$density)]
xs <- data.frame(qth=seq(min.qth,mid.qth,length=4))
ys <- cbind(1,xs$qth) %*% as.matrix(fit$coeff)

gids.1 <- sapply(1:4,function(j){
	x <- xs[[1]][j]
	y <- ys[j]
	dist.sq <- (qth-x)^2+(qth2-y)^2
	ind <- which.min(dist.sq)
	return(names(qth)[ind])
})
syms.1 <- unlist(mget(gids.1,org.Hs.egSYMBOL))

syms <- c("PSMB5","CCND1","CTNNB1","TSC2",syms.1)
gids <- c(unlist(mget(syms,org.Hs.egSYMBOL2EG)))

sc1 <- t(sc[gids,])
df <- data.frame(exp=as.vector(sc1),
	sym=factor(rep(syms,each=423),levels=syms),
	class=factor(rep(c("selective","non-selective"),each=423*4),
		levels=c("selective","non-selective")))

qth.1 <- qth[gids]
qth2.1 <- qth2[gids]


## Fig. 3A
setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")
pdf("qth_qth2_p1.pdf",width=5,height=5)
# png("qth_qth2_v2.png",width=800,height=800)
par(mar=c(5,5,3,3))
cols <- densCols(data.frame(qth,qth2),colramp = colorRampPalette(blues9[-(1:3)]))
plot(qth,qth2,pch=20,col=cols,type="n",
	xlab="1st percentile score",ylab="99th percentile score")
	#cex.lab=2,cex.axis=2)
abline(h=0,v=0,lty=1,col="grey80")
abline(a=0,b=1,lty=2)
points(qth,qth2,pch=20,col=cols,cex=.5)
abline(fit$coefficients,col=2)
cols <- brewer.pal(11,"Spectral")[c(4,7)]
points(qth.1,qth2.1,bg=rep(cols,each=4),col=1,cex=1,pch=21)
dev.off()

##
pdf("qth_qth2_p2.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
cols <- densCols(data.frame(qth,qth2),colramp = colorRampPalette(rev(blues9[-1])))
plot(qth,qth2,pch=20,col=cols,type="n",
	xlab="",ylab="",
	axes=F)
text(qth.1[1:4]-.25,qth2.1[1:4]+.2,syms[1:4])
text(qth.1[5:8]+.3,qth2.1[5:8]-c(0,.3,.2,.2),syms[5:8])
dev.off()


## Fig. 3B
setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")
pdf("violin_8syms.pdf",width=6,height=3)
cols <- brewer.pal(11,"Spectral")[c(4,7)]
p <- df %>% ggplot(aes(x=sym,y=exp,fill=class)) + 
	geom_violin(size=0.5) +
	scale_fill_manual(values=cols) +
	theme(legend.position="none") +
	theme_bw() + #scale_x_continuous(breaks = NULL) +
	xlab("") + ylab ("Dependent score (60:40)") +
	geom_segment(x = 1 - 0.25, xend = 1 + 0.25, yend = qth.1[1], y = qth.1[1], color = "grey20", linetype = 1) +
	geom_segment(x = 2 - 0.25, xend = 2 + 0.25, yend = qth.1[2], y = qth.1[2], color = "grey20", linetype = 1) +
	geom_segment(x = 3 - 0.25, xend = 3 + 0.25, yend = qth.1[3], y = qth.1[3], color = "grey20", linetype = 1) +
	geom_segment(x = 4 - 0.25, xend = 4 + 0.25, yend = qth.1[4], y = qth.1[4], color = "grey20", linetype = 1) +
	geom_segment(x = 5 - 0.25, xend = 5 + 0.25, yend = qth.1[5], y = qth.1[5], color = "grey20", linetype = 1) +
	geom_segment(x = 6 - 0.25, xend = 6 + 0.25, yend = qth.1[6], y = qth.1[6], color = "grey20", linetype = 1) +
	geom_segment(x = 7 - 0.25, xend = 7 + 0.25, yend = qth.1[7], y = qth.1[7], color = "grey20", linetype = 1) +
	geom_segment(x = 8 - 0.25, xend = 8 + 0.25, yend = qth.1[8], y = qth.1[8], color = "grey20", linetype = 1) +
	geom_segment(x = 1 - 0.25, xend = 1 + 0.25, yend = qth2.1[1], y = qth2.1[1], color = "grey20", linetype = 1) +
	geom_segment(x = 2 - 0.25, xend = 2 + 0.25, yend = qth2.1[2], y = qth2.1[2], color = "grey20", linetype = 1) +
	geom_segment(x = 3 - 0.25, xend = 3 + 0.25, yend = qth2.1[3], y = qth2.1[3], color = "grey20", linetype = 1) +
	geom_segment(x = 4 - 0.25, xend = 4 + 0.25, yend = qth2.1[4], y = qth2.1[4], color = "grey20", linetype = 1) +
	geom_segment(x = 5 - 0.25, xend = 5 + 0.25, yend = qth2.1[5], y = qth2.1[5], color = "grey20", linetype = 1) +
	geom_segment(x = 6 - 0.25, xend = 6 + 0.25, yend = qth2.1[6], y = qth2.1[6], color = "grey20", linetype = 1) +
	geom_segment(x = 7 - 0.25, xend = 7 + 0.25, yend = qth2.1[7], y = qth2.1[7], color = "grey20", linetype = 1) +
	geom_segment(x = 8 - 0.25, xend = 8 + 0.25, yend = qth2.1[8], y = qth2.1[8], color = "grey20", linetype = 1)
print(p)

dev.off()


## Fig. 3C
eff <- ef.sel$eff
sel <- ef.sel$sel
eff.1 <- eff[match(gids,names(qth))]
sel.1 <- sel[match(gids,names(qth))]

xhist <- hist(eff,breaks=500,col="pink")
yhist <- hist(sel,breaks=500,col="pink")
top <- max(xhist$counts,yhist$counts)

gids.all <- rownames(sc)

sel.th <- find.threshold(sel,
	xlab="dependency score",plot=TRUE,
	main="", ps=1e-3,
	side.fit="lower",side.thres="upper")


setwd(plot_dir);setwd("04_efficacy_selectivity_relationship")
pdf("efficacy_selectivity_p1.pdf",width=5,height=5)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

par(mar = c(4,4,1,1))
cols <- densCols(data.frame(eff,sel),colramp = colorRampPalette(blues9[-(1:3)]))
plot(eff,sel,pch=20,col=cols,type="n",
	xlab="Efficacy",ylab="Selectivity")#,cex.lab=3,cex.axis=3)
abline(h=0,v=0,lty=1,col="grey80")
abline(v=sc.ths$r0.6[["0.001"]],col="grey20",lty=2)
abline(h=sel.th,col="grey20",lty=2)
points(eff,sel,pch=20,col=cols,cex=.5)#,cex=1.5)
cols <- brewer.pal(11,"Spectral")[c(4,7)]
points(eff.1,sel.1,bg=rep(cols,each=4),col=1,cex=1,pch=21)

if(0){
	xys1 <- locator(15,type="l",col=4)
	xys2 <- locator(15,type="l",col=4)
	library(sp)
	i1 <- which(point.in.polygon(eff,sel,xys1$x,xys1$y)==1)
	i2 <- which(point.in.polygon(eff,sel,xys2$x,xys2$y)==1)	
	i12 <- c(i1,i2)
	i12 <- i12[which(eff[i12] < 0)]

	setwd(rda_dir)
	save(xys1,xys2,i1,i2,i12,file="polygons_efsel_r0.6_q1.rda")
	# saveRDS(i12,file="index_selected_on_efsel_plot_v2.rds")
}else if(0){
	setwd(rda_dir)
	load(file="polygons_efsel_r0.6_q1.rda") # xys1,xys2,i1,i2,i12
}

par(mar = c(0,4,1,1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,col=c("grey","grey50")[(xhist$breaks< sc.ths$r0.6[["0.001"]])+1],border=NA,lwd=.25)

par(mar = c(4,0,1,1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0,col=c("grey","grey50")[(yhist$breaks > sel.th)+1],border=NA,lwd=.25,horiz=TRUE)
dev.off()

##
pdf("efficacy_selectivity_p2.pdf",width=5,height=5)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

par(mar = c(4,4,1,1))
plot(eff,sel,pch=20,col=cols,type="n",
	xlab="",ylab="",
	axes=F)#,cex.lab=3,cex.axis=3)

syms12 <- unlist(mget(gids.all[i12],org.Hs.egSYMBOL))
text(eff[i12],sel[i12],syms12,cex=.6)


dev.off()


##
setwd(rda_dir)
# save(coefs,ef.sels,file="eff-sel_6thres.rda")
load("eff-sel_6thres.rda")
ef.sel <- ef.sels$r0.6$q1

for(k in 1:6){
	ef.sel <- ef.sels[[k]]$q1

	## Fig. 3C
	eff <- ef.sel$eff
	sel <- ef.sel$sel
	eff.1 <- eff[match(gids,names(qth))]
	sel.1 <- sel[match(gids,names(qth))]

	xhist <- hist(eff,breaks=500,col="pink")
	yhist <- hist(sel,breaks=500,col="pink")
	top <- max(xhist$counts,yhist$counts)

	gids.all <- rownames(sc)

	sel.th <- find.threshold(sel,
		xlab="dependency score",plot=TRUE,
		main="", ps=1e-3,
		side.fit="lower",side.thres="upper")

	setwd(plot_dir);setwd("08_comparison")
	pdf("efficacy_selectivity_p1.pdf",width=5,height=5)
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

	par(mar = c(4,4,1,1))
	cols <- densCols(data.frame(eff,sel),colramp = colorRampPalette(blues9[-(1:3)]))
	plot(eff,sel,pch=20,col=cols,type="n",
		xlab="Efficacy",ylab="Selectivity")#,cex.lab=3,cex.axis=3)
	abline(h=0,v=0,lty=1,col="grey80")
	abline(v=sc.ths$r0.6[["0.001"]],col="grey20",lty=2)
	abline(h=sel.th,col="grey20",lty=2)
	points(eff,sel,pch=20,col=cols,cex=.5)#,cex=1.5)
	cols <- brewer.pal(11,"Spectral")[c(4,7)]
	points(eff.1,sel.1,bg=rep(cols,each=4),col=1,cex=1,pch=21)

	if(0){
		xys1 <- locator(15,type="l",col=4)
		xys2 <- locator(15,type="l",col=4)
		library(sp)
		i1 <- which(point.in.polygon(eff,sel,xys1$x,xys1$y)==1)
		i2 <- which(point.in.polygon(eff,sel,xys2$x,xys2$y)==1)	
		i12 <- c(i1,i2)
		i12 <- i12[which(eff[i12] < 0)]

		setwd(rda_dir)
		save(xys1,xys2,i1,i2,i12,file="polygons_efsel_r0.6_q1.rda")
		# saveRDS(i12,file="index_selected_on_efsel_plot_v2.rds")
	}else if(0){
		setwd(rda_dir)
		load(file="polygons_efsel_r0.6_q1.rda") # xys1,xys2,i1,i2,i12
	}

	par(mar = c(0,4,1,1))
	barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,col=c("grey","grey50")[(xhist$breaks< sc.ths$r0.6[["0.001"]])+1],border=NA,lwd=.25)

	par(mar = c(4,0,1,1))
	barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0,col=c("grey","grey50")[(yhist$breaks > sel.th)+1],border=NA,lwd=.25,horiz=TRUE)
	dev.off()

	##
	pdf("efficacy_selectivity_p2.pdf",width=5,height=5)
	nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

	par(mar = c(4,4,1,1))
	plot(eff,sel,pch=20,col=cols,type="n",
		xlab="",ylab="",
		axes=F)#,cex.lab=3,cex.axis=3)

	syms12 <- unlist(mget(gids.all[i12],org.Hs.egSYMBOL))
	text(eff[i12],sel[i12],syms12,cex=.6)


	dev.off()

}






