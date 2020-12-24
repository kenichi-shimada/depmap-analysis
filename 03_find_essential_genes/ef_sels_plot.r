setwd(rda_dir)
load("eff-sel_6thres.rda")
n.dep.cells <- readRDS(file="n.dep.cells_v3.rds")
n.lin.dep <- readRDS(file="n.lin.dep_v3.rds")

setwd(src_dir)
source("functions.r")

for(k in 1:6){
	ef.sel <- ef.sels[[k]]$q1

	## Fig. 3C
	eff <- ef.sel$eff
	sel <- ef.sel$sel
	# eff.1 <- eff[match(gids,names(qth))]
	# sel.1 <- sel[match(gids,names(qth))]

	gids.all <- rownames(sc)

	eff.th <- sc.ths[[k]][["0.001"]]
	sel.th <- find.threshold(sel,
		xlab="dependency score",plot=TRUE,
		main="", ps=1e-3,
		side.fit="lower",side.thres="upper")

	if(0){
		xhist <- hist(eff,breaks=500,col="pink")
		yhist <- hist(sel,breaks=500,col="pink")
		top <- max(xhist$counts,yhist$counts)

		setwd(plot_dir);setwd("08_comparison")
		pdf(paste0("efficacy_selectivity_",k,",.pdf"),width=5,height=5)
		nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(6,1), c(1,6), TRUE)

		par(mar = c(4,4,1,1))
		cols <- densCols(data.frame(eff,sel),colramp = colorRampPalette(blues9[-(1:3)]))
		plot(eff,sel,pch=20,col=cols,type="n",
			xlab="Efficacy",ylab="Selectivity")#,cex.lab=3,cex.axis=3)
		abline(h=0,v=0,lty=1,col="grey80")
		abline(v=eff.th,col="grey20",lty=2)
		abline(h=sel.th,col="grey20",lty=2)
		points(eff,sel,pch=20,col=cols,cex=.5)#,cex=1.5)
		cols <- brewer.pal(11,"Spectral")[c(4,7)]
		# points(eff.1,sel.1,bg=rep(cols,each=4),col=1,cex=1,pch=21)

		par(mar = c(0,4,1,1))
		barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,col=c("grey","grey50")[(xhist$breaks < eff.th)+1],border=NA,lwd=.25)

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

	ndc <- n.dep.cells[,k]

	setwd(plot_dir);setwd("08_comparison")
	pdf(paste0("n_dep-cells_",k,",.pdf"),width=5,height=5)
	par(mar = c(5,5,3,3))
	min.n <- round(387*0.01)
	cols <- c(rep("grey70",min.n),
		colorRampPalette(brewer.pal(11,"Spectral"))(max(ndc)-min.n+1))
	o <- order(ndc)
	plot(ef.sel[o,],pch=20,col=cols[ndc[o]+1],type="n",
		xlab="Efficacy",ylab="Selectivity",
		main=paste0("k=",k))#,cex.lab=3,cex.axis=3)
	abline(h=0,v=0,lty=1,col="grey80")
	abline(v=eff.th,col="grey20",lty=2)
	points(ef.sel[o,],pch=20,col=cols[ndc[o]+1],cex=.5)#,cex=1.5)

	dev.off()

	##
	pdf(paste0("n_lin_dep_",k,",.pdf"),width=5,height=5)
	par(mar = c(5,5,3,3))
	cols <- c("grey70",colorRampPalette(brewer.pal(11,"Spectral"))(max(n.lin.dep[,k])))
	o <- order(n.lin.dep[,k])
	plot(ef.sel[o,],pch=20,col=cols[n.lin.dep[,k][o]+1],type="n",
		xlab="Efficacy",ylab="Selectivity",
		main=paste0("k=",k))#,cex.lab=3,cex.axis=3)
	abline(h=0,v=0,lty=1,col="grey80")
	abline(v=eff.th,col="grey20",lty=2)
	points(ef.sel[o,],pch=20,col=cols[n.lin.dep[,k][o]+1],cex=.5)#,cex=1.5)

	dev.off()

}

##
if(0){
	pdf("n_dep-cells_p2.pdf",width=5,height=5)
	par(mar=c(5,5,3,3))
	image(min.n:max(ndc),1,matrix(min.n:max(ndc),ncol=1),col=cols[-seq(min.n)],axes=F,
		xlab="",ylab="")
	box()
	axis(1,at=c(min.n+1,seq(75,300,75),387),las=1)
	dev.off()
}

##
pdf("eff_ecdf_v3.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
cols <- rev(brewer.pal(6,"Spectral"))
eff <- sapply(1:6,function(k)ef.sels[[k]]$q1$eff)
plot.ecdf(eff[,1],col=NA,verticals=T,pch=NA,
	xlab="Efficacy (1st percentile score)",
	ylab="Empirical CDF",
	main="")
for(i in 1:6){
	plot.ecdf(eff[,i],col=cols[i],verticals=T,add=T,pch=NA)	
}
for(i in 1:6){
	r <- sum(eff[,i] < sc.ths[[i]][1])/length(eff[,i])
	points(sc.ths[[i]][1],r,pch=20,col=cols[i])
}
text(-1,0.4,"Essentiality\nthreshold")
ratio <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))
legend(-3,1,ratio,lty=1,col=cols,pch=20)
dev.off()

##
pdf("sel_ecdf_v3.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
cols <- rev(brewer.pal(6,"Spectral"))
sel <- sapply(1:6,function(k)ef.sels[[k]]$q1$sel)
plot.ecdf(sel[,1],col=NA,verticals=T,pch=NA,
	xlab="Selectivity",
	ylab="Empirical CDF",
	main="")
abline(v=0,lty=2)
for(i in 1:6){
	plot.ecdf(sel[,i],col=cols[i],verticals=T,add=T,pch=NA)	
}
ratio <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))
legend(.6,0.5,ratio,lty=1,col=cols)
dev.off()

##
pdf("ecdf_n_dep-cells.pdf",width=5,height=5)
par(mar=c(5,5,3,3))
cols <- rev(brewer.pal(6,"Spectral"))
plot.ecdf(n.dep.cells[,1],col=NA,verticals=T,pch=NA,
	xlab="# dependent cell lines",
	ylab="Empirical CDF",
	main="")
for(i in 1:6){
	plot.ecdf(n.dep.cells[,i],col=cols[i],verticals=T,add=T,pch=NA)	
}
ratio <- rev(c("CRISPR","80:20","60:40","40:60","20:80","shRNA"))
legend(200,0.6,ratio,lty=1,col=cols)
dev.off()
