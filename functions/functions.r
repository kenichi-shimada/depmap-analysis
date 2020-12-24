gene.ed.plot <- function(eff=eff.1,dep=dep.1,gid,ed.lims=NULL){
	e.lims <- quantile(as.vector(eff),c(0.0005,.9995),na.rm=T)
	d.lims <- quantile(as.vector(dep),c(0.0005,.9995),na.rm=T)
	if(missing(ed.lims)){
		ed.lims <- range(e.lims,d.lims)		
	}

	if(missing(gid)){
		sym <- ""
	}else{
		e.lim1 <- range(eff[gid,],na.rm=T)
		d.lim1 <- range(dep[gid,],na.rm=T)

		ed.lims <- range(ed.lims,e.lim1,d.lim1)

		gid <- gid[1]
		sym <- unlist(mget(gid,org.Hs.eg.db::org.Hs.egSYMBOL))
	}

	ed <- data.frame(e=as.vector(eff),d=as.vector(dep))
	i.w.ed <- which(!is.na(ed$e) & !is.na(ed$d))
	ed1 <- ed[i.w.ed,] ## 5334181

	bws <- c(diff(e.lims),diff(d.lims))/50
	est <- KernSmooth::bkde2D(ed1,bandwidth=bws,gridsize=rep(500,2))
	ncols <- rep(rep(c("grey80","grey50"),c(4,1)),length.out=200)[-(1:4)]
	levs <- c(0.005,0.01,0.02,0.05,0,0.2,0.5,1:4)

	uniq.cols <- c("grey80","grey70","grey50","grey20") 
	cols <- uniq.cols[c(1,2,2,2,3,3,3,4,4,4,4)]

	smoothScatter(as.vector(eff),as.vector(dep),
		transformation=function(x)x^.5,
		colramp=colorRampPalette(c("white",blues9[3:6])),
		xlab="",ylab="",asp=1,
		main=sym,xlim=ed.lims,ylim=ed.lims,pch=NA)
	mtext("CRISPR dependency score",1,3)
	mtext("shRNA dependency score",2,3)

	contour(est$x1,est$x2,est$fhat,levels=levs,add=T,col=cols)
	abline(h=0,v=0,lty=2)
	abline(a=0,b=1,lty=2)

	if(!missing(gid)){
		points(eff[gid,],dep[gid,],pch=20,col=1,cex=.5)
	}
}

find.threshold <- function(dep.scores,#breaks=100,
	side.fit=c("upper","lower"),side.thres=c("upper","lower"),
	ps=1e-3,plot=TRUE,...){

	# xh <- hist(dep.scores,col="lightblue",border=NA,breaks=seq.br,plot=TRUE,...)
	d <- density(dep.scores)
	mod <- d$x[which.max(d$y)]

	if(side.fit[1]=="upper"){
		val <- (dep.scores-mod)
		val <- val[val > 0] 
		sd <- sqrt(sum(val^2)/(length(val)-1))
	}else{
		val <- (dep.scores-mod)
		val <- val[val < 0] 
		sd <- sqrt(sum(val^2)/(length(val)-1))
	}

	if(!is.na(ps)){
		if(side.thres[1]=="upper"){
			ths <- qnorm(1-ps,mean=mod,sd=sd)
		}else{
			ths <- mod * 2 - qnorm(1-ps,mean=mod,sd=sd)
		}
		names(ths) <- ps
	}
	if(plot){
		plot(d$x,d$y,type="n",...)
		polygon(d$x,d$y,col="lightblue",border=NA)
		ybr <- dnorm(d$x,mean=mod,sd=sd)
		f <- max(d$y)/max(ybr)
		lines(d$x,f*ybr,col=4,lwd=2)
		if(!is.na(ps))abline(v=ths,col=2)
	}

	return(ths)
}

