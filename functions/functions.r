find.threshold <- function(dep.scores,#breaks=100,
	side.fit=c("upper","lower"),side.thres=c("upper","lower"),
	ps=c(1e-2,1e-3,1e-5),
	plot=TRUE,...){

		# xh <- hist(dep.scores,col="lightblue",border=NA,breaks=seq.br,plot=TRUE,...)
		d <- density(dep.scores)
		mod <- d$x[which.max(d$y)]

		# if(0){
		# 	seq.br <- seq(xlim[1],xlim[2],length=breaks)
		# 	nl <- length(xh$breaks)
		# 	br <- (xh$breaks[-1]+xh$breaks[-nl])/2

		# 	ct <- xh$counts
		# 	names(ct) <- br
		# }else{
		# 	ct <- d$y
		# 	names(ct) <- d$x
		# }

		# mod <- as.numeric(names(which.max(ct)))

		if(side.fit[1]=="upper"){
			val <- (dep.scores-mod)
			val <- val[val > 0] 
			sd <- sqrt(sum(val^2)/(length(val)-1))
			names(ths) <- ps
		}else{
			val <- (dep.scores-mod)
			val <- val[val < 0] 
			sd <- sqrt(sum(val^2)/(length(val)-1))
			ths <- qnorm(1-ps,mean=mod,sd=sd)
			names(ths) <- ps
		}

		if(side.thres[1]=="upper"){
			ths <- qnorm(1-ps,mean=mod,sd=sd)
			col <- "pink"
		}else{
			ths <- mod * 2 - qnorm(1-ps,mean=mod,sd=sd)
			col <- "lightblue"
		}
		if(plot){
			plot(d$x,d$y,type="n",...)
			polygon(d$x,d$y,col=col,border=NA)
			ybr <- dnorm(seq.br,mean=mod,sd=sd)
			f <- max(d$y)/max(ybr)
			lines(seq.br,f*dnorm(seq.br,mean=mod,sd=sd),col=4,lwd=2)
			abline(v=ths,col=2)
		}

		return(ths)
}