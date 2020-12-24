setwd(src_dir);source("functions.r")

## updated version: allowing multiple quantile threshold (qtl) to define efficacy/selectivity

## cutoff values (quantiles) for efficacy and selectivity
qtl <- c(0.5,1,2.5,5,10,25)
qtl2 <- 100-qtl
cr.ratio <- paste0("r",seq(0,1,.2))

 ## find essential and tumor suppressor genes
sc.ths <- is.ts <- is.ess <- list()
n.ts.gns <- n.ess.gns <- array(NA,c(length(cr.ratio),length(qtl)))


for(k in 1:6){
	is.ts[[k]] <- is.ess[[k]] <- list()

	## load the dependency matrix
	c1 <- seq(0,1,.2)[k]
	setwd(rda_dir);setwd(paste0("r",c1))
	ef <- readRDS(file=paste0("scores_15847g_423c.rds"))

	## compute thresholds for essential and tumor suppressor genes
	sc <- as.vector(ef)
	sc <- sc[!is.na(sc)]

	p.th <- 1e-3
	sc.ths[[k]] <- find.threshold(sc,
		xlab="dependency score",plot=FALSE,
		main="", ps=c(p.th,1-p.th),
		side.fit="upper",side.thres="lower")

	## 
	setwd(plot_dir);setwd("03_histogram_essential_genes")

	pdf(paste0("lower_",cr.ratio[k],".pdf"),width=5,height=5)
	par(mar=c(5,5,3,3))

	for(i in 1:6){
		## efficacy -qtl[i]-th quantile score
		qth <- apply(ef,1,function(x)quantile(x,qtl[i]/100,na.rm=T))

		##
		ths <- find.threshold(qth,
			xlab="dependency score",plot=TRUE,
			main=paste0(qtl[i]," percentile"), ps=NA,
			side.fit="upper",side.thres="lower",
			xlim=c(-3,2))
		abline(v=sc.ths[[k]],col=2)

		is.ess[[k]][[i]] <- idx <- qth < sc.ths[[k]][1]
		n.ess.gns[k,i] <- n.genes <- sum(idx)

		# text(ths,as.vector(par()$usr[3:4]%*%c(.25,.75)),
		# 	paste("p< ",p.th," (",n.genes,"genes)"),srt=90,pos=3,cex=.5)
	}
	dev.off()

	##
	pdf(paste0("upper_",cr.ratio[k],".pdf"),width=5,height=5)
	par(mar=c(5,5,3,3))

	for(i in 1:6){
		qth2 <- apply(ef,1,function(x)quantile(x,qtl2[i]/100,na.rm=T))
		ths <- find.threshold(qth2,
			xlab="dependency score",plot=TRUE,
			main=paste(qtl2[i]," percentile"),ps=NA,
			side.fit="upper",side.thres="upper",
			xlim=c(-3,2))
		abline(v=sc.ths[[k]],col=2)

		is.ts[[k]][[i]] <- idx <- qth2 > sc.ths[[k]][2]
		n.ts.gns[k,i] <- n.genes <- sum(idx)

		# text(ths,as.vector(par()$usr[3:4]%*%c(.25,.75)),
		# 	paste("p<",p.th," (",n.genes,"genes)"),srt=90,pos=3,cex=.5)
	}
	dev.off()
}

names(sc.ths) <- names(is.ess) <- names(is.ts) <- cr.ratio
is.ess <- lapply(is.ess,function(x){
	names(x) <- paste0("q",qtl)
	return(x)
})
is.ts <- lapply(is.ts,function(x){
	names(x) <- paste0("q",qtl)
	return(x)
})
rownames(n.ess.gns) <- rownames(n.ts.gns) <- cr.ratio
colnames(n.ess.gns) <- colnames(n.ts.gns) <- paste0("q",qtl)

# is.ess <- ess.genes <- coef <- list()
# system.time(
# 	for(i in 1:6){
# 		cat("*")
# 		ess.genes[[i]] <- names(which(is.ess[[i]]))
# 		coef[[i]] <- cor(t(ef[is.ess[[i]],]),method="spearman",use="pairwise.complete.obs")
# 	}
# )

## Genes that are neither essential nor growth-suppressing
## take p=1e-3 as the threshold

lo.tab <- sapply(is.ess,function(x)x[[2]])
hi.tab <- sapply(is.ts,function(x)x[[2]])

is.non <- lapply(cr.ratio,function(r)!is.ess[[r]]$q1 & !is.ts[[r]]$q1)
names(is.non) <- cr.ratio
sapply(is.non,sum)

## saving objects
setwd(rda_dir)
save(is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns,
	file="find_essential_genes.rda")
