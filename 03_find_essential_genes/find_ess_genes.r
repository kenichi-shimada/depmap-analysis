k <- as.numeric(commandArgs(TRUE))[1] # 1-12

c1 <- seq(0,1,.2)[k]

setwd(rda_dir)
ef <- readRDS(file=paste0("scores_15847g_423c_",c1,".rds"))

## cutoff values (quantiles) for efficacy and selectivity
qtl <- c(0.5,1,2.5,5,10,25)
qtl2 <- 100-qtl
qths <- lapply(qtl,function(th){
	apply(ef,1,function(x)quantile(x,th/100,na.rm=T))
}) ## essential genes
qths2 <- lapply(qtl2,function(th){
	apply(ef,1,function(x)quantile(x,th/100,na.rm=T))
}) ## tumor suppressor genes
rths <- lapply(1:6,function(i)qths2[[i]]-qths[[i]])

names(qths) <- paste0("q",qtl)
names(qths2) <- paste0("q",qtl2)
names(rths) <- paste0("r",qtl2-qtl)

##
r1 <- range(unlist(qths),na.rm=T)
r2 <- range(unlist(rths),na.rm=T)

## Identification of essential and growth-suppressing genes
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

## create directories for saving plots
setwd(plot_dir);setwd("03_histogram_essential_genes")
for(c1 in paste0("r",seq(0,1,.2)))dir.create(c1)



## thresholds for essential genes
lo.ths <- lo.idxs <- lo.ngns <- list()

for(i in 1:6){
	## fig.
	qth <- qths[[i]]

	## Fig. 3A, S2
	setwd(plot_dir)
	setwd("03_histogram_essential_genes")
	setwd(paste0("r",c1))

	pdf(paste0("lower_",i,".pdf"),width=5,height=5)
	par(mar=c(5,5,3,3))

	lo.ths[[i]] <- ths <- find.threshold(qth,
		xlab="dependency score",plot=TRUE,
		main=paste(qtl[i]," percentile"),
		side.fit="upper",side.thres="lower")

	lo.ngns[[i]] <- n.genes <- sapply(ths,function(n)sum(qth < n))

	lo.idxs[[i]] <- idx <- cbind(th2=qth < ths[1],
		th3=qth < ths[2],
		th5=qth < ths[3])

	text(ths,as.vector(par()$usr[3:4]%*%c(.25,.75)),
		paste("p<",sig.ps," (",n.genes,"genes)"),srt=90,pos=3,cex=.5)
	dev.off()

}

names(lo.ths) <- names(lo.idxs) <- names(lo.ngns) <- paste0("q",qtl)

n.ess.gns <- sapply(lo.ngns,function(x){names(x) <- sig.ps;return(x)})

is.ess <- ess.genes <- coef <- list()
system.time(
	for(i in 1:6){
		cat("*")
		is.ess[[i]] <- lo.idxs[[i]][,2]
		ess.genes[[i]] <- names(which(is.ess[[i]]))
		coef[[i]] <- cor(t(ef[is.ess[[i]],]),method="spearman",use="pairwise.complete.obs")
	}
)


## thresholds for tumor suppressor genes
is.ts <- ts.genes <- list()
system.time(
	for(i in 1:6){
	# for(i in 5:6){
		cat("*")
		is.ts[[i]] <- lo.idxs[[i]][,2] # q2.5, 1e-3
		ts.genes[[i]] <- names(which(is.ts[[i]]))
	}
)


## Genes that are neither essential nor growth-suppressing
## take p=1e-3 as the threshold
lo.tab <- sapply(lo.idxs,function(x)x[,2])
sum.lo.tab <- apply(lo.tab,1,sum)

hi.tab <- sapply(hi.idxs,function(x)x[,2])
sum.hi.tab <- apply(hi.tab,1,sum)

is.non <- sum.lo.tab==0 & sum.hi.tab==0 

## saving objects
setwd(rda_dir);setwd(paste0("r",c1))

saveRDS(is.non,file="non-essential-genes-6thres.rds")
saveRDS(ess.genes,file="ess_genes_6thres.rds")
saveRDS(is.ess,file="is.ess_genes_6thres.rds")
save(lo.idxs,lo.ths,n.ess.gns,file="eff-threshold.rda")
saveRDS(coef,file="eff_coef_genes_6thres.rds")

saveRDS(ts.genes,file="ts_genes_6thres.rds")
saveRDS(is.ts,file="is.ts_genes_6thres.rds")
save(hi.idxs,hi.ths,n.ts.gns,file="ts-eff-threshold.rda")

save.image("ess.genes.rda")



if(0){
	## Essential genes, ts genes
	setwd(rda_dir);setwd(paste0("r",c1))
	ess.genes <- readRDS(file="ess_genes_6thres.rds")
	is.ess <- readRDS(file="is.ess_genes_6thres.rds")
	load(file="eff-threshold.rda") # lo.idxs,lo.ths,lo.idx.sum,n.ess.gns

	ts.genes <- readRDS(file="ts_genes_6thres.rds")
	is.ts <- readRDS(file="is.ts_genes_6thres.rds")
	x <- load(file="ts-eff-threshold.rda") # hi.idxs,hi.ths,hi.idx.sum,n.ts.gns

	coef <- readRDS(file="eff_coef_genes_6thres.rds")
}


if(0){
	## gene overlap between different quantile thresholds
	setwd(plot_dir);setwd("03_histogram_essential_genes")

	lo.mod.matrix <- sapply(lo.idxs,function(x)x[,2])
	hi.mod.matrix <- sapply(hi.idxs,function(x)x[,2])

	pdf("venn_diagram_ess_genes.pdf")
	vennDiagram(lo.mod.matrix[,1:5])
	dev.off()

	pdf("venn_diagram_ts_genes.pdf")
	vennDiagram(hi.mod.matrix[,1:5])
	dev.off()
}

