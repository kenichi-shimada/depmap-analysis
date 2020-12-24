setwd(rda_dir)
load(file="find_essential_genes.rda") # is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns

##
str(is.ess)

##
setwd(rda_dir);setwd(paste0("r",c1))
ef <- readRDS(file=paste0("scores_15847g_423c.rds"))

## compute thresholds for essential and tumor suppressor genes
c1 <- seq(0,1,.2)[k]

sc <- as.vector(ef)
sc <- sc[!is.na(sc)]

p.th <- 1e-3
sc.ths[[k]] <- find.threshold(sc,
	xlab="dependency score",plot=FALSE,
	main="", ps=c(p.th,1-p.th),
	side.fit="upper",side.thres="lower")
