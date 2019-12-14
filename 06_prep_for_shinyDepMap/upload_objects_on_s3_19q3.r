## This is for preparation of a shiny app only.
## Large objects are uploaded to Amazon's cloud. A shiny-based app (shinyDepMap)
## hosted on shinyapps.io downloads the files upon request.

library(aws.s3)
library(dplyr)

setwd(rda_dir)
xx <- load("depmap_initial_19q3_local_run.rda")
#cl.info,df,df.0,eids.ess,thres.eff,dummy,cs,p.scores,coef,graphs

eids.all <- df$eid

## compounds
system.time(
	xx <- lapply(eids.all,function(n){
		cat("*")
		this.ef <- p.scores[n,]
		s3saveRDS(this.ef,object=paste0("scores_all_19q3/",n,".rds"),bucket="depmap")
	})
)
system.time(
	xx <- lapply(eids.ess,function(n){
		cat("*")
		this.coef <- coef[n,]
		s3saveRDS(this.coef,object=paste0("coef_ess_19q3/",n,".rds"),bucket="depmap")
	})
)

##
m1 <- levels(df.0$mem1)
m2 <- levels(df.0$mem2)
m3 <- levels(df.0$mem3)

## mem1
for(this.cl in m1){
	cat("*")
	nodes <- graphs$mem1[[this.cl]]$nodes
	edges <- graphs$mem1[[this.cl]]$edges
	s3save(nodes,edges,object=paste0("graph_19q3/",this.cl,".rda"),bucket="depmap")
}

## mem2
for(this.cl in m2){
	cat("*")
	nodes <- graphs$mem2[[this.cl]]$nodes
	edges <- graphs$mem2[[this.cl]]$edges
	s3save(nodes,edges,object=paste0("graph_19q3/",this.cl,".rda"),bucket="depmap")
}

## mem3
for(this.cl in m3){
	cat("*")
	nodes <- graphs$mem3[[this.cl]]$nodes
	edges <- graphs$mem3[[this.cl]]$edges
	s3save(nodes,edges,object=paste0("graph_19q3/",this.cl,".rda"),bucket="depmap")
}

delete_object("graph_19q3/NA.rda",bucket="depmap")
