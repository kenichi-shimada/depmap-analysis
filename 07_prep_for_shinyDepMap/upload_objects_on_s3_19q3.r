## This is for preparation of a shiny app only.
## Large objects are uploaded to Amazon's cloud. A shiny-based app (shinyDepMap)
## hosted on shinyapps.io downloads the files upon request.

library(aws.s3)
library(dplyr)

if(0){
	setwd("~/Dropbox (HMS)/projects/shinyapps/shinyDepMap/shinyDepMap/data")
	xx <- load("depmap_initial_19q3_local_run.rda")
	eids.all <- df$eid

	# cl.info,df,df.0,eids.ess,thres.eff,dummy,cs,p.scores,coef,graphs

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

}else{
	for(mix.ratio in 6:1){
	    # c1 <- seq(0,1,.2)[mix.ratio]
        # setwd(rda_dir);#setwd(paste0("r",c1))
		setwd("~/Dropbox (HMS)/projects/shinyapps/shinyDepMap/shinyDepMap-online/data")
		x <- load(paste0("depmap_initial_19q3_local_run_v3_",mix.ratio,".rda"))

		# common: cl.info
		# thres dependent but consolidated: eids.ess, thres.eff
		# thres or mix.ratio dependent: dfs,mems,dummy,cs,p.scores,coef,graphs

		for(thres in 1:6){
			s3saveRDS(dfs[[thres]],object=paste0("objects_ext_19q3/df_",mix.ratio,"_",thres,".rds"),bucket="depmap")
		}
		
		if(0){
			s3saveRDS(mems,object=paste0("objects_ext_19q3_v3/mems_",mix.ratio,".rds"),bucket="depmap")
		}
		# }
		## compounds
		eids.all <- rownames(p.scores)
		system.time({
			cat(paste0(mix.ratio,"_all"))
			xx <- lapply(eids.all,function(n){
				# cat("*")
				this.ef <- p.scores[n,]
				s3saveRDS(this.ef,
					object=paste0("scores_all_19q3_v3/",mix.ratio,"/",n,".rds"),
					bucket="depmap")
			})
		})

		system.time({
			cat(paste0(mix.ratio,"_ess"))
			xx <- lapply(eids.ess,function(n){
				# cat("*")
				this.coef <- coef[n,]
				s3saveRDS(this.coef,
					object=paste0("coef_ess_19q3_v3/",mix.ratio,"/",n,".rds"),
					bucket="depmap")
			})
		})

		# m1 <- levels(mems$mem1)
		# m2 <- levels(mems$mem2)
		# m3 <- levels(mems$mem3)

		# ## mem1
		# for(this.cl in m1){
		# 	# cat("*")
		# 	nodes <- graphs$mem1[[this.cl]]$nodes
		# 	edges <- graphs$mem1[[this.cl]]$edges
		# 	s3save(nodes,edges,object=paste0("graph_19q3_v3/",mix.ratio,"/",this.cl,".rda"),bucket="depmap")
		# }
		# cat("*")

		# ## mem2
		# for(this.cl in m2){
		# 	# cat("*")
		# 	nodes <- graphs$mem2[[this.cl]]$nodes
		# 	edges <- graphs$mem2[[this.cl]]$edges
		# 	s3save(nodes,edges,object=paste0("graph_19q3_v3/",mix.ratio,"/",this.cl,".rda"),bucket="depmap")
		# }
		# cat("*")

		# ## mem3
		# for(this.cl in m3){
		# 	# cat("*")
		# 	nodes <- graphs$mem3[[this.cl]]$nodes
		# 	edges <- graphs$mem3[[this.cl]]$edges
		# 	s3save(nodes,edges,object=paste0("graph_19q3_v3/",mix.ratio,"/",this.cl,".rda"),bucket="depmap")
		# }
		# cat("*\n")
		# Sys.sleep(10)
	}
}
