library(Rtsne) # 0.15

j <- as.numeric(commandArgs(TRUE)[1]) ## i in 1:1200

ps <- expand.grid(i=1:200,k=1:6)

i <- ps$i[j]
k <- ps$k[j]
c1 <- seq(0,1,.2)[k]

setwd(rda_dir);setwd(paste0("r",c1))
coef <- readRDS(file="ess_coef.rds")
d <- 1-coef

## Rtsne
set.seed(i)
system.time(rt.outs <- Rtsne(X=d,is_distance=TRUE,perplexity=5,max_iter=20000,pca=FALSE,theta=0.5))
rt <- rt.outs$Y
rownames(rt) <- rownames(d)

setwd(rda_dir);setwd(paste0("r",c1))

if(!dir.exists("rtsne_multi_thres_1")){
	dir.create("rtsne_multi_thres_1")
}

setwd("rtsne_multi_thres_1")
saveRDS(rt,file=paste0("rtsne_ori_",i,".rds"))

