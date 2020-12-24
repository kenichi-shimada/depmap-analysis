library(ADaM2)
library(dplyr)
library(org.Hs.eg.db)
library(impute)
is.ess <- list()
eff.ths <- rep(NA,6)

## Load data
setwd(rda_dir)
load(file="eff-sel_6thres.rda") # coefs,ef.sels
load(file="find_essential_genes.rda") # is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns
cl.info <- readRDS("celllines_19q3.rds") 
%>% filter(DepMap_ID %in% colnames(is.ess.imp[[1]]))

## impute essentiality matrix
is.ess.imp <- parallel::mclapply(1:6,function(k){
  ef.sel <- ef.sels[[k]]

  c1 <- seq(0,1,.2)[k]
  setwd(rda_dir);setwd(paste0("r",c1))
  scores <- readRDS(file=paste0("scores_15847g_423c.rds"))
  imp.scores <- impute.knn(scores)$data

  ess.mat <- array(as.numeric(imp.scores < sc.ths[[k]][1]),dim(imp.scores))
  rownames(ess.mat) <- rownames(scores)
  colnames(ess.mat) <- colnames(scores)

  return(ess.mat)
},mc.cores=6)

cl.info <- cl.info %>% filter(DepMap_ID %in% colnames(scores))
colnames(is.ess.imp[[1]]))


## Remove lineages with less than 10 cell lines
lin2cell <- tapply(cl.info$DepMap_ID,cl.info$lineage,identity)
n.lin2cell <- sapply(lin2cell,length)

hist(n.lin2cell,breaks=80,col="grey50")
abline(v=5,col=2)
lin <- names(lin2cell)[n.lin2cell >= 10]

lin2cell <- lin2cell[lin]
n.lin2cell <- n.lin2cell[lin]



## Load reference genes
data(curated_BAGEL_essential)
eids.bagel.ess <- unlist(mget(curated_BAGEL_essential,org.Hs.egSYMBOL2EG,ifnotfound=NA)) # 326, 8 NAs
eids.bagel.ess <- eids.bagel.ess[!is.na(eids.bagel.ess)] # 318

# Compute crossoverpoints
setwd(plot_dir);setwd("adam2")

crossoverpoints <- list()
for(k in 1:6){
  crossoverpoints[[k]] <-  sapply(lin,function(ti){
    cat("*")
    cells <- lin2cell[[ti]]
    this.ess <- is.ess.imp[[k]][,cells,drop=F]

    pdf(paste0("bagel_",ti,"_",k,".pdf"),width=5,height=5)
    par(mar=c(5,5,3,3))
    pprofile<-ADAM2.panessprofile(depMat=this.ess)
    invisible(capture.output(nullmodel<-ADAM2.generateNullModel(depMat=this.ess,ntrials = 1000)))
    EO<-ADAM2.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )
    TPR<-ADAM2.truePositiveRate(this.ess,eids.bagel.ess)
    # TPR<-ADAM2.truePositiveRate(this.ess,eids.ess)      
    at <- ADAM2.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'essential genes')
    dev.off()

    return(at)
  })
}

crossoverpoints <- do.call(cbind,crossoverpoints)
rownames(crossoverpoints) <- lin


## number of essential genes per lineage
n.lin.ess <- parallel::mclapply(1:6,function(k){
  sapply(lin,function(ti){
    cells <- lin2cell[[ti]]
    this.ess <- is.ess.imp[[k]][,cells,drop=F]
    n.ess <- rowSums(this.ess)
    return(n.ess)
  })
},mc.cores=6)



##
setwd(plot_dir);setwd("adam2")
lin.crossoverpoints <- list()
for(k in 1:6){
  n.lin <- n.lin.ess[[k]]
  cp <- crossoverpoints[,k]
  n.ess <- t(apply(n.lin,1,function(n){
    as.numeric(n >= cp)
  }))
  colnames(n.ess) <- lin

  pdf(paste0("lin_",k,".pdf"))
  pprofile<-ADAM2.panessprofile(depMat=n.ess)
  invisible(capture.output(nullmodel<-ADAM2.generateNullModel(depMat=n.ess,ntrials = 1000)))
  EO<-ADAM2.empiricalOdds(observedCumSum = pprofile$CUMsums,simulatedCumSum =nullmodel$nullCumSUM )
  # TPR<-ADAM2.truePositiveRate(n.ess,eids.ess)
  TPR<-ADAM2.truePositiveRate(n.ess,eids.bagel.ess)  
  lin.cp <- ADAM2.tradeoffEO_TPR(EO,TPR$TPR,test_set_name = 'essential genes')
  dev.off()

  lin.crossoverpoints[[k]] <- lin.cp
}

lin.crossoverpoints <- do.call(c,lin.crossoverpoints)



## Find commonly essential genes
is.comm <- list()
for(k in 1:6){
  n.lin <- n.lin.ess[[k]]
  cp <- crossoverpoints[,k]
  n.ess <- t(apply(n.lin,1,function(n){
    as.numeric(n >= cp)
  }))
  colnames(n.ess) <- lin
  is.comm[[k]] <- as.numeric(rowSums(n.ess) >= lin.crossoverpoints[k])
}
is.comm <- do.call(cbind,is.comm)
rownames(is.comm) <- rownames(n.ess)
colnames(is.comm) <- paste0("r",seq(0,1,.2))

## Save data
setwd(rda_dir)
saveRDS(is.comm,file="is.comm.rds")
save(crossoverpoints,lin.crossoverpoints,file="crossoverpoints_bagel_ge10.rda")
save(lin,lin2cell,n.lin2cell,file="lineage_cells.rda")
