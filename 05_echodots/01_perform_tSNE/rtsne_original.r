library(Rtsne) # 0.15

setwd(rda_dir)

ess.genes <- readRDS(file="ess_2492_genes.rds")
ef <- readRDS(file="scores_15847g_423c_101119.rds")[ess.genes,]

i <- as.numeric(commandArgs(TRUE)[1]) ## i in 1:200

ef.all <- ef
nrows <- apply(ef.all,1,function(x)sum(is.na(x)))
no.na <- nrows==0

# Algorithm for ECODOTs
# coef <- cor(t(ef.1),t(ef.2),method="spearman",use="pairwise.complete.obs")
coef.0 <- cor(t(ef.all[no.na,]),t(ef.all[no.na,]),method="spearman",use="everything")# 2696 x 2696
coef.1 <- cor(t(ef.all[!no.na,]),t(ef.all[no.na,]),method="spearman",use="pairwise.complete.obs") # 45 x 2696
coef.2 <- cor(t(ef.all[no.na,]),t(ef.all[!no.na,]),method="spearman",use="pairwise.complete.obs") # 2696 x 45
coef.3 <- cor(t(ef.all[!no.na,]),t(ef.all[!no.na,]),method="spearman",use="pairwise.complete.obs") # 45 x 45

coef <- array(NA,rep(nrow(ef.all),2))
colnames(coef) <- rownames(coef) <- rownames(ef.all)

coef[no.na,no.na] <- coef.0
coef[!no.na,no.na] <- coef.1
coef[no.na,!no.na] <- coef.2
coef[!no.na,!no.na] <- coef.3

# Algorithm for ECHODOTs
n.coef <- coef
n.coef[n.coef > 1] <- 1
n.coef[n.coef < -1] <- -1
d <- 1 - n.coef

## Rtsne
set.seed(i)
system.time(rt.outs <- Rtsne(X=d,is_distance=TRUE,perplexity=5,max_iter=20000,pca=FALSE,theta=0.5))
rt <- rt.outs$Y
rownames(rt) <- rownames(d)

setwd(rda_dir);setwd("rtsne_original")
saveRDS(rt,file=paste0("rtsne_2492gns_20000_original_",i,".rds"))
