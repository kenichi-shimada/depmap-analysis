library(fgsea)
library(BiocParallel)
register(SerialParam())

setwd(rda_dir)
load(file="gsea_data_6thres.rda") # dfs,msigdb.1

i <- as.numeric(commandArgs(TRUE))[1] # 1-72

i1 <- names(dfs)
i2 <- names(dfs[[1]])
ids <- expand.grid(idx1=i1,idx2=i2,ef.sel=1:2)

ef.sel <- dfs[[ids$idx1[i]]][[ids$idx2[i]]]
col <- ids$ef.sel[i]

gns <- rownames(ef.sel)
v <- ef.sel[[col]]
names(v) <- gns

sorted.all.genes <- sort(v,decreasing=T)

## no filter
set.seed(12345)
fgseaRes <- fgsea(pathways = msigdb.1,
              stats = sorted.all.genes,
              minSize=15,
              maxSize=500,
              nperm=1e7)

setwd(rda_dir);setwd("fgsea")
fn <- paste0("fgsea_1e7_new_",i,".rds")
saveRDS(fgseaRes,file=fn)
