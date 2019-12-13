library(fgsea)
library(BiocParallel)
register(SerialParam())

setwd(rda_dir)
load(file="gsea_data.rda") # df,msigdb.1

i <- as.numeric(commandArgs(TRUE))[1] ## 1-1345

gns <- rownames(df)
v <- df[[i]]
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
fn <- paste0("fgsea_1e7_",i,".rds")
saveRDS(fgseaRes,file=fn)
