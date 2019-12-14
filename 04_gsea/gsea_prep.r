## this file is to prepare two objects to compute GSEA in the cloud, namely
# - df: efficacy-selectivity scores
# - msigdb.1: list object contains MSigDB pathway data

setwd(data_dir);setwd("msigdb")
gmts <- dir(pattern="gmt$")
## download the following files (require login credentials so manually download)
# [1] "c2.cp.biocarta.v7.0.entrez.gmt" "c2.cp.kegg.v7.0.entrez.gmt"    
# [3] "c2.cp.reactome.v7.0.entrez.gmt" "c5.all.v7.0.entrez.gmt"        
# [5] "c6.all.v7.0.entrez.gmt"         "h.all.v7.0.entrez.gmt"    	

msigdb <- lapply(gmts,function(gmt){
	tmp <- strsplit(readLines(gmt),"\t")
	nms <- sapply(tmp,function(x)x[1])
	if(gmt=="c6.all.v7.0.entrez.gmt"){
		nms <- paste0("ONCOGENIC_SIGNATURE_",nms)
	}
	tmp <- lapply(tmp,function(x)x[-(1:2)])
	names(tmp) <- nms
	return(tmp)
})
msigdb <- do.call(c,msigdb)

setwd(rda_dir)
save(msigdb,file="msigdb.v7.0.rda")


##
load(file="eff-sel.rda") # ef.sel,thres
all.genes <- unique(unlist(msigdb)) ## 19842

##
exp.genes <- rownames(ef.sel)
msigdb.ol <- lapply(msigdb,function(x)x[x %in% exp.genes])

##
n <- sapply(msigdb.ol,length)
n1 <- n >= 15 & n <= 500

msigdb.1 <- msigdb.ol[n1] ## 6551
table(sapply(names(msigdb.1),function(x)sub("_.+","",x)))

##
used.genes <- unique(unlist(msigdb.1)) ## 14831
df <- ef.sel[exp.genes %in% used.genes,] ## 14831 x 2

##
setwd(rda_dir)
save(df,msigdb.1,file="gsea_data.rda")
