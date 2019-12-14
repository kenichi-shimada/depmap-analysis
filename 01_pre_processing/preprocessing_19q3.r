## information of cell lines
cells <- read.csv("https://ndownloader.figshare.com/files/16757723",stringsAsFactors=F)


## CRISPR efficacy
eff <- read.csv("https://ndownloader.figshare.com/files/16757666",stringsAsFactors=F)
cells.ind <- eff[[1]]
rownames(eff) <- cells.ind

col.eff <- colnames(eff)[-1]
col.eff <- sub("\\.\\.([0-9]+)\\.","_\\1",col.eff) ## gene symbol doesn't contain '_'

sym.eff <- sub("_.+","",col.eff)
eids.eff <- sub(".+_","",col.eff)

eff <- eff[,-1] ## all unique genes
colnames(eff) <- eids.eff
eff <- t(eff)

eids.eff <- rownames(eff)


## shRNA efficacy
dep <- read.csv("https://ndownloader.figshare.com/files/13515395",stringsAsFactors=F,row.names=1) ## 17212 x 712

colnames(dep) <- sub("^X([0-9])","\\1",colnames(dep))
colnames(dep)[which(!colnames(dep) %in% cells$CCLE)] ## three cell lines are orphan: SW527, GISTT1, AZ521,MB157
grep("AZ521",cells$CCLE,value=T)	
colnames(dep)[grep("AZ521_STOMACH",colnames(dep))] <- "AZ521_SMALL_INTESTINE"

bid <- cells$DepMap_ID[match(colnames(dep),cells$CCLE)]
colnames(dep) <- bid
dep <- dep[,!is.na(bid)] # 710

row.dep <- rownames(dep) ## gene symbol doesn't contain '_'
row.dep <- sub(" \\(([0-9&]+)\\)",":\\1",row.dep)

eids.dep <- sub(".+:","",row.dep)
list.eids.dep <- strsplit(eids.dep,"&")
uniq.ed <- unique(unlist(list.eids.dep)) ## all unique
table(table(unlist(list.eids.dep))) ## 17731

syms.dep <- sub(":.+","",row.dep) 
list.syms.dep <- strsplit(syms.dep,"&")
uniq.sd <- unique(unlist(list.syms.dep)) ## all unique
table(table(unlist(list.syms.dep))) 

table(sapply(list.syms.dep,length)) ## 17309

rownames(dep) <- eids.dep
dep <- as.matrix(dep) ## 17309 x 710
grep("&",rownames(dep)) ## what should I do?

eids.dep <- rownames(dep)


## merge eff and dep: keep genes and cell lines in which CRISPR and shRNA data are tested
table(table(eids.dep)) # 17309

eids.dep.1 <- unlist(strsplit(eids.dep,"&"))
table(table(eids.dep.1)) # 17731
eids.1 <- eids.dep.1[!eids.dep.1 %in% eids.dep] ## 708
comm.eids.1 <- intersect(eids.eff,eids.1) 

##
comm.eids <- intersect(eids.eff,eids.dep) ## 15847
comm.cells <- intersect(colnames(eff),colnames(dep)) ## 423

eff.1 <- eff[comm.eids,comm.cells] ## 15847 x 423
dep.1 <- dep[comm.eids,comm.cells] ## 15847 x 423


## save RData files
setwd(rda_dir)
saveRDS(cells,file="celllines_19q3.rds")	
saveRDS(eff,file="gene_effect_19q3.rds")
saveRDS(dep,file="gene_dep_19q3.rds")
save(eff.1,dep.1,file="eff_dep_19q3.rda")

