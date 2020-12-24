library(dplyr)
library(org.Hs.eg.db)

setwd(rda_dir)
x <- load(file="find_essential_genes.rda") 
sig.idx <- readRDS(file="sig.idx_6thres_v3.rds") # sig.idx

for(k in 1:6){
	c1 <- seq(0,1,.2)[k]

	setwd(rda_dir);setwd(paste0("r",c1));
	x <- load(file="ecodots_ori_v3.rda")	# mems,liks,ds,ratios,optims,types
	
	# ess.genes.all <- readRDS(file="ess_genes_6thres.rds")
	# if(k==1){
	# 	ess.genes <- unique(unlist(ess.genes.all[1:5]))
	# }else{
	# 	ess.genes <- unique(unlist(ess.genes.all))
	# }

	ess.genes <- names(which(is.ess[[k]]$q1))

	syms <- unlist(mget(ess.genes,org.Hs.egSYMBOL))

	th <- .5

	nclust.1 <-sapply(mems,function(x){
		x1 <- table(x)
		sum(x1>1)
	})
	names(nclust.1) <- ds

	##
	sig.i <- sig.idx[[k]]
	ds1 <- rev(ds[sig.i]) # c("114","80","50") # c("131", "87", "53")

	df <- data.frame(
		gid=ess.genes,sym=syms,
		m1=mems[[ds1[1]]],
		m2=mems[[ds1[2]]],
		m3=mems[[ds1[3]]],
		l1=liks[[ds1[1]]],
		l2=liks[[ds1[2]]],
		l3=liks[[ds1[3]]],
		stringsAsFactors=F)

	sum(is.na(df$m1)) # 122
	sum(is.na(df$m2)) # 18
	sum(is.na(df$m3)) # 5

	df$m1[is.na(df$m1)] <- max(df$m1,na.rm=T) + seq(sum(is.na(df$m1)))
	df$m2[is.na(df$m2)] <- max(df$m2,na.rm=T) + seq(sum(is.na(df$m2)))
	df$m3[is.na(df$m3)] <- max(df$m3,na.rm=T) + seq(sum(is.na(df$m3)))

	df$l1[is.na(df$l1)] <- 0
	df$l2[is.na(df$l2)] <- 0
	df$l3[is.na(df$l3)] <- 0

	## 1. one cluster (lo) => one cluster (hi), if not, split the lower cluster, because of the scheme, should do 2=>3, and 1=>2
	## 2. one cluster (lo) => none (hi) is not taken care of at this point

	## filling from hi to lo one lo should lead to one hi

	#########################
	## 2 => 3 ("M" => "L")
	#########################
	lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..
	table(sapply(lst23,length)) ## 690 has 1, 12 has 2, 1 has 16

	##
	m2dup.i <- names(which(sapply(lst23,length)>1))
	m2 <- df$m2
	m3 <- df$m3

	for(x in m2dup.i){
		idx <- which(m2==x)
		tmp <- paste0(m2[idx],".",m3[idx])
		m2[idx] <- tmp
	}

	tab <- sort(table(m2),decreasing=T)
	new.lab <- names(tab)
	m2.1 <- as.numeric(factor(m2,levels=new.lab))
	m2.1[m2.1 %in% as.numeric(names(which(table(m2.1)==1)))] <- NA
	m2.1[is.na(m2.1)] <- max(m2.1,na.rm=T) + seq(sum(is.na(m2.1)))

	df$m2 <- m2.1

	lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..
	table(sapply(lst23,length))


	#########################
	## 1 => 2 ("S" => "M")
	#########################
	lst12 <- tapply(df$m2,df$m1,table) # note, 0 is counted..
	table(sapply(lst12,length)) ## 690 has 1, 12 has 2, 1 has 16

	##
	m1dup.i <- names(which(sapply(lst12,length)>1))
	m1 <- df$m1
	m2 <- df$m2

	for(x in m1dup.i){
		idx <- which(m1==x)
		tmp <- paste0(m1[idx],".",m2[idx])
		m1[idx] <- tmp
	}

	tab <- sort(table(m1),decreasing=T)
	new.lab <- names(tab)
	m1.1 <- as.numeric(factor(m1,levels=new.lab))
	m1.1[m1.1 %in% as.numeric(names(which(table(m1.1)==1)))] <- NA
	m1.1[is.na(m1.1)] <- max(m1.1,na.rm=T) + seq(sum(is.na(m1.1)))

	df$m1 <- m1.1

	lst12 <- tapply(df$m2,df$m1,table) # note, 0 is counted..
	table(sapply(lst12,length))

	##
	lst12 <- tapply(df$m2,df$m1,table) # note, 0 is counted..
	lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..

	# sum(table(m1)>1)
	# sum(table(m2)>1)
	# sum(table(m3)>1)


	## cluster with one element => treated as zero (not in cluster) - can be done at the end.
	m1 <- df$m1
	old1 <- sort(table(m1),decreasing=T)# old
	# table(diff(as.numeric(names(old1))))
	new1 <- paste0("S",names(old1))
	names(new1) <- names(old1)

	##
	m2 <- df$m2
	old2 <- sort(table(m2),decreasing=T)
	# table(diff(as.numeric(names(old2))))
	new2 <- paste0("M",names(old2))
	names(new2) <- names(old2)

	##
	m3 <- df$m3
	old3 <- sort(table(m3),decreasing=T)
	# table(diff(as.numeric(names(old3))))
	new3 <- paste0("L",names(old3))
	names(new3) <- names(old3)

	##
	m1 <- factor(new1[m1],levels=new1)
	m2 <- factor(new2[m2],levels=new2)
	m3 <- factor(new3[m3],levels=new3)

	# table(m1)
	# table(m2)
	# table(m3)

	df$m1 <- m1
	df$m2 <- m2
	df$m3 <- m3

	##
	setwd(rda_dir);setwd(paste0("r",c1));
	saveRDS(df,file="ori_hierarchy_cluster_v3.rds")
}

setwd(rda_dir);setwd("r0");dim(readRDS("ori_hierarchy_cluster_v3.rds"))
setwd(rda_dir);setwd("r0.2");dim(readRDS("ori_hierarchy_cluster_v3.rds"))
setwd(rda_dir);setwd("r0.4");dim(readRDS("ori_hierarchy_cluster_v3.rds"))
setwd(rda_dir);setwd("r0.6");dim(readRDS("ori_hierarchy_cluster_v3.rds"))
setwd(rda_dir);setwd("r0.8");dim(readRDS("ori_hierarchy_cluster_v3.rds"))
setwd(rda_dir);setwd("r1");dim(readRDS("ori_hierarchy_cluster_v3.rds"))

