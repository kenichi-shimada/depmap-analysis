library(dplyr)
library(org.Hs.eg.db)

setwd(rda_dir)
# load(file="clue/ecodots_various_d.rda")
load(file="ecodots_original.rda")	 # mems,liks,ds,ratios,optims,types
ess.genes <- readRDS(file="ess_2492_genes.rds")
syms <- unlist(mget(ess.genes,org.Hs.egSYMBOL))

th <- .5

ds <- c("131", "87", "53")

df <- data.frame(
	gid=ess.genes,sym=syms,
	m1=as.character(mems[[ds[1]]]),
	m2=as.character(mems[[ds[2]]]),
	m3=as.character(mems[[ds[3]]]),
	l1=liks[[ds[1]]],
	l2=liks[[ds[2]]],
	l3=liks[[ds[3]]],
	stringsAsFactors=F)

df$l1[is.na(df$l1)] <- 0
df$l2[is.na(df$l2)] <- 0
df$l3[is.na(df$l3)] <- 0

sum(is.na(df$m1)) # 122
sum(is.na(df$m2)) # 18
sum(is.na(df$m3)) # 5
## task: splitting clusters if all the cluster members at lower level don't correspond to the same cluster at higher level.
## checking likelihood distribution
if(0){
	sapply(df[c("l1","l2","l3")],function(x)table(x >= th))
	plot.ecdf(df$l1,pch=NA,verticals=T,col=1,xlim=c(0,1))
	plot.ecdf(df$l2,add=T,pch=NA,verticals=T,col=2)
	plot.ecdf(df$l3,add=T,pch=NA,verticals=T,col=3)
	abline(v=0.5,lty=2)
}

## 1. one cluster (lo) => one cluster (hi), if not, split the lower cluster, because of the scheme, should do 2=>3, and 1=>2
## 2. one cluster (lo) => none (hi) is not taken care of at this point

## filling from hi to lo one lo should lead to one hi
#########################
## 3 => 2 ("54" => "88")
#########################
df.0 <- df
# df <- df.0

lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..
table(sapply(lst23,length)) ## 1 has 0, 382 has 1, 21 has 2

## 1 has 0...
which(sapply(lst23,length)==0)
mm3 <- max(as.numeric(df$m3),na.rm=T) # 203
df$m3[df$m2==358] <- mm3+1 # no m3.

## 21 has 2
m2dup.i <- names(which(sapply(lst23,length)==2))

## when length(outlier) == 1...
m2dup.1 <- m2dup.i[sapply(lst23[m2dup.i],function(x){
	tmp <- any(x==1)
})] ## tricky: "62",11th => m3:18->32

## m2 -> NA: 1,2,3,4,6,7,9,10,
m2.mis.syms <- sapply(m2dup.1,function(m){
	this <- lst23[[m]]
	tmp <- names(which(this==1))
	(df %>% filter(is.na(m1) & m2==m & m3==tmp))$sym
})
non.null <- which(sapply(m2.mis.syms,length)==1)# c(1:4,6:7,9:10)
df$m2[df$sym %in% unlist(m2.mis.syms[non.null])] <- NA

lapply(m2dup.1,function(x){
	df %>% filter(m2==x)
})[!non.null]

m2.mis.syms <- sapply(m2dup.1,function(m){
	this <- lst23[[m]]
	tmp <- names(which(this==1))
	(df %>% filter(m2==m & m3==tmp))$sym
})[-non.null]

## m2 -> majoirty (guessed from m1): 12 (93->201)
df$m2[df$sym=="RAB1B"] <- 201

## m3 -> majority: 5 (67->26),8 (11->2),11 (18->32),
df$m3[df$sym=="PRMT1"] <- 26 
df$m3[df$sym=="LIN52"] <- 2
df$m3[df$sym=="CFDP1"] <- 32

##
lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..
table(sapply(lst23,length)) ## 395 has 1, 9 has 2

## 9 has 2...
m2dup.i <- names(which(sapply(lst23,length)==2))
m2dup.2 <- m2dup.i[sapply(lst23[m2dup.i],function(x){
	tmp <- !any(x==1)
})] ## tricky: "62",11th

tab2 <- lapply(m2dup.2,function(i){
	tmp <- df %>% filter(m2 == i)	
	tmp1 <- df %>% filter(m1 %in% tmp$m1 | m3 %in% tmp$m3)
	table(m1=tmp1$m1,m3=tmp1$m3)
})

lapply(m2dup.2,function(x){
	df %>% filter(m2==x)
})

#1,3,4,6,7,8,9 - split two (m3)
mm2 <- max(as.numeric(df$m2),na.rm=T) ## 404
for(i in c(1,3:4,6:9)){
	idx <- m2dup.2[i]
	tmp3 <- lst23[[idx]]
	tmp3i <- names(tmp3)[1]
	df$m2[which(df$m2==idx & df$m3==tmp3i)] <- mm2+i
}

#2,5 - most complicated: m2:1->0
mm3 <- max(as.numeric(df$m3),na.rm=T) ## 204
df$m3[which(df$m2=="119")] <- "150" ## 2: m3(34 => 150)
df$m3[which(df$m2=="18")] <- "15" ## 5: m3(1 => 15)

#########################
## 2 => 1 ("88" => "132")
#########################
lst12 <- tapply(df$m2,df$m1,table) # note, 0 is counted..
table(sapply(lst12,length)) ## 591 has 1, 8 has 2

## 8 has 2
m1dup.i <- names(which(sapply(lst12,length)==2))
m1dup.1 <- m1dup.i[sapply(lst12[m1dup.i],function(x){
	tmp <- any(x==1)
})] ## tricky: "62",11th => m2:18->22

lapply(m1dup.1,function(x){
	df %>% filter(m1==x)
})

## m2 -> NA: 1,2,3,4,6,7,9,10,
m1.mis.syms <- sapply(m1dup.1,function(m){
	this <- lst12[[m]]
	tmp <- names(which(this==1))
	(df %>% filter(m1==m & m2==tmp))$sym
})
df$m1[df$sym %in% m1.mis.syms] <- NA

##
lst12 <- tapply(df$m2,df$m1,table) # note, 0 is counted..
table(sapply(lst12,length)) ## 594 has 1, 5 has 2

## 5 has 2...
m1dup.i <- names(which(sapply(lst12,length)==2))
m1dup.2 <- m1dup.i[sapply(lst12[m1dup.i],function(x){
	tmp <- !any(x==1)
})] ## tricky: "62",11th

lapply(m1dup.2,function(x){
	df %>% filter(m1==x)
})

##
mm1 <- max(as.numeric(df$m1),na.rm=T) ## 599
for(i in 1:5){
	idx <- m1dup.2[i]
	tmp2 <- lst12[[idx]]
	tmp2i <- names(tmp2)[1]
	df$m1[which(df$m1==idx & df$m2==tmp2i)] <- mm1+i
}


lst12 <- tapply(df$m2,df$m1,table) 
lst23 <- tapply(df$m3,df$m2,table) 
table(sapply(lst12,length)) ## all ones 604
table(sapply(lst23,length)) ## all ones 411

if(0){
	setwd(rda_dir)
	saveRDS(df,file="members_prelim.rds")
}else{
	setwd(rda_dir)
	df <- readRDS(file="members_prelim.rds") # df
}

if(0){
	## filling from lo to hi
	#########################
	## 1 ("132")
	#########################
	th <- .5
	kept.0 <- df$l1 >= th ## this is redundant.
	df$m1[!kept.0] <- "0" # 35
	table(table(df$m1))
	df$m1[df$m1 %in% names(which(table(df$m1)==1))] <- "0" ## 'm1=="448"'' has only one member <- should be zero now.
	table(table(df$m1))
	sum(df$m1=="0") ## 157 genes have "0".

	## add back clusters at higher level when 0.
	#########################
	## 1 => 2 ("132" => "88")
	#########################
	# table(df$m1,df$m2)[1:5,1:5]

	# tab2 <- table(df$m2) # not now
	# ones2 <- names(which(tab2==1))
	# df$m2[df$m2 %in% ones2] <- "0"

	mm2 <- max(as.numeric(df$m2),na.rm=T)
	# kept.1 <- df$l2 >= th ## this is redundant.
	# df$m2[!kept.1] <- "0" # one more
	df$m2[is.na(df$m2)] <- "0"

	# table(table(df$m2))
	# df$m2[df$m2 %in% names(which(table(df$m2)==1))] <- "0"  # none
	# table(table(df$m2))

	tab12 <- tapply(df$m2,df$m1,function(x){
		y <- unique(x[!is.na(x)])
		if(length(y)==0)y <- NA
		return(y)
	})
	table(sapply(tab12,length)) ## 598 has 1, 4 has 2. 14 contain "0".
	which(sapply(tab12,length)==14)
	cl1.has.0 <- unique(df$m1[df$m2=="0"])
	cl1.has.0 <- cl1.has.0[cl1.has.0!="0"]

	for(i in seq(cl1.has.0)){
		df$m2[df$m1==cl1.has.0[i]] <- mm2+i
	}
	mm2 <- max(as.numeric(df$m2),na.rm=T)
	sum(df$m2=="0") # 135

	#########################
	## 2 => 3 ("88" => "54")
	#########################
	mm3 <- max(as.numeric(df$m3),na.rm=T)
	# kept.2 <- df$l3 >= th ## this is redundant.
	# df$m3[!kept.2] <- "0" # one more
	df$m3[is.na(df$m3)] <- "0"

	table(table(df$m3))
	df$m3[df$m3 %in% names(which(table(df$m3)==1))] <- "0"  # none
	table(table(df$m3))

	tab23 <- tapply(df$m3,df$m2,function(x){
		y <- unique(x[!is.na(x)])
		if(length(y)==0)y <- NA
		return(y)
	})
	table(sapply(tab23,length)) ## 405 has 1, 16 has 2, 36 contain "0".
	which(sapply(tab23,length)==36) ## 412 has 1, 8 has 2, 41 contain "0".

	cl2.has.0 <- unique(df$m2[df$m3=="0"])
	cl2.has.0 <- cl2.has.0[cl2.has.0!="0"]

	for(i in seq(cl2.has.0)){
		df$m3[df$m2==cl2.has.0[i]] <- mm3+i
	}
	mm3 <- max(as.numeric(df$m3),na.rm=T) # 298
	sum(df$m3=="0") # 34

	# df[df$m3=="0",]
}

## sort the number
df$m1[is.na(df$m1)] <- "0"
df$m2[is.na(df$m2)] <- "0"
df$m3[is.na(df$m3)] <- "0"

## cluster with one element => treated as zero (not in cluster) - can be done at the end.
m1 <- df$m1
old1 <- sort(table(m1),decreasing=T)# old
ones1 <- names(which(old1==1)) #none
m1[m1 %in% ones1] <- "0"
old1 <- sort(table(m1),decreasing=T)# old
is.0 <- which(names(old1)=="0")
old1 <- c(old1[-is.0],old1[is.0])

new1 <- c(paste0("S",as.character(c(seq(length(old1)-1)))),"NA")
names(new1) <- names(old1)

##
m2 <- df$m2
old2 <- sort(table(m2),decreasing=T)
ones2 <- names(which(old2==1))
m2[m2 %in% ones2] <- "0"
old2 <- sort(table(m2),decreasing=T)# old
is.0 <- which(names(old2)=="0")
old2 <- c(old2[-is.0],old2[is.0])

new2 <- c(paste0("M",as.character(c(seq(length(old2)-1)))),"NA")
names(new2) <- names(old2)

##
m3 <- df$m3
old3 <- sort(table(m3),decreasing=T)
ones3 <- names(which(old3==1))
m3[m3 %in% ones3] <- "0"
old3 <- sort(table(m3),decreasing=T)# old
is.0 <- which(names(old3)=="0")
old3 <- c(old3[-is.0],old3[is.0])

new3 <- c(paste0("L",as.character(c(seq(length(old3)-1)))),"NA")
names(new3) <- names(old3)

m1 <- factor(new1[m1],levels=new1)
m2 <- factor(new2[m2],levels=new2)
m3 <- factor(new3[m3],levels=new3)

table(m1)
table(m2)
table(m3)

mems <- df
mems$m1 <- m1
mems$m2 <- m2
mems$m3 <- m3

##
if(0){
	setwd(rda_dir)
	saveRDS(mems,file="hierarchy_cluster.rds")
}else{
	setwd(rda_dir)
	mems <- readRDS(file="hierarchy_cluster.rds")
}
##
if(0){
	grep("RAN",syms,value=T)
	x <- "RAN"
	unlist(mget(df$gid[df$m1==df$m1[df$sym==x]],org.Hs.egSYMBOL))
	unlist(mget(df$gid[df$m2==df$m2[df$sym==x]],org.Hs.egSYMBOL))
	unlist(mget(df$gid[df$m3==df$m3[df$sym==x]],org.Hs.egSYMBOL))

	df$l1[df$m1==df$m1[df$sym==x]]
	df$l2[df$m2==df$m2[df$sym==x]]
	df$l3[df$m3==df$m3[df$sym==x]]
}
