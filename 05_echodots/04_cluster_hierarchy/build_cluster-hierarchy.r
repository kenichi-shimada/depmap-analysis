library(dplyr)
library(org.Hs.eg.db)

setwd(rda_dir)
load(file="ecodots.rda")	 # mems,liks,ds,ratios,optims,types
ess.genes <- readRDS(file="ess_2492_genes.rds")
syms <- unlist(mget(ess.genes,org.Hs.egSYMBOL))

ds <- rev(readRDS("three_ds.rds"))

df <- data.frame(
	gid=ess.genes,sym=syms,
	m1=as.character(mems[[ds[1]]]),
	m2=as.character(mems[[ds[2]]]),
	m3=as.character(mems[[ds[3]]]),
	stringsAsFactors=F)

sum(is.na(df$m1)) # 119
sum(is.na(df$m2)) # 15
sum(is.na(df$m3)) # 7

## task: splitting clusters if all the cluster members at lower level don't correspond to the same cluster at higher level.

## 1. one cluster (lo) => one cluster (hi), if not, split the lower cluster, because of the scheme, should do 2=>3, and 1=>2
## 2. one cluster (lo) => none (hi) is not taken care of at this point

## filling from hi to lo one lo should lead to one hi
#########################
## 3 => 2 ("53" => "87")
#########################
lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..
table(sapply(lst23,length)) ## 1 has 0, 382 has 1, 21 has 2

## 1 has 0...
m2.wo.m3 <- names(which(sapply(lst23,length)==0))
df %>% filter(m2==m2.wo.m3)
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
df %>% filter(m2=="93")
df$m2[df$sym=="RAB1B"] <- 201

## m3 -> majority: 5 (67->26),8 (11->2),11 (18->32),
df %>% filter(m2=="227")
df$m3[df$sym=="PRMT1"] <- 26 
df %>% filter(m2=="297")
df$m3[df$sym=="LIN52"] <- 2
df %>% filter(m2=="62")
df$m3[df$sym=="CFDP1"] <- 32

##
lst23 <- tapply(df$m3,df$m2,table) # note, 0 is counted..
table(sapply(lst23,length)) ## 395 has 1, 9 has 2

## 9 has 2...
m2dup.i <- names(which(sapply(lst23,length)==2))
m2dup.2 <- m2dup.i[sapply(lst23[m2dup.i],function(x){
	tmp <- !any(x==1)
})] 

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
## 2 => 1 ("87" => "131")
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
table(sapply(lst12,length)) ## 595 has 1, 4 has 2

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



## rename cluster labels after they are sorted in the decreasing order.
## missing values are set as "0"

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
setwd(rda_dir)
saveRDS(mems,file="hierarchy_cluster.rds")
write.csv(mems,"SD3_cluster_hierarchy.csv")
