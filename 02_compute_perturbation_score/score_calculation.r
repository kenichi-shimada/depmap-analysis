setwd(rda_dir)
load(file="eff_dep_19q3.rda") 

ed <- cbind(e=as.vector(eff.1),d=as.vector(dep.1))
i.e.wo.d <- which(!is.na(ed[,1]) & is.na(ed[,2]))
i.d.wo.e <- which(is.na(ed[,1]) & !is.na(ed[,2]))
i.w.ed <- which(!is.na(ed[,1]) & !is.na(ed[,2]))
i.missing <- which(is.na(ed[,1]) & is.na(ed[,2]))

e.wo.d <- ed[i.e.wo.d,1] ## 1345642
d.wo.e <- ed[i.d.wo.e,2] ## 19001

ed1 <- ed[i.w.ed,] ## 5334181

s <- svd(ed1,1,0)
f <- s$u %*% t(s$d) 
contrib <- s$d/sum(s$d)
proj <- ed1 %*% contrib

##
n <- 1000
i.min.max <- i.w.ed[unique(c(head(order(ed1[,1],decreasing=T),n),
		head(order(ed1[,1],decreasing=F),n),
		head(order(ed1[,2],decreasing=T),n),
		head(order(ed1[,2],decreasing=F),n)))] # 3933

set.seed(123)
row.i <- unique(c(sample(i.w.ed,1000000),i.min.max))

df <- data.frame(ed[row.i,])
names(df) <- c("e","d")

if(1){
	system.time({
		mod.1 <- loess(as.vector(dep.1) ~ as.vector(eff.1), data=df)
		pred.d <- predict(mod.1,newdata=e.wo.d)
		mod.2 <- loess(as.vector(eff.1) ~ as.vector(dep.1), data=df)
		pred.e <- predict(mod.2,newdata=d.wo.e)
	}) ## 1,000,000 -> 8506 seconds (110x)
	save(mod.1,pred.d,mod.2,pred.e,file="impute_missings_101119.rda")	
}else{
	setwd(rda_dir)
	x <-load("impute_missings_101119.rda")	
}

ed[i.e.wo.d,"d"] <- pred.d
ed[i.d.wo.e,"e"] <- pred.e

ed.imp <- ed
proj.all <- ed.imp %*% contrib
scores <- eff.1
scores[,] <- proj.all

if(0){
	## View missing values
	par(mfrow=c(1,2))
	plot(table(rowSums(is.na(scores))),main="# missed cell lines (of 423) per gene",las=2)
	plot(table(colSums(is.na(scores))),main="# missing genes (of 15847) per cell line",las=2)
}

setwd(rda_dir)
saveRDS(scores,file="scores_15847g_423c_101119.rds")
