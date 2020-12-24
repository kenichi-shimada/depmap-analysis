outs <- list()
for (r in paste0("r",seq(0,1,.2))){
	setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3_1/")
	setwd(r)
	setwd("clue_multi_thres")
	txt <- system("ls -l|grep \"Sep 1[23456]\"",intern=T)
	outs[[r]] <- sort(as.numeric(sub(".+cons_([0-9]+)\\.rds","\\1",txt)))
}

missed <- sapply(1:6,function(i){
	x <- outs[[i]]
	idx <- which(!30:200 %in% x) + 171*(i-1)
	paste(idx,collapse=",")
})

query <- paste(missed[2:3],collapse=",")

171-sapply(outs,length)

plot(c(30,200),c(1,6),type="n")
for(i in 1:6){
j <- which(!30:200 %in% outs[[i]])
	points((30:200)[j],rep(i,length(j)),col=i,pch=20)
}

