
setwd(rda_dir);setwd("r0/clue_multi_thres")
length((30:200)[paste0("cons_",30:200,".rds") %in% dir()]) # 91 -> 115

(30:200)[c(152,160)] # 181, 189
system("squeue -u ks389")


setwd(rda_dir);setwd("r0.2/clue_multi_thres")
length((30:200)[paste0("cons_",30:200,".rds") %in% dir()]) # 171
30:200 # 171

setwd(rda_dir);setwd("r0.4/clue_multi_thres")
length((30:200)[paste0("cons_",30:200,".rds") %in% dir()]) # 171
((30:200)[!paste0("cons_",30:200,".rds") %in% dir()])

171 + 171 + match(188,30:200) # 501
system("squeue -u ks389|grep 11180052_501")

setwd(rda_dir);setwd("r0.6/clue_multi_thres")
length((30:200)[paste0("cons_",30:200,".rds") %in% dir()]) # 171
(30:200)[!paste0("cons_",30:200,".rds") %in% dir()]
171*3 + match(173,30:200)

setwd(rda_dir);setwd("r0.8/clue_multi_thres")
length((30:200)[paste0("cons_",30:200,".rds") %in% dir()]) # 171
idx1 <- ((30:200)[!paste0("cons_",30:200,".rds") %in% dir()]) # 147
idx1 <- ((1:171)[!paste0("cons_",30:200,".rds") %in% dir()]) # 147
idx1 <- idx1 + 171 * 4

setwd(rda_dir);setwd("r1/clue_multi_thres")
length((30:200)[paste0("cons_",30:200,".rds") %in% dir()]) # 171
idx2 <- ((30:200)[!paste0("cons_",30:200,".rds") %in% dir()]) # 147
idx2 <- ((1:171)[!paste0("cons_",30:200,".rds") %in% dir()]) # 147
idx2 <- idx2 + 171 * 5
length((1:171)[paste0("cons_",30:200,".rds") %in% dir()]) # 129
paste((1:171)[!paste0("cons_",30:200,".rds") %in% dir()],collapse=",")

idx <- c(idx1,idx2)
diff(idx)
paste(c(idx1,idx2),collapse=",")

ps <- expand.grid(d=30:200,k=1:6)
ps[idx,]