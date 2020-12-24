setwd(rda_dir)
x <- load("depmap_initial_19q3_local_run_v3_4.rda")

##
library(dplyr)
library(igraph)
library(RColorBrewer)

##
draw.graph <- function(g="BCL2L1",clust=c("S","M","L"),seed=123){
	if(clust=="S"){
		m1 <- (mems %>% filter(sym==g))$mem1
		m1s <- 1
		g <- graphs$mem1[[m1]]
	}else if(clust=="M"){
		m2 <- (mems %>% filter(sym==g))$mem2
		m1s <- as.numeric(factor((mems %>% filter(mem2==m2))$mem1))
		uniq.cols <- brewer.pal(8,"Pastel1")
		cols <- uniq.cols[m1s]
		g <- graphs$mem2[[m2]]
	}else if(clust=="L"){
		m3 <- (mems %>% filter(sym==g))$mem3
		m1s <- as.numeric(factor((mems %>% filter(mem3==m3))$mem1))
		uniq.cols <- brewer.pal(8,"Pastel1")
		cols <- uniq.cols[m1s]
		g <- graphs$mem3[[m3]]
	}

	uniq.cols <- brewer.pal(8,"Pastel1")
	cols <- uniq.cols[m1s]

	edges <- g$edges
	nodes <- g$nodes

	net <- graph_from_data_frame(edges,vertices = nodes,directed=F)

	V(net)$color <- cols
	V(net)$frame.color <- NA
	V(net)$label.family <- "Helvetica"
	# V(net)$label.font <- "helvetica"
	E(net)$width <- E(net)$width
	E(net)$color <- "grey70"

	set.seed(seed)
	x <- net %>% plot(layout=layout_nicely)
}

setwd(plot_dir);setwd("07_echodots_prep")

pdf("S152_cul3.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
draw.graph(g="CUL3",clust="S",seed=1)
dev.off()

pdf("L91_selenoproteins.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
draw.graph(g="GPX4",clust="L",seed=1)
dev.off()

pdf("L119_kras.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
draw.graph(g="KRAS",clust="L",seed=2)
dev.off()

##
kp <- unlist(mget(c("KEAP1","CUL3","KCTD10"),org.Hs.egSYMBOL2EG))
coef.kp <- coef[kp,kp]
colnames(coef.kp) <- rownames(coef.kp) <- c("KEAP1","CUL3","KCTD10")

mget(c("KEAP1","CUL3","KCTD10"),org.Hs.egSYMBOL2EG)
(df.0 %>% filter(mem3=="L33") %>% arrange(mem2,mem1))
(df.0 %>% filter(mem3=="L78") %>% arrange(mem2,mem1))
mems %>% filter(mem3=="L91")