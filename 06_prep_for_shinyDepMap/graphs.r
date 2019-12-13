setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")
load("depmap_initial_19q3_local_run.rda")

##
library(dplyr)
library(igraph)
library(RColorBrewer)

##
draw.graph <- function(g="BCL2L1"){
	m3 <- (df.0 %>% filter(sym==g))$mem3

	m1s <- as.numeric(factor((df.0 %>% filter(mem3==m3))$mem1))
	uniq.cols <- brewer.pal(8,"Pastel1")
	cols <- uniq.cols[m1s]

	g <- graphs$mem3[[m3]]
	edges <- g$edges
	nodes <- g$nodes

	net <- graph_from_data_frame(edges,vertices = nodes,directed=F)

	V(net)$color <- cols
	V(net)$frame.color <- NA
	V(net)$label.family <- "Helvetica"
	# V(net)$label.font <- "helvetica"
	E(net)$width <- E(net)$width
	E(net)$color <- "grey70"

	set.seed(1)
	x <- net %>% plot(layout=layout_nicely)
}
setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3/graphs")
pdf("L33_apoptosis_nfkb.pdf",width=5,height=5)
par(mar=c(0,0,0,0))
draw.graph(g="BCL2L1")
dev.off()

pdf("L78_selenoproteins.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
draw.graph(g="GPX4")
dev.off()

pdf("L137_keap1.pdf",width=4,height=4)
par(mar=c(0,0,0,0))
draw.graph(g="KEAP1")
dev.off()

##
kp <- unlist(mget(c("KEAP1","CUL3","TSC2","KCTD10"),org.Hs.egSYMBOL2EG))
coef[kp,kp]

(df.0 %>% filter(mem3=="L33") %>% arrange(mem2,mem1))
(df.0 %>% filter(mem3=="L78") %>% arrange(mem2,mem1))
