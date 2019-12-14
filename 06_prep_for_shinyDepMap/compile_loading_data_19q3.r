library(org.Hs.eg.db)
library(dplyr)
library(tibble)

setwd(rda_dir)

## perturbation scores
p.scores <- readRDS("scores_15847g_423c_101119.rds")
p.scores <- round(p.scores,4)

used.cells <- colnames(p.scores)

## cell line info
cl.info <- readRDS("celllines_19q3.rds")
names(cl.info)[3] <- "CCLE_Name"
cl.info <- data.frame(cl.info,ind=as.numeric(sub("ACH-0","",cl.info$DepMap_ID)))
cl.info <- cl.info[match(used.cells,cl.info$DepMap_ID),] ## cl.info done

## essential gene ids (names are symbols)
eids.ess <- readRDS("ess_2492_genes.rds")
ess.syms <- unlist(mget(eids.ess,org.Hs.egSYMBOL))
names(eids.ess) <- paste0(ess.syms," (",eids.ess,")")## eids.ess done

## efficacy vs selectivity
load("eff-sel.rda") # ef.sel, thres
names(ef.sel) <- c("efficacy","selectivity")
es <- round(ef.sel,3)

eids.all <- rownames(ef.sel)
is.ess <- eids.all %in% eids.ess
syms.all <- unlist(mget(eids.all,org.Hs.egSYMBOL,ifnotfound=NA))
syms.all.1 <- paste0(syms.all," (",eids.all,")")## eids.ess done

## cluster membership
mems <- readRDS("hierarchy_cluster.rds")[-2] # symbol is removed
rownames(mems) <- c()
colnames(mems) <- c("eid","mem1","mem2","mem3")

## compute density on efficacy-selectivity plot
smoothScatterCalcDensity <-	function (x, nbin, bandwidth, range.x) {
    if (length(nbin) == 1) 
        nbin <- c(nbin, nbin)
    if (!is.numeric(nbin) || length(nbin) != 2) 
        stop("'nbin' must be numeric of length 1 or 2")
    if (missing(bandwidth)) {
        bandwidth <- diff(apply(x, 2, stats::quantile, probs = c(0.05, 
            0.95), na.rm = TRUE, names = FALSE))/25
        bandwidth[bandwidth == 0] <- 1
    }
    else {
        if (!is.numeric(bandwidth)) 
            stop("'bandwidth' must be numeric")
        if (any(bandwidth <= 0)) 
            stop("'bandwidth' must be positive")
    }
    rv <- KernSmooth::bkde2D(x, bandwidth = bandwidth, gridsize = nbin, 
        range.x = range.x)
    rv$bandwidth <- bandwidth
    rv
}

xy <- xy.coords(ef.sel, setLab = FALSE)
select <- is.finite(xy$x) & is.finite(xy$y)
x <- cbind(xy$x, xy$y)[select, ]
map <- smoothScatterCalcDensity(x, nbin=1024)
mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
dens <- map$fhat[cbind(xbin, ybin)]
dens[is.na(dens)] <- 0


##
df <- data.frame(eid=eids.all,efficacy=es[,1],selectivity=es[,2],
	is.ess=is.ess,sym=syms.all,sym.all=syms.all.1,stringsAsFactors=F) %>%
    mutate(dens=round(dens,3)) %>%
    mutate(text = syms.all)

## rtsne coordinates
setwd(rda_dir);setwd("rtsne")
rt <- readRDS(file="rtsne_1.rds") %>%`colnames<-`(c("x","y"))
rt <- round(rt,3)

## df.0, subsetted data for essential genes only
df.0 <- df %>% filter(is.ess) %>% 
    left_join(mems,by="eid") %>%
    mutate(x=rt[,1],y=rt[,2]) %>%
    mutate(text=paste0(sym,
        "<br />eff:",efficacy,", sel:",selectivity,
        "<br />clust:",mem1,",",mem2,",",mem3))

##
thres.eff <- thres$eff

##
dummy <- data.frame(efficacy=c(-1,-.9),selectivity=c(0,0.1),x=c(0,0.1),y=c(0,0.1),text="")
cs <- c("Small","Medium","Large")

if(1){
    setwd(rda_dir)
    save(cl.info,df,df.0,eids.ess,thres.eff,dummy,cs,
    	file="depmap_initial_19q3.rda")
}else{
    setwd(rda_dir)
    load(file="depmap_initial_19q3.rda")
}

##
setwd(rda_dir)
coef <- readRDS("eff_coef_2492_genes.rds")
coef <- round(coef,4)

##
m1 <- levels(df.0$mem1)
m2 <- levels(df.0$mem2)
m3 <- levels(df.0$mem3)

##
levs <- m1

this.mem <- "mem1"
m1s <- lapply(m1,function(this.cl){
    cat("*")
    df.1 <- df.0[df.0[[this.mem]] %in% this.cl,] %>%
        arrange(mem3,mem2,mem1)
    eids <- df.1$eid
    syms <- df.1$sym
    names(syms) <- eids

    coefs <- coef[eids,eids]
    th <- .1
    idx <- which(coefs > th & lower.tri(coefs),arr.ind=T)

    lev1 <- levs[levs %in% unique(df.1$mem1)]
    nodes <- df.0[df.0[[this.mem]] %in% this.cl,] %>%
        mutate(mem1 = mem1,exclude=NULL) %>%
        mutate(mem1=factor(mem1,levels=lev1,exclude=NULL)) %>%
        mutate(title = sub("<br />.+<br />clust:","<br />",text)) %>%
        dplyr::select(eid,sym,mem1,title) %>%
        `colnames<-`(c("id","label","group","title")) 

    edges <- data.frame(from=eids[idx[,"row"]],to=eids[idx[,"col"]],
        stringsAsFactors=F) %>%
        mutate(coef=coefs[idx])%>%
        mutate(width=round(coef*10)) %>%
        mutate(title=paste0(syms[from],"-",syms[to],"<br />",coef)) %>%
        dplyr::select(from,to,width,title)
    return(list(nodes=nodes,edges=edges))
})

this.mem <- "mem2"
m2s <- lapply(m2,function(this.cl){
    cat("*")
    df.1 <- df.0[df.0[[this.mem]] %in% this.cl,] %>%
        arrange(mem3,mem2,mem1)
    eids <- df.1$eid
    syms <- df.1$sym
    names(syms) <- eids

    coefs <- coef[eids,eids]
    th <- .1
    idx <- which(coefs > th & lower.tri(coefs),arr.ind=T)

    lev1 <- levs[levs %in% unique(df.1$mem1)]
    nodes <- df.0[df.0[[this.mem]] %in% this.cl,] %>%
        mutate(mem1 = mem1,exclude=NULL) %>%
        mutate(mem1=factor(mem1,levels=lev1,exclude=NULL)) %>%
        mutate(title = sub("<br />.+<br />clust:","<br />",text)) %>%
        dplyr::select(eid,sym,mem1,title) %>%
        `colnames<-`(c("id","label","group","title")) 

    edges <- data.frame(from=eids[idx[,"row"]],to=eids[idx[,"col"]],
        stringsAsFactors=F) %>%
        mutate(coef=coefs[idx])%>%
        mutate(width=round(coef*10)) %>%
        mutate(title=paste0(syms[from],"-",syms[to],"<br />",coef)) %>%
        dplyr::select(from,to,width,title)
    return(list(nodes=nodes,edges=edges))
})

this.mem <- "mem3"
m3s <- lapply(m3,function(this.cl){
    cat("*")
    df.1 <- df.0[df.0[[this.mem]] %in% this.cl,] %>%
        arrange(mem3,mem2,mem1)
    eids <- df.1$eid
    syms <- df.1$sym
    names(syms) <- eids

    coefs <- coef[eids,eids]
    th <- .1
    idx <- which(coefs > th & lower.tri(coefs),arr.ind=T)

    lev1 <- levs[levs %in% unique(df.1$mem1)]
    nodes <- df.0[df.0[[this.mem]] %in% this.cl,] %>%
        mutate(mem1 = mem1,exclude=NULL) %>%
        mutate(mem1=factor(mem1,levels=lev1,exclude=NULL)) %>%
        mutate(title = sub("<br />.+<br />clust:","<br />",text)) %>%
        dplyr::select(eid,sym,mem1,title) %>%
        `colnames<-`(c("id","label","group","title")) 

    edges <- data.frame(from=eids[idx[,"row"]],to=eids[idx[,"col"]],
        stringsAsFactors=F) %>%
        mutate(coef=coefs[idx])%>%
        mutate(width=round(coef*10)) %>%
        mutate(title=paste0(syms[from],"-",syms[to],"<br />",coef)) %>%
        dplyr::select(from,to,width,title)
    return(list(nodes=nodes,edges=edges))
})

names(m1s) <- m1
names(m2s) <- m2
names(m3s) <- m3
m1s <- m1s[!grepl("NA",m1)]
m2s <- m2s[!grepl("NA",m2)]
m3s <- m3s[!grepl("NA",m3)]

graphs <- list(mem1=m1s,mem2=m2s,mem3=m3s)

setwd(rda_dir)
save(cl.info,df,df.0,eids.ess,thres.eff,dummy,cs,
    p.scores,coef,graphs,
    file="depmap_initial_19q3_local_run.rda")
