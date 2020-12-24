library(org.Hs.eg.db)
library(dplyr)
library(tibble)

# setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/rdas/19q3")

## density
smoothScatterCalcDensity <- function (x, nbin, bandwidth, range.x) {
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

mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2

## cell line info
setwd(rda_dir)
load(file="find_essential_genes.rda") # is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns
cl.info <- readRDS("celllines_19q3.rds")
names(cl.info)[3] <- "CCLE_Name"
cl.info <- data.frame(cl.info,ind=as.numeric(sub("ACH-0","",cl.info$DepMap_ID)))

p.scores <- readRDS("r0/scores_15847g_423c.rds")
used.cells <- colnames(p.scores)

cl.info <- cl.info[match(used.cells,cl.info$DepMap_ID),] ## cl.info done

x <- load("eff-sel_6thres.rda") # ef.sel, thres

##
for(k in 1:6){ # mix.ratio between CRISPR and shRNA
    c1 <- seq(0,1,.2)[k]

    setwd(rda_dir);setwd(paste0("r",c1));
    x <- load(file="ecodots_ori_v3.rda")   # mems,liks,ds,ratios,optims,types

    ## perturbation scores
    p.scores <- readRDS("scores_15847g_423c.rds")
    p.scores <- round(p.scores,4)

    used.cells <- colnames(p.scores)

    eids.ess <- names(which(is.ess[[k]]$q1))
    ess.syms <- unlist(mget(eids.ess,org.Hs.egSYMBOL))

    names(eids.ess) <- paste0(ess.syms," (",eids.ess,")")## eids.ess done

    ## efficacy vs selectivity
    ef.sel <- ef.sels[[k]]
    ef.sel <- lapply(ef.sel,function(es){
        names(es) <- c("efficacy","selectivity")
        es <- round(es,3)
        return(es)
    })

    eids.all <- rownames(ef.sel[[1]])

    is.ess.1 <- is.ess[[k]]$q1
    syms.all <- unlist(mget(eids.all,org.Hs.egSYMBOL,ifnotfound=NA))
    syms.all.1 <- paste0(syms.all," (",eids.all,")")## eids.ess done

    dfs <- lapply(ef.sel,function(es){
        xy <- xy.coords(es, setLab = FALSE)
        select <- is.finite(xy$x) & is.finite(xy$y)
        x <- cbind(xy$x, xy$y)[select, ]
        map <- smoothScatterCalcDensity(x, nbin=1024)

        xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
        ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
        dens <- map$fhat[cbind(xbin, ybin)]
        dens[is.na(dens)] <- 0
        if(0){
            ## confirmation that density is computed properly
            dens.thres <- 0.1
            plot(es,col=(dens > dens.thres)+1,pch=20)
            smoothScatter(es)
            hist(dens,breaks=1000)
            abline(v=dens.thres,col=2)
            sum(dens > dens.thres)
        }

        ##
        df <- data.frame(eid=eids.all,efficacy=es[,1],selectivity=es[,2],
            is.ess=is.ess.1,sym=syms.all,sym.all=syms.all.1,stringsAsFactors=F) %>%
            mutate(dens=round(dens,3)) %>%
            mutate(text = syms.all)
        return(df)
    })

    # head(dfs[[1]])
    # head(dfs[[2]])
    # head(dfs[[3]])

    ## rtsne coordinates
    setwd(rda_dir);setwd(paste0("r",c1));
    setwd("rtsne_multi_thres_1")
    rt <- readRDS(file="rtsne_ori_1.rds") %>%
        `colnames<-`(c("x","y"))
    rt <- round(rt,3)

    ## cluster membership
    setwd(rda_dir);setwd(paste0("r",c1));
    mems <- readRDS("ori_hierarchy_cluster_v3.rds")
    rownames(mems) <- c()
    colnames(mems) <- c("eid","sym","mem1","mem2","mem3","lik1","lik2","lik3")
    mems <- mems %>%
        mutate(lik1=round(lik1,3)) %>%
        mutate(lik2=round(lik2,3)) %>%
        mutate(lik3=round(lik3,3)) %>%
        mutate(x=rt[,1],y=rt[,2])

    if(0){
        ## df.0, subsetted data for essential genes only
        df.0s <- lapply(dfs,function(df){
            df.0 <- df %>% filter(is.ess) %>% 
                left_join(mems,by="eid") %>%
                mutate(text=paste0(sym,
                    "<br />eff:",efficacy,", sel:",selectivity,
                    "<br />clust:",mem1,",",mem2,",",mem3))
            return(df.0)
        })
    }

    ##
    thres.eff <- sc.ths[[paste0("r",c1)]][["0.001"]]

    ##
    dummy <- data.frame(efficacy=c(-1,-.9),selectivity=c(0,0.1),x=c(0,0.1),y=c(0,0.1),text="")
    cs <- c("Small","Medium","Large")

    ##
    setwd(rda_dir); setwd(paste0("r",c1))
    coef <- readRDS("ess_coef.rds")
    coef <- round(coef,4)

    ##
    m1 <- levels(mems$mem1)
    m2 <- levels(mems$mem2)
    m3 <- levels(mems$mem3)

    ##
    levs <- m1

    this.mem <- "mem1"
    m1s <- lapply(m1,function(this.cl){
        cat("*")
        df.1 <- mems[mems[[this.mem]] %in% this.cl,] %>%
            arrange(mem3,mem2,mem1)
        eids <- df.1$eid
        syms <- df.1$sym
        names(syms) <- eids

        coefs <- coef[eids,eids]
        th <- .1
        idx <- which(coefs > th & lower.tri(coefs),arr.ind=T)

        lev1 <- levs[levs %in% unique(df.1$mem1)]
        nodes <- mems[mems[[this.mem]] %in% this.cl,]%>%
            mutate(mem1 = mem1,exclude=NULL) %>%
            mutate(mem1=factor(mem1,levels=lev1,exclude=NULL)) %>%
            mutate(title=paste0(sym,"<br />",mem1,",",mem2,",",mem3)) %>%
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
        df.1 <- mems[mems[[this.mem]] %in% this.cl,] %>%
            arrange(mem3,mem2,mem1)
        eids <- df.1$eid
        syms <- df.1$sym
        names(syms) <- eids

        coefs <- coef[eids,eids]
        th <- .1
        idx <- which(coefs > th & lower.tri(coefs),arr.ind=T)

        lev1 <- levs[levs %in% unique(df.1$mem1)]
        nodes <- mems[mems[[this.mem]] %in% this.cl,] %>%
            mutate(mem1 = mem1,exclude=NULL) %>%
            mutate(mem1=factor(mem1,levels=lev1,exclude=NULL)) %>%
            mutate(title=paste0(sym,"<br />",mem1,",",mem2,",",mem3)) %>%
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
        df.1 <- mems[mems[[this.mem]] %in% this.cl,] %>%
            arrange(mem3,mem2,mem1)
        eids <- df.1$eid
        syms <- df.1$sym
        names(syms) <- eids

        coefs <- coef[eids,eids]
        th <- .1
        idx <- which(coefs > th & lower.tri(coefs),arr.ind=T)

        lev1 <- levs[levs %in% unique(df.1$mem1)]
        nodes <- mems[mems[[this.mem]] %in% this.cl,] %>%
            mutate(mem1 = mem1,exclude=NULL) %>%
            mutate(mem1=factor(mem1,levels=lev1,exclude=NULL)) %>%
            mutate(title=paste0(sym,"<br />",mem1,",",mem2,",",mem3)) %>%
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

    # m1s <- m1s[!grepl("NA",m1)]
    # m2s <- m2s[!grepl("NA",m2)]
    # m3s <- m3s[!grepl("NA",m3)]

    graphs <- list(mem1=m1s,mem2=m2s,mem3=m3s)

    setwd(rda_dir);# setwd(paste0("r",c1))
    save(cl.info,dfs,mems,eids.ess,thres.eff,dummy,cs,
        file=paste0("depmap_initial_19q3_v3_",k,".rda"))

    save(cl.info,dfs,mems,eids.ess,thres.eff,dummy,cs,p.scores,coef,graphs,
        file=paste0("depmap_initial_19q3_local_run_v3_",k,".rda"))
}
## compile

# common: cl.info, dummy,cs
# mix_ratio: mems, eids.ess, thres.eff
# mix_ratio, threshod: dfs
dfs,
mems,
eids.ess,
thres.eff,

# p.scores,
# coef,
# graphs,

eids.ess.all <- thres.eff.all <- mems.all <- dfs.com <- dfs.all <- list()
for(mix.ratio in 1:6){
    setwd(rda_dir)
    file <- paste0("depmap_initial_19q3_local_run_v3_",mix.ratio,".rda")
    x <- load(file)

    dfs.all[[mix.ratio]] <- lapply(dfs,function(x)x[c(1:4,7)])
    dfs.com[[mix.ratio]] <- lapply(dfs,function(x)x[c(1,5:6)])
    mems.all[[mix.ratio]] <- mems

    thres.eff.all[[mix.ratio]] <- thres.eff
    eids.ess.all[[mix.ratio]] <- eids.ess
}

## dfs
com <- dfs.com[[1]][[1]]
sapply(dfs.com,sapply,function(y)identical(com,y)) #com are all identical
dfs <- dfs.all

dfs <- lapply(dfs,function(x){
    lapply(x,function(y){
        names(y)[4] <- "is.ess"
        return(y)
    })
})

## mems
mems <- mems.all

## thres.eff
thres.eff <- lapply(thres.eff.all,function(x)round(x,3))

## eids.ess
eids.ess <- eids.ess.all

## 
names(mems) <- names(eids.ess) <- names(thres.eff) <- names(dfs) <- 
    rev(c("CRISPR","80:20","60:40 (default)","40:60","20:80","shRNA"))
lapply(list(cl.info,dummy,cs,mems,eids.ess,thres.eff,dfs,com),function(x)format(object.size(x),units="MiB",standard="IEC"))

##
if(0){
    setwd(rda_dir)
    save(cl.info,dummy,cs,mems,eids.ess,thres.eff,dfs,com,file="depmap_initial_19q3_v3.rda")
}else{
    setwd(rda_dir)
    x <- load("depmap_initial_19q3_v3.rda")
}
n1 <- names(dfs)
n2 <- names(dfs[[1]])

nlayer <- 31
n.dthres <- 3000

if(0){
    library(ggplot2)
    library(plotly)
    library(dplyr)

    all.scatter.plots <- list()
    for(i in n1){
        for(j in n2){
            df <- dfs[[i]][[j]] %>% left_join(com,by="eid")
            dthres <- sort(dfs[[i]][[j]]$dens)[n.dthres]

            p <- ggplot(df,aes(efficacy,selectivity)) +
            geom_vline(xintercept=c(0,thres.eff[[i]][j]), col=c("grey50","grey70")) +
            geom_hline(yintercept=0,col="grey50") +
            stat_density_2d(aes(fill = stat(level)),geom="polygon",bins=nlayer) +
            scale_fill_gradientn(colors=rev(blues9)) + 
            geom_point(data=(df %>% filter(dens <= dthres)),
              aes(efficacy, selectivity,text=sym), 
              color=blues9[9], size = 1) +
            theme_bw() +
              theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"lines"),
                axis.text=element_text(size=10),axis.title=element_text(size=10),
                plot.title=element_text(size=10),
                axis.text.x = element_text(angle = 0)) +
            geom_point(data=dummy[1,],aes(efficacy,selectivity),col=NA,size=1) + # sym.query
            geom_point(data=dummy[2,],aes(efficacy,selectivity),col=NA,size=1) # sym1

            ply <- plotly::ggplotly(p,tooltip="text")
            ndata <- length(ply$x$data)
            # sapply(ply$x$data,function(x)length(x$text))  
            for(k in seq(nlayer+3)){
              ply$x$data[[k]]$hoverinfo <- "none"
            }
            all.scatter.plots[[i]][[j]] <- ply
        }
    }

    setwd(rda_dir)
    saveRDS(all.scatter.plots,file="depmap_all_scatter_plots_v3.rds")
}

##
if(0){
    x <- load(file="depmap_initial_19q3_v3.rda")
    # lapply(x,function(y)format(object.size(get(y)),units="MiB",standard="IEC"))

    # setwd("/n/groups/mitchison/Kenichi/projects/09 depmap/plots/19q3_1/adam2")
    # # crossoverpoints <- readRDS(file="crossoverpoints.rds")
    # crossoverpoints <- readRDS(file="crossoverpoints_bagel.rds")

    # x <- load(file="essentiality_matrix.rda") # is.ess,eff.ths

    setwd(rda_dir)
    save(cl.info,dummy,cs,mems,eids.ess,thres.eff,dfs,com,
        eff.ths,#crossoverpoints,
        file="depmap_initial_19q3_v3.rda")
}
