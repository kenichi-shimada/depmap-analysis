## derived from dbscan::hullplot
my.hullplot <- function (x, cl, col = NULL, cex = 0.5, hull_lwd = 1, hull_lty = 1, 
    solid = TRUE, alpha = 0.2, main = "Convex Cluster Hulls", pch=20,
    ...) 
{
    if (ncol(x) > 2) 
        x <- prcomp(x)$x
    if (is(cl, "xics") || "clusters_xi" %in% names(cl)) {
        clusters_xi <- cl$clusters_xi
        cl_order <- cl$order
    }
    else clusters_xi <- NULL
    if (is.list(cl)) 
        cl <- cl$cluster
    if (!is.numeric(cl)) 
        stop("Could not get cluster assignment vector from cl.")
    if (is.null(col)) 
        col <- palette()
    if (max(cl) + 1L > length(col)) 
        warning("Not enough colors. Some colors will be reused.")
    idx <- tapply(seq(cl),cl,identity)
    # return(idx)
    plot(x[idx[[1]], 1:2], col = col[min(cl)%%length(col) + 1L], cex = cex, 
        main = main, xlim=range(x[,1]),ylim=range(x[,2]),pch=pch,...)
    sapply(names(idx)[-1],function(i){
		points(x[idx[[i]],1:2],col=col[as.numeric(i)%%length(col) + 1L],cex=.8,pch=pch)    	
    })
    col_poly <- adjustcolor(col, alpha.f = alpha)
    border <- col
    if (is.null(hull_lwd) || is.na(hull_lwd) || hull_lwd == 0) {
        hull_lwd <- 1
        border <- NA
    }
    if (is(cl, "xics") || "clusters_xi" %in% names(cl)) {
        clusters_xi <- clusters_xi[order(-(clusters_xi$end - 
            clusters_xi$start)), ]
        ci_order <- clusters_xi$cluster_id
    }
    else {
        ci_order <- 1:max(cl)
    }
    for (i in 1:length(ci_order)) {
        if (is.null(clusters_xi)) {
            d <- x[cl == i, ]
        }
        else {
            d <- x[cl_order[clusters_xi$start[i]:clusters_xi$end[i]], 
                ]
        }
        ch <- chull(d)
        ch <- c(ch, ch[1])
        if (!solid) {
            lines(d[ch, ], col = col[ci_order[i]%%length(col) + 
                1L], lwd = hull_lwd, lty = hull_lty)
        }
        else {
            polygon(d[ch, ], col = col_poly[ci_order[i]%%length(col_poly) + 
                1L], lwd = hull_lwd, lty = hull_lty, border = border[ci_order[i]%%length(col_poly) + 
                1L])
        }
    }
}
