plot.stepfun <- function (x, xval, xlim, ylim = range(c(y, Fn.kn)), xlab = "x", 
    ylab = "f(x)", main = NULL, add = FALSE, verticals = TRUE, 
    do.points = (n < 1000), pch = par("pch"), col = par("col"), 
    col.points = col, cex.points = par("cex"), col.hor = col, 
    col.vert = col, lty = par("lty"), lwd = par("lwd"), horiz=FALSE, ...) 
{
    if (!is.stepfun(x)) {
        if (is.numeric(x)) {
            sarg <- substitute(x)
            x <- ecdf(x)
            attr(x, "call") <- call("ecdf", sarg)
        }
        else stop("'plot.stepfun' called with wrong type of argument 'x'")
    }
    if (missing(main)) 
        main <- {
            cl <- attr(x, "call")
            deparse(if (!is.null(cl)) 
                cl
            else sys.call())
        }
    knF <- knots(x)
    xval <- if (missing(xval)) 
        knF
    else sort(xval)
    if (missing(xlim)) {
        rx <- range(xval)
        dr <- if (length(xval) > 1L) 
            max(0.08 * diff(rx), median(diff(xval)))
        else abs(xval)/16
        xlim <- rx + dr * c(-1, 1)
    }
    else dr <- diff(xlim)
    xval <- xval[xlim[1L] - dr <= xval & xval <= xlim[2L] + dr]
    ti <- c(xlim[1L] - dr, xval, xlim[2L] + dr)
    ti.l <- ti[-length(ti)]
    ti.r <- ti[-1L]
    y <- x(0.5 * (ti.l + ti.r))
    n <- length(y)
    Fn.kn <- x(xval)
    dev.hold()
    on.exit(dev.flush())
    if (horiz){
	    if (add) 
	        segments(y,ti.l, y, ti.r, col = col.hor, lty = lty, 
	            lwd = lwd, ...)
	    else {
	        if (missing(ylim)) 
	            ylim <- range(c(y, Fn.kn))
	        plot(NA, NA, type = "n", xlim = ylim, ylim = xlim, xlab = ylab, 
	            ylab = xlab, main = main, ...)
	        segments(y,ti.l, y, ti.r, col = col.hor, lty = lty, 
	            lwd = lwd)
	    }
	    if (do.points) 
	        points(Fn.kn, xval, pch = pch, col = col.points, cex = cex.points)
	    if (verticals) 
	        segments(y[-n], xval, y[-1L], xval, col = col.vert, lty = lty, 
	            lwd = lwd)
	    invisible(list(t = ti, y = y))
    }else{
	    if (add) 
	        segments(ti.l, y, ti.r, y, col = col.hor, lty = lty, 
	            lwd = lwd, ...)
	    else {
	        if (missing(ylim)) 
	            ylim <- range(c(y, Fn.kn))
	        plot(NA, NA, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, 
	            ylab = ylab, main = main, ...)
	        segments(ti.l, y, ti.r, y, col = col.hor, lty = lty, 
	            lwd = lwd)
	    }
	    if (do.points) 
	        points(xval, Fn.kn, pch = pch, col = col.points, cex = cex.points)
	    if (verticals) 
	        segments(xval, y[-n], xval, y[-1L], col = col.vert, lty = lty, 
	            lwd = lwd)
	    invisible(list(t = ti, y = y))
    }
}
