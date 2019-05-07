#' \code{manhattan plot}
#'
#'
#' @param inputdf An input data.frame to plot.
#'                If a cex weight column "cw" presents, dots will be bigger for larger value.
#'                plot_mht: [df, =col2plot, "chr", "pos", "cw"(optional)]
#' @param hdf A data.frame contains regions to be highlighted in the manhanttan plot. [df, =NULL, required columns: "chr", "pos"]
#' @param col2plot The column name to plot. [chr, ="ModelFreq"]
#' @param jmph Plot joint mht. [logical, =FALSE], if =TURE, must have cex weight column (cw) in inputdf: [,cw].
#' @param pl Plot lines instead of dots. [logical, =FALSE]
#' @param cl_file The full path of the chromosome length file. [chr, ="~/bin/zmSNPtools/Rcodes/chr_length_B73v3.csv"]
#' @param GAP The physical gap value between chrs. [num, =5e+06]
#' @param ... Pass parameters to plot.
#'
#' @export
plot_mht <- function(inputdf, hdf=NULL, col2plot="ModelFreq", jmph=FALSE, pl=FALSE,
                     cl_file = "~/bin/zmSNPtools/Rcodes/chr_length_B73v3.csv",
                     cex=.9, pch=16, col=rep(c("slateblue", "cyan4"), 5),
                     GAP=5e+06, yaxis=NULL,
                      ... ){

    res <- inputdf
    cl <- read.csv(cl_file)

  res <- newpos(d=res, GAP, cl)
  chrtick <- chrline_tick(GAP, cl)


  #### setup the cavon
  if(is.null(yaxis)){
    plot(x=-1000, y=-1000,  type="p", xaxt="n", xlab="",
         xlim=c(0, max(chrtick$chrlines)), ylim=c(0, max(res[, col2plot], na.rm=TRUE)*1.3 ),
         ...)
  }else{
    plot(x=-1000, y=-1000,  type="p", xaxt="n", yaxt="n", xlab="",
         xlim=c(0, max(chrtick$chrlines)),
         ...)
    axis(side=2, at=NULL, labels=yaxis)
  }
  axis(side=1, at=chrtick$ticks, labels=c("chr1", "chr2", "chr3", "chr4", "chr5",
                                          "chr6", "chr7", "chr8", "chr9", "chr10"))
  abline(v=chrtick$chrlines, col="grey")

  for(i in 1:10){
      if(jmph){
          points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
                 pch = pch, col=col[i], cex=res[res$chr==i, ]$cw*cex);
      }else{
          if(pl){
              segments(x0=subset(res, chr==i)$newpos, y0=rep(0, nrow(res[res$chr==i,])),
                        x1=subset(res, chr==i)$newpos, y1=res[res$chr==i, col2plot],
                    col=col[i])
          }else{
            if(sum(names(inputdf) %in% "cw") == 1){
              points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
                     pch = pch, col=col[i], cex=res[res$chr==i, ]$cw*cex);
            }else{
              points(x=subset(res, chr==i)$newpos, y=res[res$chr==i, col2plot],
                     pch = pch, col=col[i], cex=cex);
            }
          }
      }
  }
  ### add gene annotation information
  if(!is.null(hdf)){
    hdf <- newpos(d=hdf, GAP, cl)
    abline(v=hdf$newpos, col="red", lty=2, lwd=0.6)
  }

}


#' @rdname plot_mht
#' @export
newpos <- function (d, GAP, cl)
{
  if (!("chr" %in% names(d) & "pos" %in% names(d))){
    stop("Make sure your data frame contains columns chr and pos")
  }

  cl$accumpos <- cl$BP
  cl <- cl[order(cl$CHR), ]
  d$newpos <- d$pos;
  for (i in 2:10) {
    cl[cl$CHR == i, ]$accumpos <- cl[cl$CHR == (i - 1), ]$accumpos + cl[cl$CHR == i, ]$accumpos + GAP
    d[d$chr == i, ]$newpos <- d[d$chr == i, ]$pos + cl[cl$CHR == (i - 1), ]$accumpos + GAP
  }
  return(d)
}

#' @rdname plot_mht
#' @export
chrline_tick <- function(GAP, cl){
  #xscale:

    cl$accumpos <- cl$BP
    cl <- cl[order(cl$CHR), ]
    for (i in 2:10) {
        cl[cl$CHR == i, ]$accumpos <- cl[cl$CHR == (i - 1), ]$accumpos + cl[cl$CHR == i, ]$accumpos + GAP
    }

  cl$ticks <- cl$accumpos[1]/2
  cl$chrlines <- cl$accumpos[1]+GAP/2
  for(i in 2:10){
    cl$ticks[i] <- cl$accumpos[i-1] + (cl$accumpos[i]-cl$accumpos[i-1])/2;
    cl$chrlines[i] <- cl$accumpos[i]+ GAP/2;
  }
  return(cl)
}

