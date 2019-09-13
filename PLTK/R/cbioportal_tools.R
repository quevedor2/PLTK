#' Plots Mutual exclusivity heatmaps
#' @description Plots the Mutual exclusivity table directly from
#' cBioportal, or the matrix that was assembled using the sister
#' function
#'
#' @param mat The matrix from genMEMatrix, or the table from table.tsv
#' @param min Min log2 odds ratio to plot
#' @param max Max log2 odds ratio
#' @param step Increment in steps from min to max
#' @param lo Colour fr the low LOR (mutual exclusive)
#' @param mid Colour for the mid LOR
#' @param hi Colour for the high LOR (co-occurence)
#'
#' @return Plots a heatmap or bubble-heatmap plot for the given ME matrix
#' @importFrom gplots heatmap.2
#' ggplot2 ggplot
#' @export
#'
#' @examples
#' plotMEMatrix(or.mat)
#'
#'m.or <- melt(or.mat, varnames = c("A", "B"), value.name = "Log2 Odds Ratio")
#' m.p <- melt(p.mat, varnames = c("A", "B"), value.name = "p-Value")
#' tbl2 <- merge(m.or, m.p, by=c("A", "B"))
#' plotMEMatrix(tbl2)
plotMEMatrix <- function(mat, min=-0.5, max=0.5, step=0.01,
                         lo="#2c7bb6", mid="white", hi="#1a9641"){
  if(class(mat) == 'matrix'){
    breaks <- seq(min, max, by=step)
    my_palette <- colorRampPalette(c(lo, mid, hi)) (n=length(breaks)-1)

    heatmap.2(mat, trace="none", na.color = "grey", scale="none",
              dendrogram='none', Rowv = FALSE, Colv=FALSE,
              col = my_palette, breaks=breaks,
              margins=c(10,10), cexRow=0.9, cexCol=0.9,
              keysize=1.3, key.title=NA , key.ylab=NA, density.info='none', key.xlab='Log2 Odds Ratio')
  } else if(class(mat) == 'data.frame'){
    ggplot(mat, aes(B, A)) +
      geom_point(aes(colour = `Log2 Odds Ratio`,
                     size = `p-Value`))  +
      scale_colour_gradient2(low = lo, mid=mid,
                             high = hi, midpoint=0,
                             limits=c(min,max)) +
      scale_size(range = c(30,5)) +
      theme_bw() +
      theme(axis.title=element_blank())


  }
}

#' Genreates Mutual exclusivitiy matrices
#' @description Takes the table.tsv downloadable file from cBioportal
#' and outputs a matrix of values for A and B columns
#'
#' @param tbl The table.tsv downloaded from cBioportal Mutual Exclusivity
#' tab and imported via read.table with headers
#' @param val Value that you want to encode in the matrix
#'
#' @return Returns a matrix of A by B for given Val
#' @export
#'
#' @examples
#' or.mat <- genMEMatrix(tbl, val='Log2 Odds Ratio')
genMEMatrix <- function(tbl, val='Log2 Odds Ratio'){
  require(reshape2)
  tbl$A <- gsub(":.*$", "", tbl$A)
  tbl$B <- gsub(":.*$", "", tbl$B)

  # Invert tbl to make symmetrical
  .inverseTbl <- function(tbl){
    tbl.inv <- tbl
    tbl.inv$A <- tbl$B
    tbl.inv$B <- tbl$A
    tbl.inv
  }
  itbl <- .inverseTbl(tbl)

  tbl <- rbind(tbl, itbl)

  tbl.mat <- dcast(tbl, A~B, value.var=val)
  tbl.mat[upper.tri(tbl.mat)] <- NA

  rownames(tbl.mat) <- tbl.mat[,1]
  as.matrix(tbl.mat[,-1])
}
