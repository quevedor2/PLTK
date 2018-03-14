#' Vis.demo: demo.plotLikelihoodRatio
#' @description An example piece of code to run the plotLikelihoodRatio() function
#' @return
#' @export
#'
#' @examples demo.plotLikelihoodRatio()
demo.plotLikelihoodRatio <- function(){
  z <- data.frame("Group A"=rnorm(n = 1000, mean = 10, sd = 5),
                  "Group B"=rpois(n = 1000, lambda = 2),
                  "Group C"=rpois(n = 1000, lambda = 2))
  z[1,1] <- NA
  z[5, 2] <- NA
  rownames(z) <- sample(sapply(letters, function(x) paste0(x, c(1:100))), size = 1000, replace=F)
  
  samples <- matrix(c(0.4, 15), nrow=1, dimnames = list('paclitaxel', c('MCF-7', "HeLa")))
  print(samples)
  
  plotLikelihoodRatio(z, samples, bw=1, xlab="Delta AAC", col.idx=c(1,2))
}

#' Visualization: plotLikelihoodRatio
#' @description Creates a KDE curve for two datasets and calculates the likelihood ratio between the curves by calculating the natural log-ratio of point-probability estimates between the two datasets. Originally used to quantify the differences in drug-sensitivity for a given cell line (REH) between two datasets (CCLE and GDSC) for paclitaxel.
#'
#' @param groups A list of numerical vectors (groups) or a matrix/dataframe with each column being its own group
#' @param xmat A 1-row matrix or vector with identification of x-values to calculate point likelihood ratios
#' @param f1col Group 1 colour
#' @param f2col Group 2 colour
#' @param gen.plot Generates plots
#' @param upp.y Scaling for plotting the y-axis; keeps all y-values proportional
#' @param low.y Scaling for plotting the y-axis; keeps all y-values proportional
#' @param parmar par()$mar values if you want to change default
#' @param addSegments Boolean on whether to add rectangular segments showing heights of likelihoods
#' @param col.idx The index of the two columns to use if given a dataframe/matrix with more than 2 columns
#' @param ... Passed into density() and plot() function
#'
#' @return
#' @export
#' @importFrom scales alpha
#' 
#' @examples
#' demo.plotLikelihoodRatio()
plotLikelihoodRatio <- function(groups, xmat, 
                                f1col='black', f2col='red', 
                                gen.plot=TRUE, upp.y=2, low.y=-0.2,
                                parmar=c(5.1, 4.1, 4.1, 4.1),
                                addSegments=T, col.idx=NA,
                                ...){
  require(scales)
  par(mar=parmar)
  # Converts matrices/dataframes to a list
  if(is.matrix(groups) || is.data.frame(groups)){
    groups <- as.matrix(groups)
    if(!is.na(col.idx)) groups <- groups[,col.idx]
    tmp <- list()
    for(i in seq(1:ncol(groups))){
      tmp[[colnames(groups)[i]]] <- groups[,i]
    }
    groups <- tmp
  }
  if(length(groups) != 2) stop("Only two groups can be passed into this function; try specifying col.idx if you are using a dataframe or matrix")
  
  # Setting up the KDE and approximated functions
  seq1 <- na.omit(groups[[1]])
  seq2 <- na.omit(groups[[2]])
  dseq1 <- density(na.omit(seq1), ...)
  dseq2 <- density(na.omit(seq2), ...)
  fn1.pdf <- approxfun(dseq1)
  fn2.pdf <- approxfun(dseq2)
  
  ks.pval <- ks.test(seq1, seq2)$p.val
  if(gen.plot){
    max.p <- max(c(dseq1$y, dseq2$y), na.rm=TRUE)
    scale.val <- upp.y / max.p
    
    xrange <- summary(c(dseq1$x, dseq2$x))
    yrange <- summary(c(dseq1$y, dseq2$y))
    xposQuant <- function(xsummary, quant){
      xsummary['Min.'] + ((xsummary['Max.'] - xsummary['Min.']) * quant)
    } 
    
    
    
    plot(x=c(xrange['Min.'], xrange['Max.']), y=c(yrange['Min.'], yrange['Max.']), 
         type='n', ylim=c(low.y, upp.y+0.2), axes=FALSE,
         ylab='', yaxt='n', ...)
    mtext(paste(names(groups)[1], "density"), side = 2, line=3)
    mtext(paste(names(groups)[2], "density"), side = 4, line=3)
    axis(side = 1)
    axis(side = 2, at = seq(0, upp.y, by=upp.y/5), labels=round(seq(0, max(dseq1$y), by=max(dseq1$y)/5), 2), las=2)
    axis(side = 4, at = seq(0, upp.y, by=upp.y/5), labels=round(seq(0, max(dseq2$y), by=max(dseq2$y)/5), 2), las=2)
    lines(x = dseq1$x, y = dseq1$y * scale.val, col=scales::alpha(f1col, 0.6), yaxt='n', ...)
    lines(dseq2$x, dseq2$y * scale.val, col=scales::alpha(f2col, 0.6), yaxt='n', ...)
    text(x=rep(xrange['Min.'], 2), y=c(upp.y+0.2, low.y), adj=0,
         labels=names(groups), col=c(f1col, f2col))
    
    #KS-Statistic
    text(x=xposQuant(xrange, 0.87), y=(upp.y+0.1), labels="KS-Statistic:", adj=1, cex=0.9)
    text(x=xposQuant(xrange, 0.9), y=(upp.y+0.1), labels=paste0("p= ", round(ks.pval,4)), adj=0, cex=0.9)
  }
  
  # Bayes factor for each set of points
  bf.vals <- list()
  row.cnt <- 1
  if(is.vector(xmat)) {xmat <- t(as.matrix(xmat)); colnames(xmat) <- as.character(xmat[1,])}
  for(each.col in seq(ncol(xmat))){
    # Obtain post estimates for probability of x in both curves
    x <- xmat[,each.col]
    p2 <- fn2.pdf(x)
    p1 <- fn1.pdf(x)
    if(!is.na(x) && is.na(p2)) p2 <- 0
    if(!is.na(x) && is.na(p1)) p1 <- 0
    bf <- log(p2/p1)
    bf.id <- colnames(xmat)[each.col]
    bf.vals[[bf.id]] <- bf
    
    # Plot the point prob estimates and Likelihood Ratios
    if(!gen.plot || is.na(x)) next
    text(x=xposQuant(xrange, 0.87), y=(upp.y + 0.1 - row.cnt*0.1), adj=1,
         labels=paste0(bf.id, ":"), cex=0.9)
    text(x=xposQuant(xrange, 0.9), y=(upp.y + 0.1 - row.cnt*0.1), adj=0, cex=0.9,
         labels=paste0("LR= ", if(is.na(x)) 'NA' else round(bf, 2)))
    
    if(is.na(x)) next
    points(rep(x, 2), c((p1 * scale.val), (p2 * scale.val)), col=c(f1col, f2col), pch=19)
    text(x=rep(x, 2), y=c(upp.y+0.2, low.y), col=c(f1col, f2col),
         labels=c(round(p1, 2), round(p2, 2)), cex=0.9)
    mtext(text = bf.id, side = 3, at=x, line=0)
    row.cnt <- row.cnt + 1
    
    if(addSegments){
      x.width <- (xrange['Max.'] - xrange['Min.']) * 0.01
      rect(xleft = x - x.width, ybottom = 0, xright = x, ytop = (p1 * scale.val), 
           col = scales::alpha(f1col,0.5), border = NA)
      rect(xleft = x, ybottom = 0, xright = x + x.width, ytop = (p2 * scale.val), 
           col = scales::alpha(f2col,0.5), border = NA)
    }
  }
  
  return(list('ks'=ks.pval,
              'bf'=do.call("cbind", bf.vals)))
}


#' Vis.demo: demo.plotScatterLine
#' @description An example piece of code to run the plotScatterLine() function
#' @return
#' @export
#'
#' @examples demo.plotScatterLine()
demo.plotScatterLine <- function(){
  z <- data.frame("A"=sample(x = seq(0,1,by=0.001), size = 100,replace = T),
                  "B"=sample(x = seq(0,1,by=0.001), size = 100,replace = T),
                  "C"=sample(x = seq(0,1,by=0.001), size = 100,replace = T))
  z[5, 3] <- NA
  rownames(z) <- sample(sapply(letters, function(x) paste0(x, c(1:5))), size = 100, replace=F)
  plotScatterLine(z, vio.col=c("red", "blue", "green"),
                  withViolin=T, connect=T, 
                  main="Main title", xlab="X label", ylab="Y label")
}

#' Visualization: plotScatterLine
#' @description Creates a stripchart, boxplot, or violin plot for each group and then offers the option to link changes between individual elements between groups. If rownames or a labelled vector in a list is provided, it will annotate the top X differences between groups based on taking the max difference for a given element across all groups.
#'
#' @param groups a labeled matrix, dataframe or list containing multiple groups (cols)
#' @param targ.pnt Individual targets based on rownames or labelled vectors in list
#' @param top.anno Annotate the top X differences between groups [default=10]
#' @param y.lim provide only if specific min/max required c(min,max)
#' @param withViolin boolean if violin contour to be drawn
#' @param withBox boolean if box contour to be drawn
#' @param vio.col vector containing custom colors for violin fill
#' @param connect boolean if dots should be connected
#' @param show.all Boolean to show all connecting segments if targ.pnt is provided
#' @param ... 
#'
#' @return
#' @export
#' @import vioplot
#' @importFrom scales alpha
#' 
#' @examples demo.plotScatterLine()
plotScatterLine <- function (groups,
                             targ.pnt=NA,
                             top.anno=10,
                             y.lim = NA, 
                             withViolin = F, 
                             withBox = F, 
                             vio.col = "lightblue", 
                             connect = T,
                             show.all=T,
                             ...){
  require(scales)
  if(is.matrix(groups) || is.data.frame(groups)){
    groups <- as.matrix(groups)
    tmp <- list()
    for(i in seq(1:ncol(groups))){
      tmp[[colnames(groups)[i]]] <- groups[,i]
    }
    groups <- tmp
  }
  if (length(y.lim)<2) { y.lim = c(min(unlist(groups), na.rm=T)-0.1,
                                   max (unlist(groups), na.rm = T)+ (max (unlist(groups), na.rm = T)*.20)) }
  
  par(las=1,bty="l")  ## my preferred setting
  plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0, (length(groups)+2)), ylim = y.lim, ...)
  title (font = 2, cex.main = 1.5, cex.lab = 1.5, ...)
  axis (side = 1, lwd = 2, labels = F, lwd.tick=0)
  axis (side = 1, at = seq(1, length(groups)), labels = names(groups), cex.axis = 1.5, font = 2, lwd = 2)
  axis (side = 2, lwd = 2, las = 1, font = 2, cex.axis = 1.2)
  
  groups <- lapply(groups, na.omit)
  
  addMedianSeg <- function(x, y){
    segments(x0 = x-0.2, x1 =x+0.2 , y0 = median(y), y1 = median(y) , col = "black", lwd = 2)
  }
  
  if (withViolin){
    require(vioplot)
    #source("~/PughLab/INSPIRE/CCRC/Poster_Analysis/vioplot.R") # change this for custom colors
    lapply(seq_along(groups), function(cnt.x){
      require(scales)
      vioplot(groups[[cnt.x]], at=cnt.x, col = adjustcolor(vio.col[cnt.x], alpha = 0.3), 
              horizontal = F, add = T,frame=F, drawRect = F)
      addMedianSeg(cnt.x, groups[[cnt.x]])
    })
  }
  
  if (withBox){
    boxplot (groups, add = T, outline = F, boxwex = 0.5, axes=FALSE, frame=F, col = adjustcolor(vio.col, alpha = 0.5))
  }
  stripchart (groups, at =seq(1, length(groups)), pch = 16, vertical = T, add = T, cex = 1.5)
  lapply(seq_along(groups), function(cnt.x) addMedianSeg(cnt.x, groups[[cnt.x]]))
  
  if (connect){
    # Connects the segments
    segConnect <- function(groups){
      seg.connect <- do.call(cbind, groups)
      seg.connect <- seg.connect[which(complete.cases(seg.connect)), , drop=FALSE]
      seg.connect
    }
    g.seg.connect <- segConnect(groups)
    
    # Formats links between groups for specific targets or not
    if(!any(is.na(targ.pnt))){
      groups <- lapply(groups, function(x) x[targ.pnt])
      seg.connect <- segConnect(groups)
    } else {
      seg.connect <- g.seg.connect
    }
    
    # Link all segments together using a line, returns the top sets with biggest deltas
    linkSegs <- function(s.conn, ...){
      apply(s.conn, 1, function(e.row){
        for(e.pnt in seq(2, length(e.row))){
          segments(x0=e.pnt-1, y0=as.numeric(e.row[e.pnt-1]),
                   x1=e.pnt, y1=as.numeric(e.row[e.pnt]), ...)
        }
        max(abs(diff(as.numeric(e.row))))
      })
    }
    if(show.all) max.diff <- linkSegs(g.seg.connect, lty=1, col=scales::alpha("black", 0.20))
    max.diff <- linkSegs(seg.connect, lty=1, col=scales::alpha("black", 0.50))
    
    # Annotate the Top X changes
    if(!any(is.na(top.anno))){
      if(length(max.diff) < top.anno) top.anno <- length(max.diff)
      top.seg.connect <- seg.connect[names(sort(max.diff, decreasing = T)[1:top.anno]), ,drop=FALSE]
      
      annoSegs <- function(s.conn, groups){
        spacer <- round(1 / nrow(s.conn), 2)
        linkSegs(s.conn, lwd=2, col="red")
        
        
        lapply(seq(nrow(s.conn)), function(e.cnt){
          anno.ypos <- (0 + (spacer * e.cnt))
          anno.xpos <- length(groups) + 0.75
          point.ypos <- s.conn[e.cnt, length(groups)]
          
          text(x = anno.xpos,  y=anno.ypos,
               labels = rownames(s.conn)[e.cnt], adj = 0)
          segments(x0 =  length(groups), y0 = point.ypos,
                   x1 = anno.xpos, y1 = anno.ypos, col = "dark grey", lty=2)
        })
      }
      top.seg.connect <- top.seg.connect[order(top.seg.connect[,ncol(top.seg.connect)]),,drop=FALSE]
      if(nrow(top.seg.connect) != 0){
        annoSegs(top.seg.connect, groups)
      } else {
        warning("No annotations provided.  Did you remember to label your rows?")
      }
      
      
    }
  }
  return (y.lim)
}
