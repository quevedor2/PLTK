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

plotLikelihoodRatio <- function(groups, xmat, 
                            f1col='black', f2col='red', 
                            gen.plot=TRUE, upp.y=2, low.y=-0.2,
                            parmar=c(5.1, 4.1, 4.1, 4.1),
                            addSegments=T, seg.thresh=2, col.idx=NA,
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
    lines(dseq1$x, dseq1$y * scale.val, col=alpha(f1col, 0.6), yaxt='n', ...)
    lines(dseq2$x, dseq2$y * scale.val, col=alpha(f2col, 0.6), yaxt='n', ...)
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
           col = alpha(f1col,0.5), border = NA)
      rect(xleft = x, ybottom = 0, xright = x + x.width, ytop = (p2 * scale.val), 
           col = alpha(f2col,0.5), border = NA)
    }
  }
  
  return(list('ks'=ks.pval,
              'bf'=do.call("cbind", bf.vals)))
}
demo.plotLikelihoodRatio()
