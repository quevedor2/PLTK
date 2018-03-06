
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


# Plot scatterplot/boxplot/violin plot with two groups with or without connecting lines
# groups = a labeled matrix, dataframe or list containing multiple groups (cols)
# targ.pnt = Individual targets based on rownames or labelled vectors in list
# top.anno = Annotate the top X differences between groups [default=10]
# y.lim = provide only if specific min/max required c(min,max)
# withVioin = boolean if violin contour to be drawn
# withBox = boolean if box contour to be drawn
# vio.col = vector containing custom colors for violin fill
# connect = boolean if dots should be connected
plotScatterLine <- function (groups,
                             targ.pnt=NA,
                             top.anno=10,
                             y.lim = NA, 
                             withViolin = F, 
                             withBox = F, 
                             vio.col = "lightblue", 
                             connect = T,
                             ...){
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
  
  plot(0,type='n',axes=FALSE,ann=FALSE, xlim = c(0, (length(groups)+2)), ylim = y.lim)
  title (font = 2, cex.main = 1.5, cex.lab = 1.5, ...)
  axis (side = 1, lwd = 2, labels = F, lwd.tick=0)
  axis (side = 1, at = seq(1, length(groups)), labels = names(groups), cex.axis = 1.5, font = 2, lwd = 2)
  axis (side = 2, lwd = 2, las = 1, font = 2, cex.axis = 1.2)
  
  groups <- lapply(groups, na.omit)
  
  addMedianSeg <- function(x, y){
    segments(x0 = x-0.2, x1 =x+0.2 , y0 = median(y), y1 = median(y) , col = "black", lwd = 2)
  }
  
  if (withViolin){
    library(vioplot)
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
    # Formats links between groups for specific targets or not
    if(!any(is.na(targ.pnt))){
      groups <- lapply(groups, function(x) x[targ.pnt])
    }
    seg.connect <- do.call(cbind, groups)
    seg.connect <- seg.connect[which(complete.cases(seg.connect)), , drop=FALSE]
    
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
    max.diff <- linkSegs(seg.connect, lty=1, col=alpha("black", 0.50))

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
      top.seg.connect <- top.seg.connect[order(top.seg.connect[,ncol(top.seg.connect)]),]
      if(nrow(top.seg.connect) != 0){
        annoSegs(top.seg.connect, groups)
      } else {
        warning("No annotations provided.  Did you remember to label your rows?")
      }
      
      
    }
  }
  return (y.lim)
}
