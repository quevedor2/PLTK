#----------------------------------------------------------------------------------------
#' GR Viz: Main for plotting GRanges objects
#' @description A general wrapper for plotting GRanges objects assuming the elementMetadata() contains the numeric values of plotting.  Specifies the analysis/visualization type and it will generate a plot with the aforementioned data
#'
#' @param gr [GRanges]: Input GRanges object
#' @param plot.settings [List]: Output of initializeGrPlot():  Meant to keep all GRanges plots aligned
#' @param data.type [Character]: Data type:  IGV Copy-number Tracks ("cn"), "expr", "mutation"
#' @param target.chr [Character]: Target chromosome to plot.  Must be the same one specified in initializeGrPlot()
#' @param add.axis [Boolean]: Whether to add an x-axis for the plot
#' @param side [Integer]: What side to put the axis (1=below, 3=above) [Default:1]
#' @param axis.marks [Integer]: How many axis tick marks to add [Default:5]
#' @param yrange [Integer Vector]: Your ylims for plot() function.  Everything will scale to fit within this range [Default: c(0,1)]
#' @param y.spacer [Numeric]: The amount of space to add between CN tracks; for example, a spacer of 0.1 for one sample in a yrange of c(0,1) will add a 0.1 gap and a 0.05 gap for 2 samples [Default: 0.1]
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plotGrMetadata <- function(gr,  plot.settings, data.type='cn', target.chr=NULL, 
                           add.axis=TRUE, side=1, axis.marks=5, yrange=c(0,1), y.spacer=0.1,
                           anno.track=0, add.annotations=FALSE, ...){
  if(!is.null(target.chr)) if(!grepl("^chr", target.chr)) target.chr <- paste0('chr', target.chr)
  #if(!is.null(target.chr)) gr <- gr[seqnames(gr) == target.chr]
  
  chr.ends <- plot.settings[['cumchr']]
  xrange <- plot.settings[['xlims']]
  
  all.chrs <- as.character(seqnames(gr)@values)
  
  if(is.null(target.chr)) {
    chr.ords <- seq_along(all.chrs)
  } else {
    chr.ords <- grep(paste0("^", target.chr, "$"), all.chrs)
  }
  
  if(add.annotations && data.type=='expr'){
    anno.yrange  <- c(min(yrange), max(yrange) + (diff(yrange) * anno.track))
    plot(0, type='n', xlim=xrange, ylim=anno.yrange, axes=FALSE, xlab='', ylab='', ...)
  } else {
    plot(0, type='n', xlim=xrange, ylim=yrange, axes=FALSE, xlab='', ylab='', ...)
  }
  
  for(each.chr.idx in chr.ords){
    chr.end <- chr.ends[each.chr.idx,'ends']
    gr.chr <- gr[seqnames(gr) == all.chrs[each.chr.idx]]
    if(add.axis & !is.null(target.chr)) axis(side = side, 
                                             at=seq(xrange[1], xrange[2], by=diff(xrange)/axis.marks),
                                             labels = round((seq(1, chr.end, 
                                                                 by=diff(c(1, chr.end))/axis.marks)/1000000), 0))
    
    # Calculate the adjusted start/end index
    adj.x <- chr.ends[each.chr.idx, 'cumStarts']
    print(adj.x)
    switch(data.type,
           cn=addCnSegs(gr.chr=gr.chr, adj.x=adj.x, ...),
           expr=addExprScores(gr.chr=gr.chr, yrange=yrange, adj.x=adj.x, ...),
           mutation=addMut())
  }
}

#----------------------------------------------------------------------------------------
#' GR Viz: Plotting function to add copy-number segments
#' @description Works in tandem with plotGrMetadata() to add coloured CN rectangles to the existing plot
#'
#' @param gr.chr [GRanges]: Copy-number GRanges object with CN data in elementMetadata
#' @param col.ids [Vector]: Which columns/samples to plot and in what order
#' @param adj.x [Integer]: Pulled from plot.settings to adjust the x-coordinates for the right chromosome
#' @param cn.colors [Char. Vector]: Range of colours to feed into ColorRampPalette [Default: c("blue", "white", "red")]
#' @param cn.range [Num. Vector]: Numeric vector of colour ranges [Default: seq(-3, 3, by=0.5)]
#' @param chr.lines [Boolean]: Whether to add dashed grey lines indicating edges of chromosomes
#' @param yrange [Int. Vector]: The ylim plotting space
#' @param y.spacer [Numeric]: How much space to add between tracks
#'
#' @return
#' @export
#'
#' @examples
addCnSegs <- function(gr.chr, col.ids=c(1), adj.x=0, 
                      cn.colors = c("blue", "white", "red"),
                      cn.range=NULL, chr.lines=FALSE,
                      yrange=c(0, 1), y.spacer=0.1){
  mapCol <- function(x){
    cols <- as.character(cn.colors[match(x, cn.colors$range),]$colors)
    cols[is.na(cols)] <- 'grey'
    cols
  }
  # Range of copy number values for colours
  if(is.null(cn.range)) cn.range <- seq(-3, 3, by=0.5)
  cn.colors <- data.frame("range"=as.numeric(cn.range),
                          "colors"=colorRampPalette(cn.colors)(length(cn.range)))
  
  # Gets the sgement start/end indices
  s.idx <- start(gr.chr) + adj.x
  e.idx <- end(gr.chr) + adj.x
  
  # Calculates the scaling factor to stack X number of samples in the ylim range
  y.space.sample <- diff(yrange) / length(col.ids)
  y.spacer <- y.space.sample * y.spacer
  
  for(col.idx in seq_along(col.ids)){
    each.col <- col.ids[col.idx]
    
    # Calculate the y-coordinates for the copy-number rectangles
    sample.pos <- (y.space.sample * (col.idx - 1)) # Stack samples on top of each other
    bot.idx <- min(yrange) + sample.pos + y.spacer
    top.idx <- y.space.sample + sample.pos - y.spacer
    
    print(paste0("Plotting CN for: ", each.col))
    score <- unlist(elementMetadata(gr.chr[,each.col])@listData)
    rect(xleft = s.idx, ybottom = rep(bot.idx, length(s.idx)), 
         xright = e.idx, ytop = rep(top.idx, length(s.idx)),
         col = mapCol(score), border = NA)
    
    if(chr.lines) abline(v = min(s.idx), lty=2, col="grey")
    if(chr.lines) abline(v = max(e.idx), lty=2, col="grey")
  }
}


#----------------------------------------------------------------------------------------
#' GR Viz: Initializes the parameters for GRanges visualization and Cytoband tracks
#' @description Main function is to initialize the adjusted x coordinates for all chromosomes but also optionally plots a track chromosome and cytoband/centromere positions
#'
#' @param cytoband.gr [GRanges]: A Granges object containing all UCSC cytobands, for example, PLTK::hg19.cytoband
#' @param target.chr [Character]: An optional target chromosome to plot.  If no chromosome is given, the entire genome is plotted [Default: NULL]
#' @param plot.chrom [Boolean]: Whether to generate a plot or not
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
initializeGrPlot <- function(cytoband.gr, target.chr=NULL, 
                             plot.chrom=TRUE, ...){
  cytoband.gr <- sort(cytoband.gr)
  genome.l <- sum(as.numeric(width(cytoband.gr)))
  
  # Calculate the cumulative starts for all chromosomes in the cytoband matrix
  all.chrs <- as.character(seqnames(cytoband.gr)@values)
  chr.ends <- sapply(all.chrs, function(x) max(end(cytoband.gr[seqnames(cytoband.gr)==x])))
  cum.chr.ends <- c(0, cumsum(as.numeric(chr.ends)))
  cum.chr.ends <- cum.chr.ends[-length(cum.chr.ends)]
  chr.ends <- matrix(c(chr.ends, cum.chr.ends), ncol=2, dimnames=list(all.chrs, c("ends", "cumStarts")))
  
  if(is.null(target.chr)) {
    xrange <- c(0, genome.l)
    chr.ords <- seq_along(all.chrs)
  } else {
    if(!grepl("^chr", target.chr)) target.chr <- paste0('chr', target.chr)
    
    row.idx <- grep(paste0("^(chr)?", target.chr, "$"), rownames(chr.ends))
    xrange <- c(chr.ends[row.idx, 'cumStarts'], 
                (chr.ends[row.idx, 'ends'] + chr.ends[row.idx, 'cumStarts']))
    chr.ords <- grep(paste0("^", target.chr, "$"), all.chrs)
  }
  
  if(plot.chrom){
    plot(0, type='n', xlim=xrange, ylim=c(0,1), axes=FALSE, xlab='', ylab='')
    for(each.chr.idx in chr.ords){
      cytoband.chr <- cytoband.gr[seqnames(cytoband.gr) == all.chrs[each.chr.idx]]
      
      # Calculate the adjusted start/end index
      adj.x <- chr.ends[each.chr.idx, 'cumStarts']
      s.idx <- start(cytoband.chr) + adj.x
      e.idx <- end(cytoband.chr) + adj.x
      
      addCytobands(s.idx, e.idx, cytoband.chr=cytoband.chr, chr.id=all.chrs[each.chr.idx], ...)
    }
  }
  return(list(xlims=xrange,
              cumchr=chr.ends))
}

#----------------------------------------------------------------------------------------
#' GR Viz: Generates the colours for cytobands
#'
#' @param gie [Char. Vector]: A vector containing all the ordered UCSC cytoband Giemsa colour
#' @param gie.col [Character]: Colour for giemsa positive segments
#' @param cen.col [Character]: Colour for centromere segments
#' @param plot.cband [Boolean]: Toggle for plotting cytobands or not; will always plot centromeres
#' @param alpha.factor [Integer]: How many times to dilute the intensity of the colours, i.e. 2 will dilute the colours 2x [Default: 2]
#' @param allow.alpha [Boolean]: scales::alpha doesn't work on Samwise, so this allows the trigger to remove the opacity factor in the hex code
#'
#' @return
#' @export
#'
#' @examples
getGieCol <- function(gie, gie.col='black', cen.col='red', 
                      plot.cband=TRUE, alpha.factor=2,
                      allow.alpha=TRUE){
  require(scales)
  cen.idx <- grep("acen|gvar", gie)
  gie <-gsub("acen|gvar|stalk", 50, gie)
  
  if(plot.cband){
    gie <- gsub("gneg", 0, gie)
    gie <- gsub("gpos", "", gie)
  } else {
    gie <- gsub("gneg|gpos.*$", 10, gie)
  }
  
  
  gie <- (as.numeric(gie) / 100) / alpha.factor
  
  cols <- rep(gie.col, length(gie)) 
  cols[cen.idx] <- cen.col
  cols <- scales::alpha(cols, gie)
  if(!allow.alpha) cols <- gsub("..$", "", cols)
  if(!allow.alpha) cols[c(1:length(cols))[-grep("#FF0000", cols)]] <- "white"
  cols
}

#----------------------------------------------------------------------------------------
#' GR Viz: Plotting function to add chromosome cytobands
#' @description Adds the rectangles for each chromosomal cytobadn and centromeres
#'
#' @param s.idx [Int. Vector]: All the start x-coords for cytoband/centromere segments
#' @param e.idx [Int. Vector]: All the end x-coords for cytoband/centromere segments
#' @param bot.idx [Numeric]: Y coordinate for the bottom of the plotted cytoband rectangle [Default:0.1]
#' @param top.idx [Numeric]: Y coordinate for the top of the plotted cytoband rectangle [Default:0.9]
#' @param chr.id [Character]: Chromosome ID
#' @param label.side [Character]: Whether to add the chromosome label above ("top") or below ("bottom") the cytobands
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
addCytobands <- function(s.idx, e.idx, bot.idx=0.1, top.idx=0.9, 
                         chr.id="NA", label.side='top', cytoband.chr=NA, ...){
  min.sidx <- min(s.idx)
  max.eidx <- max(e.idx)
  
  switch(label.side,
         top=lbl.idx <- (top.idx + (mean(c(bot.idx, top.idx))/5)),
         bottom=lbl.idx <- (bot.idx - (mean(c(bot.idx, top.idx))/5)))
  
  text(x = mean(c(max.eidx, min.sidx)), y = lbl.idx, 
       labels = gsub("^chr", "", chr.id), adj=1, ...)
  rect(xleft = min.sidx, ybottom = bot.idx, xright = max.eidx, ytop = top.idx, 
       col = "white", ...)
  rect(xleft = s.idx, ybottom = bot.idx, xright = e.idx, ytop = top.idx, 
       col = getGieCol(cytoband.chr$gieStain, ...), border = NA)
}

#----------------------------------------------------------------------------------------
#' GR Viz: Plotting function to add Expression ZScores
#' @description Still in the developmental build: Adds zscores of matrix for a given single sample where elementMetadata contains both symbol for Gene HUGO symbol and zscore columns
#'
#' @param gr.chr [GRanges]: Input GRanges.  Contains elementMEtadata() columns of "zscore" and "symbol"
#' @param yrange [Int. Vector]: The y-axis range to plot on; i.e. c(0, 5)
#' @param adj.x [Integer]: Pulled from plot.settings to adjust the x-coordinates for the right chromosome
#' @param anno.track [Numeric]: The fraction of the plot to allocate for an extra annotation track
#' @param add.y.axis [Boolean]: Whether to add an y-axis with z-score values
#' @param add.annotations [Boolean]: Whether to add annotations for each zscore range
#'
#' @return
#' @export
#'
#' @examples
addExprScores <- function(gr.chr, yrange, adj.x, anno.track=0.25,
                          add.y.axis=FALSE, add.annotations=TRUE, ...){
  if(add.annotations) {
    anno.yrange  <- c(min(yrange), max(yrange) + (diff(yrange) * anno.track))
  } else {
    anno.yrange <- yrange
  }
  # Gets the sgement start/end indices
  s.idx <- start(gr.chr) + adj.x
  e.idx <- end(gr.chr) + adj.x
  
  mid.idx <- apply(matrix(c(s.idx, e.idx), ncol=2), 1, mean)
  print(gr.chr$zscore)
  points(x = mid.idx, y=gr.chr$zscore, 
         pch=19, col="black")
  
  abline(h = 0, col="grey", lty=1)
  if(add.annotations){
    abline(h = max(yrange), lty=3, col="grey")
    text(x = mid.idx, y=(max(yrange) + anno.track / 10), 
         labels = gr.chr$symbol, adj=0, srt=90, ...)
  }
  if(add.y.axis) axis(side =2, las=2)
}
