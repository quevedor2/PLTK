#' cnSignature: sigClusterBreakpoints
#' @description Takes a list of copy-number segments and tries to identify regions where there are consecutive segments less than a pre-designed segment size. Within these regions, it finds the longest string of consecutive segments that are less than the pre-designed segment size, annotates it, and reports them back in a list for downstream analysis.
#'
#' @param gr A granges object
#' @param binsize A set binsize (bp)
#'
#' @return \code{segs.list}: A list of genomicRanges objects for genomic regions with continuous segments smaller than binsize
#' @examples
#' sigClusterBreakpoints(demo, 50)
sigClusterBreakpoints <- function(gr, binsize){
  require(GenomicRanges)
  # Subsets the granges object given the specifications of st/end ranges
  mkSubGranges <- function(){
    grx <- gr[(st+1):(en+2),]
    bin.id <- rep(paste0("bin", each.idx), length(seqnames(grx)))
    elementMetadata(grx) <- data.frame("Bin"=bin.id,
                                       "sub"=rep(paste0("sub", bin.cnt), length(bin.id)))
    grx
  }
  
  # RLE of segment sizes
  segsizes <- diff(end(gr))
  posidx <- getRleIdx(segsizes < binsize)
  segs.list <- list()
  
  # Cycle through all consecutive segments less than a set binsize [default=50]
  for(each.idx in which(as.logical(posidx$values))){
    tbin <- 100   # Ridiculous large number to instantiate
    bin.cnt <- 1  # Labelling of sub-fragments with larger bin
    b.id <- paste0("bin", each.idx)
    segs.per.bin <- c()
    
    st <- posidx$start.idx[each.idx]
    en <- posidx$end.idx[each.idx]
    END.NOT.REACHED <- TRUE
    while(END.NOT.REACHED){
      # Finds the largest number of segments that are in total, less than binsize
      while(tbin > binsize){
        tbin <- sum(segsizes[st:en])
        en <- en - 1
      }
      if((en + 1) == posidx$end.idx[each.idx]) END.NOT.REACHED <- FALSE
      segs.per.bin <- c(segs.per.bin, (diff(c(st,en+1)) + 1)) # Printing the count
      
      # Using these indices, finds the original segments
      if(is.null(segs.list[[b.id]])){
        segs.list[[b.id]] <- mkSubGranges()
      } else {
        segs.list[[b.id]] <- append(segs.list[[b.id]],
                                    mkSubGranges())
      }
      
      # Cycles through the next list of segments if not at the end of the bin
      st <- en + 2
      en <- posidx$end.idx[each.idx]
      bin.cnt <- bin.cnt + 1
      tbin <- 100
    }
    
    print(segs.per.bin)
  }
  
  return(segs.list)
}