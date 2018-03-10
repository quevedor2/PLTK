demo.sigClusterBreakpoints <- function(){
  require(GenomicRanges)
  if(!exists("path.to.git")) stop("Please define your path.to.git variable (i.e. path.to.git <- '~/git')")
  source(file.path(path.to.git, "/PLTK/preprocessing_tools/utils.R"))
  
  #Generate random data
  ends <- c(runif(n = 10, min = 0, max = 1000),
            rnorm(10, mean=1200, sd=100),
            runif(n = 10, min = 0, max = 1000))
  ends <- sort(as.integer(ends))
  starts <- c(1, ends[-length(ends)]+1)
  intervals <- data.frame("chr"=rep("chr1", length(starts)),
                          "start"=starts,
                          "ends"=ends)
  
  #Dataframe to Granges object
  intervals.gr <- dataframeToGranges(intervals)
  
  # Analysis
  binsize <- 50
  sigClusterBreakpoints(intervals.gr, binsize)
}


sigClusterBreakpoints <- function(gr, binsize){
  require(GenomicRanges)
  require(rlist)
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
