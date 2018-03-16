#' cnSignature: Wrapper to run all CN Signatures
#' @description A Wrapper function to run all copy-number GRanges objects through the CN signature functions
#'
#' @param gr [GRanges]: GRanges object with copy-number in the elementMetadata().  Tested on output from convertToGr()
#' @param binsize [Integer]: Maximum size of segments to be considered in a small-CN cluster [Default: 1000000]
#' @param bins [GRanges]: A pre-made GRanges object where the genome is binned into smaller segment [Default: PLTK::bins]
#'
#' @return [List]: List of all signatures List[[each sample]][[Signature1]]
#' @export
#'
#' @examples
runCnSignatures <- function(gr, binsize=1000000, bins=PLTK::bins,
                            gap=PLTK::hg19.centromeres, gap.type='centromeres'){
  sig.list <- lapply(seq_along(elementMetadata(gr)), function(sample.idx){
    sample.sig.list <- list()
    sample.gr <- collapseSample(gr, sample.idx)
    sample.id <- colnames(elementMetadata(gr))[sample.idx]
    sample.sig.list[['cluster.bp']] <- sigClusterBreakpoints(sample.gr, binsize)
    sample.sig.list[['binned.bp']] <- sigBinBreakpoints(sample.gr, bins)
    sample.sig.list[['gap.dist']] <- sigGapDist(sample.gr, gap = gap, gap.type = gap.type)
    sample.sig.list[['seg.size']] <- sigSegSize(sample.gr)
    sample.sig.list
  })
  sig.list
}

#----------------------------------------------------------------------------------------
#' cnSignature: sigClusterBreakpoints
#' @description Takes a list of copy-number segments and tries to identify regions where there are consecutive segments less than a pre-designed segment size. Within these regions, it finds the longest string of consecutive segments that are less than the pre-designed segment size, annotates it, and reports them back in a list for downstream analysis.
#'
#' @param gr A granges object
#' @param binsize A set binsize (bp)
#' @import GenomicRanges
#' @export
#' 
#' @return \code{segs.list}: A list of genomicRanges objects for genomic regions with continuous segments smaller than binsize
#' @examples
#' sigClusterBreakpoints(PLTK::genDemoData(), 50)
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
    tbin <- binsize*10000   # Ridiculous large number to instantiate
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
      tbin <- binsize*10000
    }
    
    print(segs.per.bin)
  }
  
  return(segs.list)
}


#----------------------------------------------------------------------------------------
#' cnSignature: sigBinBreakpoints
#' @description Takes a GRanges object of segments and counts the number of breakpoints found within pre-designed genomic bins

#' @param gr GRanges object of your segments
#' @param bins GRanges object of the genome binned into set segment sizes [i.e. PLTK::bins]
#'
#' @return A list containing two elements:
#'   \code{segs}: A list of granges object for the original \code{gr} placed into each \code{bins}
#'   \code{bins}: The original \code{bins} object with an additional column in the metadata indicating number of breakpoints
#' @export
#' @import GenomicRanges
#' 
#' @examples
#'  sigBinBreakpoints(PLTK::genDemoData(), PLTK::bins)
sigBinBreakpoints <- function(gr, bins){
  require(GenomicRanges)
  olaps <- GenomicRanges::findOverlaps(query = gr, subject = bins)
  idx <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
  
  split.bins <- splitAsList(gr[queryHits(olaps)], idx)
  elementMetadata(bins)$binnedBP <- sapply(split.bins, length)
  return(list("segs"=split.bins,
              "bins"=bins))
}


#----------------------------------------------------------------------------------------
#' cnSignature: Distance from Centromeres or Telomeres
#' @description Calculates the distances from the closest breakpoint end to the telomere or centromere.  Offers the option to normalize by the size of the chromosome.
#'
#' @param gr [GRanges]: GRanges object
#' @param gap.type [Character]: Either "centromeres" or "telomeres"
#' @param gap [GRanges]: GRanges gap data, either PLTK::centromeres or PLTK::telomeres
#' @param normalize [Boolean]: Normalize by chromosome lengths; uses PLTK::hg19.cytobands as default chromosome sizes
#' @param verbose [Boolean]: Prints out extra debug info
#'
#' @return List of distances from the telomere or centromere
#' @export
#'
#' @examples sigGapDist(demo, gap.type = "telomeres", gap = PLTK::hg19.telomeres, normalize=TRUE)
sigGapDist <- function(gr, gap=PLTK::hg19.centromeres, gap.type='centromeres', normalize=FALSE, verbose=FALSE){
  gap.dist <- lapply(as.character(seqnames(gr)@values), function(each.chr){
    gr.chr <- gr[seqnames(gr) == each.chr]
    gap.chr <- gap[seqnames(gap) == each.chr]
    
    centromere.chr <- PLTK::hg19.centromeres[seqnames(PLTK::hg19.centromeres) == each.chr]
    cytoband.chr <- PLTK::hg19.cytobands[seqnames(PLTK::hg19.cytobands) == each.chr]
    
    p.arm <- which((end(gr.chr) - start(centromere.chr)) < 0)
    q.arm <- which((start(gr.chr) - end(centromere.chr)) > 0)
    switch(gap.type,
           centromeres={
             if(verbose) print("Telomeres")
             p.dist.from.breakpoint <- start(gap.chr) - end(gr.chr[p.arm,])
             q.dist.from.breakpoint <- start(gr.chr[q.arm,])  - end(gap.chr)
           },
           telomeres={
             if(verbose) print("Telomeres")
             p.dist.from.breakpoint <- start(gr.chr[p.arm,]) - end(sort(gap.chr))[1]
             q.dist.from.breakpoint <- start(sort(gap.chr))[2] - end(gr.chr[q.arm,])
           },
           stop="Unknown gap type.  Must be either 'centromeres' or 'telomeres'")
    all.dist <- c(p.dist.from.breakpoint, q.dist.from.breakpoint)
    if(normalize) all.dist <- (all.dist / max(end(cytoband.chr)))
    
    return(all.dist)
  })
 names(gap.dist) <- seqnames(gr)@values
 gap.dist
}


#----------------------------------------------------------------------------------------
#' cnSignature: Segment sizes
#' @description Calculates the size of all the segments in the GRanges object.  OFfers the option to normalize by chromosome size.
#'
#' @param gr [GRanges]: GRanges object
#' @param normalize [Boolean]: Normalize by chromosome lengths; uses PLTK::hg19.cytobands as default chromosome sizes
#'
#' @return List of segment sizes
#' @export
#'
#' @examples sigSegSize(demo, normalize=TRUE)
sigSegSize <- function(gr, normalize=FALSE){
  segl <- lapply(as.character(seqnames(gr)@values), function(each.chr){
    cytoband.chr <- PLTK::hg19.cytobands[seqnames(PLTK::hg19.cytobands) == each.chr]
    
    gr.chr <- gr[seqnames(gr) == each.chr]
    segl.chr <- diff(end(gr.chr))
    
    if(normalize) segl.chr <- (segl.chr / max(end(cytoband.chr)))
    segl.chr
  })
  names(segl) <- seqnames(gr)@values
  return(segl)
}

#----------------------------------------------------------------------------------------
#' cnSignature: All copy numbers
#' @description Calculates the copy numbers in the sample copy-number data with the option to scale it based on the size of the segment
#'
#' @param gr [GRanges]: GRanges object
#' @param weight [Boolean]: Weight the copy-number for a segment based on the size of the segment
#' @param normalize [Boolean]: Normalize by chromosome lengths; uses PLTK::hg19.cytobands as default chromosome sizes
#'
#' @return Dataframe of copy-number values for each chromosome
#' @export
#' @import plyr
#'
#' @examples sigCopyNumber(gr, weight=TRUE, normalize=TRUE)
sigCopyNumber <- function(gr, weight=TRUE, normalize=FALSE){
  suppressPackageStartupMessages(require(plyr))
  if(length(elementMetadata(gr)) == 0) stop("Copy-number data per sample is required to be stored in elementMetadata() of the GRanges object")
  if(!all(elementMetadata(gr)[,1] %% 1 == 0, na.rm = TRUE)) stop("This signature requires all copy-number data to be absolute copy-states/integers.  Please ensure no log2ratios are used.")
  
  copystate <- lapply(as.character(seqnames(gr)@values), function(each.chr){
    cytoband.chr <- PLTK::hg19.cytobands[seqnames(PLTK::hg19.cytobands) == each.chr]
    gr.chr <- gr[seqnames(gr) == each.chr]
    
    # Get all copy-states int he elementMetadata
    copy.states <- na.omit(unique(elementMetadata(gr.chr)[,1]))
    # Cycles through and counts (+ weights) each copy-state based on number of segments (+ segment length)
    copy.states <- sapply(copy.states, function(each.copy.state){
      copy.idx <- (elementMetadata(gr.chr)[,1] == each.copy.state)
      cs.df <- data.frame("copy.state"=each.copy.state, "count"=sum(copy.idx, na.rm=TRUE))
      if(weight) cs.df <- data.frame("copy.state"=each.copy.state, "count"=sum(width(gr.chr[which(copy.idx), ])))
      cs.df
    })
    
    # Optional normalizing based on size of chromosome
    if(normalize && !(weight)) warning("It is STRONGLY advised not to normalize without weighing by the segment size first.")
    if(normalize) copy.states['count',] <- sapply(copy.states['count',], function(x) x / max(end(cytoband.chr)))
    
    colnames(copy.states) <- as.character(copy.states['copy.state',])
    copy.states <- copy.states[-grep('copy.state', rownames(copy.states)), ,drop=FALSE]
    as.data.frame(copy.states)
  })
  copystate <- rbind.fill(copystate)
  rownames(copystate) <- seqnames(gr)@values
  
  return(copystate)
}
