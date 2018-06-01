#----------------------------------------------------------------------------------------
#' cnTools: get Genes in TxDb.Hsapiens.UCSC.hg19.knownGene
#' @description Gets the genes from UCSC hg19 TxDb knownGene
#'
#' @return A Granges object containing strand-specific genes with EntrezIDs
#' @export
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene GenomicRanges
#' 
#' @examples getGenes()
getGenes <- function(genome.build="hg19"){
    switch(genome.build,
         hg19={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg19.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg19.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg19.knownGene
           },
         hg38={
           suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg38.knownGene))
           if(!exists("TxDb.Hsapiens.UCSC.hg38.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg38.knownGene")
           package <- TxDb.Hsapiens.UCSC.hg38.knownGene
           },
         stop("genome must be 'hg19' or 'hg38'"))
  
  genes0 <- genes(package)
  idx <- rep(seq_along(genes0), elementNROWS(genes0$gene_id))
  genes <- granges(genes0)[idx]
  genes$gene_id = unlist(genes0$gene_id)
  genes
}

#----------------------------------------------------------------------------------------
#' cnTools: annotate GRanges segments
#'
#' @param cn.data [Data.frame]: A dataframe that can be converted to GRanges object easily, or a granges object
#' @param genes [GRanges]: A GRanges object of genes with gene_ids housing annotation data. Easiest as the output from getGenes()
#'
#' @return Annotated GRanges object with gene ids for the input GRanges
#' @export
#' @import org.Hs.eg.db GenomicRanges
#' @importFrom AnnotationDbi mapIds
#' 
#' @examples 
#' annotateSegments(PLTK::genDemoData(), PLTK::getGenes())
annotateSegments <- function(cn.data, genes){
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  suppressPackageStartupMessages(require(GenomicRanges))
  suppressPackageStartupMessages(require(AnnotationDbi))
  gr0 <- cn.data
  if(class(cn.data) != 'GRanges') gr0 <- makeGRangesFromDataFrame(cn.data,keep.extra.columns=TRUE)
  if(class(cn.data) != 'GRanges') stop("Input data could not be converted to a GRanges object")
  
  olaps <- findOverlaps(genes, gr0, type="within")
  idx <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
  gr0$gene_ids <- splitAsList(genes$gene_id[queryHits(olaps)], idx)
  gr0$gene_ids <- lapply(gr0$gene_ids, function(input.id) {
    if(length(input.id) > 0){ 
      tryCatch({
        mapIds(org.Hs.eg.db,
               keys=input.id,
               column="SYMBOL",
               keytype="ENTREZID",
               multiVals="first")
      }, error=function(e){NULL})
    } else { NA }
  })
  return(gr0)
}

#----------------------------------------------------------------------------------------
#' cnTools: Aggregate Genomic Ranges
#' @description Because different samples have different segmentations, this function attempts to create a unified matrix with the fewest number of segments needed to represent all copy-number segments from all samples.  This allows for the elementMetadata() to store the associate copy-number values for those samples all in one matrix.
#'
#' @param list.gr [GRangesList]: a list of granges objects or GRangesList
#'
#' @return A single GRanges object with each row of the elementMetadata() containing all the samples copy-number values for that segment
#' @import GenomicRanges
#' @export
#'
#' @examples aggregateGr(list(gr1, gr2))
aggregateGr <- function(list.gr){
  # Loop to combine the first two GRanges object of the list, save it to the first element of the list, pop out the second element
  while(length(list.gr) >= 2){
    x <- list.gr[[1]]
    y <- list.gr[[2]]
    # Make a composite GRanges object that contains all possible segments
    int.gr <- disjoin(sort(c(granges(x), granges(y))))
    int.gr <- sort(c(int.gr, gaps(int.gr)))
    
    # Function to create an nrow(int.gr) matrix containing all the mapped elements from x and y
    xy.list <- lapply(list(x, y), function(z, int.gr){
      mat.blank <- matrix(nrow=length(int.gr),
                          ncol=1)
      olaps <- findOverlaps(int.gr, z, type = 'within')
      
      mat.fill <- apply(elementMetadata(z), 2, function(meta.z){
        mat.blank[queryHits(olaps), ] <- meta.z[subjectHits(olaps)]
        mat.blank
      })
      as.data.frame(mat.fill)
    }, int.gr=int.gr)
    
    # Combining the x-y elements back into the intersected GRanges object
    elementMetadata(int.gr) <- cbind(elementMetadata(int.gr), do.call("cbind", xy.list))
    
    list.gr[[2]] <- NULL
    list.gr[[1]] <- int.gr
  }
  list.gr[[1]]
}

#----------------------------------------------------------------------------------------
#' cnTools: Converts the .seg file to a GRanges object
#' @description Uses the standard seg inputs (https://software.broadinstitute.org/software/igv/sites/cancerinformatics.org.igv/files/linked_files/example.seg) and converts it into a list of GRanges object with the seg.mean column being the copy-number data stored in elementMetadata()
#'
#' @param seg A dataframe containing the seg file, with headers
#' @param col.id Name of sample [Default: SampleX]
#'
#' @return A Granges object with one column in the element Metadata() corresponding to the sample
#' @import GenomicRanges
#'
#' @examples
segfileToGr <- function(seg, col.id='SampleX'){
  suppressPackageStartupMessages(require(GenomicRanges))
  gr.tmp <- makeGRangesFromDataFrame(seg, keep.extra.columns=FALSE)
  elementMetadata(gr.tmp)$seg.mean <- seg$seg.mean
  colnames(elementMetadata(gr.tmp)) <- col.id
  gr.tmp
}

#----------------------------------------------------------------------------------------
#' cnTools: Convert Log2Ratio to Amp/Del
#' @description A temporary helper function to truncate log2ratios at 0.5 and turn them into Amp (1) or Del (-1)
#'
#' @param seg.gr [GRanges]:  A GRanges object with copy-number log2ratios in elementMetadata() columns
#' @param cn.thresh [Integer]: A log2ratio cutoff to indicate gain or loss
#' @param cn.scale [Numeric]: X scale 
#'
#' @return
#'
#' @examples
assignAmpDel <- function(seg.gr, cn.thresh=0.5, cn.scale=0){
  seg.gr.meta <- apply(elementMetadata(seg.gr), 2, function(x){
    del <- (-1 * cn.thresh) + cn.scale
    amp <- (1 * cn.thresh) + cn.scale
    amp.idx <-  x > amp
    del.idx <-  x < del
    neutral.idx <- (x >= del) & (x <= amp)
    x[amp.idx] <- 1
    x[del.idx] <- -1
    x[neutral.idx] <- 0
    x
  })
  elementMetadata(seg.gr) <- as.data.frame(seg.gr.meta)
  seg.gr
}
  
#----------------------------------------------------------------------------------------
#' cnTools: Convert to GRanges Wrapper
#' @description A wrapper to take any kind of copy-number data and convert it to a GRanges object with elementMetadata() storing a matrix of all copy-number values
#'
#' @param cnsegs Copy-number data: QDNAseq object, or .seg data frame
#' @param type [Character]: Specification of data type: "Qcalls", "Qsegmented", or "segfile"
#'
#' @return A Granges object with all samples combined into one singular matrix.  All copy-number values are stored in elementMetadata()
#' @import GenomicRanges QDNAseq
#' @export
#'
#' @examples
convertToGr <- function(cnsegs, type='Unknown'){
  if(class(cnsegs) == 'QDNAseqCopyNumbers'){
    suppressPackageStartupMessages(require(QDNAseq))
    gr <- makeGRangesFromDataFrame(cnsegs@featureData@data, keep.extra.columns=FALSE)
    if(type == 'Qcalls'){
      elementMetadata(gr) <- QDNAseq:::calls(cnsegs)
    } else if(type == 'Qsegmented'){
      elementMetadata(gr) <- QDNAseq:::segmented(cnsegs)
    }
    
  } else if(is.data.frame(cnsegs) && type == 'segfile'){
    suppressPackageStartupMessages(require(GenomicRanges))
    warning("SEGFILE: Assuming segfile standard outlined at https://software.broadinstitute.org/software/igv/sites/cancerinformatics.org.igv/files/linked_files/example.seg")
    cnsegs.id <- split(cnsegs, f=cnsegs$ID)
    
    # Convert segfiles dataframes to GRanges objects
    all.segfiles <- lapply(seq_along(cnsegs.id), function(seg.idx) {
      seg <- cnsegs.id[[seg.idx]]
      col.id <- names(cnsegs.id)[seg.idx]
      segfileToGr(seg, col.id)
    })
    
    # Aggregate all granges objects into one single GRanges
    if(length(all.segfiles) > 1) {
      segfile.gr <- aggregateGr(all.segfiles) 
    } else {
      segfile.gr <- all.segfiles[[1]]
    }
    
    # Classify as AMP or DEL based on log2 ratio cutoffs
    gr <- assignAmpDel(segfile.gr)
    gr
  } else {
    warning("Could not type of input segs")
  }
  
  seqlevelsStyle(gr) <- 'UCSC'
  seq.ids <- gsub("chr24", "chrY", gsub("chr23", "chrX", seqlevels(gr)))
  gr <- renameSeqlevels(gr,seq.ids)
  gr
}

#----------------------------------------------------------------------------------------
#' cnTools: Wrapper for copy-number metrics
#' @description A wrapper to run copy-number analysis on a GRanges copy-number dataset with copy-number values stored in elementMetadata().  These metrics include calculating genomic fraction and wGII scores (bp) for gains, losses, or any CN-abberation.
#'
#' @param gr [GRanges]:  GRanges object
#' @param cn.stat [Character]:  'all' for all CN-aberrations or 'gain' or 'loss'
#' @param copy.neutral [Integer]:  Integer specifying what a copy-neutral value is [default = 0]
#' @param cn.variance [Numeric]:  Numeric specifying how much range, if any, should be applied to copy.neutral[default = 0]
#' @param analysis [Character]:  The analysis to perform: "gf" genomic fraction, "wgii" for wGII scores (genomic fraction normalized for chromosome)
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
cnMetrics <- function(analysis=NA, gr=NULL, cn.stat='all', copy.neutral=0, cn.variance=0){
  #if(!validateGr(gr)) stop("Copy-number GRanges object failed validation checks.")
  switch(analysis,
         gf=PLTK:::cnGenomeFraction(analysis, gr, cn.stat, copy.neutral, cn.variance),
         wgii=PLTK:::cnGenomeFraction(analysis, gr, cn.stat, copy.neutral, cn.variance),
         stop("analysis not recognized")
  )
}

#----------------------------------------------------------------------------------------
#' cnTools: Wrapper for copy-number metrics
#' @description Calculates Genomic Fraction or wGII scores. Refer to cnMetrics for more detail.
#'
#' @param ... analysis, gr pass in from cnMetrics
#' @param analysis 
#' @param gr 
#' @param cn.stat 
#' @param copy.neutral 
#'
#' @return
#'
#' @examples
cnGenomeFraction <- function(analysis, gr, cn.stat='all', copy.neutral, cn.variance){
  # Gets copy-number breakdown of genome in basepairs (gains, losses, NA, etc)
  getGFdata <- function(each.cn, copy.neutral=0, ...){
    na.idx <- (is.na(each.cn))
    gain.idx <- (each.cn > (copy.neutral + cn.variance))
    loss.idx <- (each.cn < (copy.neutral - cn.variance))
    neutral.idx <- (each.cn >= (copy.neutral - cn.variance)) & (each.cn <= (copy.neutral + cn.variance))
    
    total.genome.size <- sum(as.numeric(width(gr)))
    non.na.genome.size <- sum(as.numeric(width(gr[which(!na.idx),])))
    gain.genome.size <- sum(as.numeric(width(gr[which(gain.idx),])))
    loss.genome.size <- sum(as.numeric(width(gr[which(loss.idx),])))
    neutral.genome.size <- sum(as.numeric(width(gr[which(neutral.idx),])))
    matrix(c("total"=total.genome.size,
             "non.na"=non.na.genome.size,
             "gain"=gain.genome.size,
             "loss"=loss.genome.size,
             "all"=(gain.genome.size + loss.genome.size),
             "neutral"=neutral.genome.size), 
           ncol=1)
  }

  # Cycles through each chromosome to get all the CN data for all samples
  chr.gf.data <- lapply(as.character(seqnames(gr)@values), function(chr.id){
    gr.chr <- gr[seqnames(gr) == chr.id]
    gf.chr <- apply(elementMetadata(gr.chr), 2, getGFdata)
    rownames(gf.chr) <- c("total", "non.na", "gain", "loss", "all", "neutral")
    gf.chr
  })
  names(chr.gf.data) <- as.character(seqnames(gr)@values)
  
  # Goes through all CN data to give back Genomic Fraction or wGII scores
  cn.row.idx <- grep(cn.stat, rownames(chr.gf.data[[1]]))
  total.row.idx <- grep("non.na", rownames(chr.gf.data[[1]]))
  switch(analysis,
         gf={
           cnstat.mat <- sapply(chr.gf.data, function(x) x[cn.row.idx,,drop=FALSE])
           total.mat <- sapply(chr.gf.data, function(x) x[total.row.idx,,drop=FALSE])
           gf.scores <- apply(cnstat.mat, 1, sum, na.rm=TRUE) / apply(total.mat, 1, sum, na.rm=TRUE)
           names(gf.scores) <- colnames(chr.gf.data[[1]])
           gf.scores
         },
         wgii={
           wgii.mat <- sapply(chr.gf.data, function(x) x[cn.row.idx,,drop=FALSE] / x[total.row.idx,,drop=FALSE])
           wgii.scores <- apply(wgii.mat, 1, mean, na.rm=TRUE)
           names(wgii.scores) <- colnames(chr.gf.data[[1]])
           wgii.scores
         },
         stop("analysis incorrectly specified"))

}

#----------------------------------------------------------------------------------------
#' cnTools: Wrapper for copy-number metrics
#' @description Takes an individual sample and reduces GRanges ranges where copy-number does not change between subsequent ranges
#'
#' @param gr [GRanges]: GRanges object with copy-number in elementMetadata()
#' @param sample.idx [Integer]: the sample index in the metadata
#' @param na.rm [Boolean]: Remove NA or keep NAs in the reduction
#'
#' @return
#' @export
#'
#' @examples
collapseSample <- function(gr, sample.idx, na.rm=TRUE){
  if(sample.idx > ncol(elementMetadata(gr))) stop("Please specify the index of a sample in elementMetadata()")
  if(is.null(sample.idx)) stop("Cannot collapse GRanges without the index of your sample in elementMetadta()")
  
  sample.id <- colnames(elementMetadata(gr))[sample.idx]
  print(sample.id)
  
  # Parse out CN information for the sample
  reduce.list <- lapply(as.character(seqnames(gr)@values), function(each.chr){
    gr.chr <- gr[seqnames(gr)==each.chr]
    each.sample <- elementMetadata(gr.chr)[,sample.idx]
    
    # Reduce continuous segments into a single segment based on rle
    rle.sample <- PLTK::getRleIdx(each.sample)
    gr.intervals <- matrix(c(start(gr.chr), end(gr.chr)),ncol=2)
    reduce.gr <- GRanges(seqnames=seqnames(gr.chr)[rle.sample$start.idx], 
                         ranges = IRanges(start=gr.intervals[rle.sample$start.idx, 1],
                                          end=gr.intervals[rle.sample$end.idx, 2]), 
                         strand=strand(gr.chr)[rle.sample$start.idx])
    if(any(!rle.sample$na.stat)) reduce.gr <- reduce.gr[-which(!rle.sample$na.stat),]

    # Reassign copy-number values to the reduced GRanges object
    #reduce.gr <-  Reduce(c, reduce.gr[!sapply(reduce.gr, is.null)])
    cn.values <- rle.sample$values
    if(na.rm) cn.values <- as.integer(na.omit(cn.values))
    elementMetadata(reduce.gr)$x <- cn.values
    colnames(elementMetadata(reduce.gr)) <- sample.id
    
    reduce.gr
  })
  Reduce(c, reduce.list)
}

#----------------------------------------------------------------------------------------
#' cnTools: Maps one GRanges objec to another (Intersect)
#'
#' @param gr [GRanges]: Input GRanges object with CNdata stored in elementMetadata()
#' @param ref.gr [GRanges]: Target GRanges object to map to
#' @param overlap [Character]: How to handle when multiple ranges from input_GRanges are within a single bin of the target_GRanges.  Options: "mode", "mean", "null" (all are weighted by segment sizes)
#' @param mode.type [Character]: "normal" or  "quick"
#'
#' @return The same Target GRanges object with elementMetadata() filled in
#' @export
#'
#' @examples
mapGrToReference <- function(gr, ref.gr, overlap='mode', mode.type='normal'){
  # Initial set-ups
  n <- ncol(elementMetadata(gr))
  sample.ids <- colnames(elementMetadata(gr))
  gr.st <- start(ranges(gr))
  gr.end <- end(ranges(gr))
  ref.st <- start(ranges(ref.gr))
  ref.end <- end(ranges(ref.gr))
  tmp.mat <- as.matrix(elementMetadata(gr)) # Speeds up analysis
  
  # Setting up overlaps
  olap <- findOverlaps(ref.gr, gr, type='any', select='all', ignore.strand=TRUE)
  qh.rle <- PLTK::getRleIdx(duplicated(queryHits(olap)))
  adj.cn.mat <- matrix(ncol=n, nrow=length(ref.gr), 
                       dimnames = list(NULL,sample.ids))
  simplifyRleMat <- function(rle.x, ana.type='single'){
    # Removes all instances where the start.idx and end.idx are the same since those are duplicates
    switch(ana.type,
           single= rle.x.type <- which(!as.logical(rle.x$values)),
           duplicate= rle.x.type <- which(as.logical(rle.x$values)))
    rle.x.mat <- (do.call("cbind.data.frame", rle.x))[rle.x.type,]
    rle.x.mat$start.idx <- as.numeric(rle.x.mat$start.idx)
    rle.x.mat$end.idx <- as.numeric(rle.x.mat$end.idx)
    
    switch(ana.type,
           single={
             diff.int <- (rle.x.mat$end.idx - rle.x.mat$start.idx)
             if(all(diff.int==0)) NA else rle.x.mat[diff.int!=0,]
           },
           duplicate=rle.x.mat)
  }  
  
  # SINGLE: Fill all reference to target GRanges 1:1 mappings
  qh.rle.single <- simplifyRleMat(qh.rle, ana.type='single')
  if(!is.na(qh.rle.single)) {
    apply(qh.rle.single, 1, function(each.s){
      olap.s <- olap[as.numeric(each.s['start.idx']):(as.numeric(each.s['end.idx'])-1),]
      adj.cn.mat[queryHits(olap.s),] <<- tmp.mat[subjectHits(olap.s),,drop=FALSE]
      return(NA)
    })
  }
  
  
  
  # Handles all instances where a reference/target range maps to multiple input ranges
  qh.rle.dup <- simplifyRleMat(qh.rle, ana.type='duplicate')
  apply(qh.rle.dup, 1, function(each.dup){
    print(each.dup)
    s.idx <- as.numeric(each.dup['start.idx'])
    e.idx <- as.numeric(each.dup['end.idx'])
    
    # Extract the query to subject duplicate mappings
    olap.dup <- olap[(s.idx-1):e.idx,]
    ref.idx <- unique(queryHits(olap.dup))
    ref.gr.dup <- ref.gr[ref.idx,]
    gr.dup <- tmp.mat[unique(subjectHits(olap.dup)), ,drop=FALSE]
    
    # Summarize the input gr elementMetadata to condense into the reference/target GR
    # Note: intersect() on GRanges and IRanges are both VERY lengthy and suboptimal for large datasets
    #int.dup <- sapply(seq_along(gr.dup), function(x) intersect(ref.gr.dup, gr.dup[x,]))
    int.dup <- interval_intersection(Intervals(c(ref.st[ref.idx],
                                                 ref.end[ref.idx])), 
                                     Intervals(c(gr.st[subjectHits(olap.dup)],
                                                 gr.end[subjectHits(olap.dup)])))
    int.width <- (int.dup[,2]-int.dup[,1])
    int.width[int.width==0] <- min(int.width[int.width!=0])
    
    switch(overlap,
           mode={
             adj.cn <- apply(gr.dup, 2, function(each.n){
               if(length(unique(each.n, na.rm=TRUE))>1){
                 if(mode.type=='normal'){
                   table.n <- table(na.omit(rep(each.n, 
                                                (int.width / min(int.width))))) 
                 } else if(mode.type=='quick'){
                   table.n <- table(na.omit(each.n))
                 }
                 mode.table <- as.numeric(names(table.n[table.n == max(table.n)]))
               } else {
                 mode.table <- unique(each.n, na.rm=TRUE)
               }
               mean(mode.table)
             })
             adj.cn <- t(as.matrix(adj.cn, nrow=1))
           },
           mean={
             adj.cn <- apply(gr.dup, 2, function(each.n){
               weighted.mean(each.n, int.width, na.rm = TRUE)
             })
             adj.cn <- t(as.matrix(adj.cn, nrow=1))
           },
           null={
             adj.cn <- matrix(rep(NA, n), nrow=1)
             colnames(adj.cn) <- colnames(elementMetadata(gr.dup))
           })
    adj.cn.mat[ref.idx,] <<-  adj.cn
  })
  elementMetadata(ref.gr) <- adj.cn.mat
  ref.gr
}

#----------------------------------------------------------------------------------------
#' cnTools: Maps one GRanges objec to another (Overlap)
#' @description A less senstive method compared to mapGrToReference() that maps one GRanges object (input) to a reference GRanges object (Target).  This works through the intervals package rather than the GRanges and IRanges to greatly increase the throughput.  It also operates using interval_overlaps rather than intersect to take a general "any overlap" method as compared to the aforementioned method
#'
#' @param gr [GRanges]: Input GRanges object containing copy-number elementMetadata()
#' @param ref.gr [GRanges]: Target GRanges object
#' @param int_ids [Characters]: The start and end column, shouldn't have to change this as its grabbed from getGrCoords() [default: c("start", "end")]
#' @param chr_id [Characters]: The chromosome id, shouldn't have to change this as its grabbed from getGrCoords() [default: "chr"]
#' @param overlap [Character]: The type of overlap summary metric: "mean", "median", "mode" (organized from quickest to slowest)
#'
#' @return A GRanges object of the Target GRanges object with the summarized elementMetadata() from input GRanges
#' @export
#' @import intervals
#'
#' @examples
overlapGrToReference <- function(gr, ref.gr,  int_ids=c("start", "end"), 
                                 chr_id='chr', overlap){
  gr.chrs <- as.character(seqnames(sort(ref.gr))@values)
  all.ref.gr.chr <- lapply(gr.chrs, function(each_chr){
    print(paste0("Running overlap ", each_chr, "..."))
    target.chr <- getGrCoords(ref.gr[seqnames(ref.gr) == each_chr])
    h_chr <- getGrCoords(gr[seqnames(gr) == each_chr])
    cn_mat <- as.matrix(elementMetadata(gr[seqnames(gr) == each_chr]))
    
    tryCatch({
      getIntervals <- function(xdf, colids){
        Intervals(as.matrix(xdf[,colids]))
      }
      overlap.idx <- interval_overlap(getIntervals(target.chr, int_ids), 
                                      getIntervals(h_chr, int_ids))
      
      # For all overlapping segments, calculates the mean, mode, or median
      getTargetVal <- function(cn_mat, overlap.idx, summfun){
        mean_cn_target <- sapply(overlap.idx, function(each_target){
          apply(cn_mat[each_target,,drop=FALSE], 2, summfun, na.rm=TRUE)
        })
        mean_cn_target <- as.data.frame(t(mean_cn_target))
        ### Developmental, may not work
        na_idx <- which(apply(mean_cn_target, 1, function(x) all(is.na(x))))
        if(length(na_idx) > 0) mean_cn_target <- mean_cn_target[-na_idx,]
        mean_cn_target
        ###
      }
      print(paste0("Summarizing overlap ", each_chr, "..."))
      switch(overlap,
             mean=summfun <- mean,
             mode=summfun <- getMode,
             median=summfun <- median)
      mean_cn_target <- getTargetVal(cn_mat, overlap.idx, summfun)
      print(paste0("Collapsing ", each_chr))
      
      
      #target.chr_idx <- sapply(overlap.idx, length) > 0  # reference dataset intervals that have a match
      ref.gr.chr <- ref.gr[seqnames(ref.gr) == each_chr] 
      elementMetadata(ref.gr.chr) <- as.matrix(mean_cn_target)
    }, error=function(e){
      ref.gr.chr <- ref.gr[seqnames(ref.gr) == each_chr] 
    })
    ref.gr.chr
  })
  
  samples.per.chr <- sapply(all.ref.gr.chr, function(x) ncol(x@elementMetadata))
  chr.rm.idx <- which(samples.per.chr != ncol(gr@elementMetadata))
  
  if(length(chr.rm.idx) != 0){
    warning(paste0(paste(gr.chrs[chr.rm.idx], collapse=","), " did not return matching samples"))
    all.ref.gr.chr <- all.ref.gr.chr[-chr.rm.idx]
  }
  
  
  Reduce(c, all.ref.gr.chr)
  
}

#----------------------------------------------------------------------------------------
#' cnTools: Get GRanges Coordinates
#'
#' @param keep.extra.columns [Boolean]: Whether to append the elementMetadata() or not
#' @param gr [GRanges]: Granges object
#'
#' @description Takes a GRanges object and converts the seqnames, range, and strand into a dataframe
#'
#' @return 4-column dataframe of "chr", "start", "end", "strand"
#' @export
#'
#' @examples getGrCoords(genDemoData())
getGrCoords <- function(gr, keep.extra.columns=FALSE){
  gr.df <- data.frame("chr"=rep(seqnames(gr)@values, 
                                seqnames(gr)@lengths),
                      "start"=start(gr),
                      "end"=end(gr),
                      "strand"=strand(gr))
  gr.df <- cbind(gr.df, as.matrix(elementMetadata(gr)))
  as.data.frame(gr.df)
}

cnDist <- function(gr, method="euclidean"){
  egr <- as.matrix(gr@elementMetadata)
  
  switch(method,
         "mhamming"={
           dist.mat <- apply(egr, 2, function(x){
             mat.x <- abs(matrix(rep(x, ncol(egr)), ncol=ncol(egr)) - egr)
             apply(mat.x, 2, sum, na.rm=TRUE)
           })
         },
         "pearson"={
           dist.mat <- round(cor(as.matrix(egr), use="pairwise.complete.obs"), 2)
           dist.mat <- dist.mat + 1
         },
         stop("Please provide a valid distance metric: 'euclidean', 'pearson', 'mhamming'")
         )
  dist.mat
}


