#----------------------------------------------------------------------------------------
#' cnTools: get Genes in TxDb.Hsapiens.UCSC.hg19.knownGene
#' @description Gets the genes from UCSC hg19 TxDb knownGene
#'
#' @return A Granges object containing strand-specific genes with EntrezIDs
#' @export
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene GenomicRanges
#' 
#' @examples getGenes()
getGenes <- function(){
  suppressPackageStartupMessages(require(TxDb.Hsapiens.UCSC.hg19.knownGene))
  suppressPackageStartupMessages(require(GenomicRanges))
  if(!exists("TxDb.Hsapiens.UCSC.hg19.knownGene")) stop("Requires TxDb.Hsapiens.UCSC.hg19.knownGene")
  genes0 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
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
#' @param [GRanges]: seg.gr A GRanges object with copy-number log2ratios in elementMetadata() columns
#' @param [Integer]: cn.thresh A log2ratio cutoff to indicate gain or loss
#'
#' @return
#'
#' @examples
assignAmpDel <- function(seg.gr, cn.thresh=0.5){
  seg.gr.meta <- apply(elementMetadata(seg.gr), 2, function(x){
    del <- (-1 * cn.thresh)
    amp <- (1 * cn.thresh)
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
#' @param [Character]: type Specification of data type: "QDNAseq", or "segfile"
#'
#' @return A Granges object with all samples combined into one singular matrix.  All copy-number values are stored in elementMetadata()
#' @import GenomicRanges QDNAseq
#' @export
#'
#' @examples
convertToGr <- function(cnsegs, type='Unknown'){
  if(class(cnsegs) == 'QDNAseqCopyNumbers' || type == 'QDNAseq'){
    suppressPackageStartupMessages(require(QDNAseq))
    gr <- makeGRangesFromDataFrame(cnsegs@featureData@data, keep.extra.columns=FALSE)
    elementMetadata(gr) <- QDNAseq:::calls(cnsegs)[,1:4]
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
#' @param [GRanges]: gr GRanges object
#' @param [Character]: cn.stat 'all' for all CN-aberrations or 'gain' or 'loss'
#' @param [Integer]: copy.neutral Integer specifying what a copy-neutral value is [default = 0]
#' @param [Character]: analysis The analysis to perform: "gf" genomic fraction, "wgii" for wGII scores (genomic fraction normalized for chromosome)
#'
#' @return
#' @export
#'
#' @examples
cnMetrics <- function(analysis=NA, gr, cn.stat='all', copy.neutral=0){
  #if(!validateGr(gr)) stop("Copy-number GRanges object failed validation checks.")
  
  switch(analysis,
         gf=cnGenomeFraction(analysis, ...),
         wgii=cnGenomeFraction(analysis, ...),
         stop("analysis not recognized")
  )
}

#----------------------------------------------------------------------------------------
#' cnTools: Wrapper for copy-number metrics
#' @description Calculates Genomic Fraction or wGII scores. Refer to cnMetrics for more detail.
#'
#' @param analysis [Character]: passed in from cnMetrics
#' @param gr [GRanges]: passed in from cnMetrics
#' @param [Character]: cn.stat passed in from cnMetrics
#' @param ... 
#'
#' @return
#'
#' @examples
cnGenomeFraction <- function(analysis, gr, cn.stat='all', ...){
  # Gets copy-number breakdown of genome in basepairs (gains, losses, NA, etc)
  getGFdata <- function(each.cn, copy.neutral=0, ...){
    na.idx <- (is.na(each.cn))
    gain.idx <- (each.cn > copy.neutral)
    loss.idx <- (each.cn < copy.neutral)
    neutral.idx <- (each.cn == copy.neutral)
    
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
  chr.gf.data <- lapply(seqnames(gr)@values, function(chr.id){
    gr.chr <- gr[which(seqnames(gr) == chr.id),]
    gf.chr <- apply(elementMetadata(gr.chr), 2, getGFdata)
    rownames(gf.chr) <- c("total", "non.na", "gain", "loss", "all", "neutral")
    gf.chr
  })
  names(chr.gf.data) <- seqnames(gr)@values
  
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
#'
#' @examples
collapseSample <- function(gr, sample.idx, na.rm=TRUE){
  if(sample.idx > ncol(elementMetadata(gr))) stop("Please specify the index of a sample in elementMetadata()")
  if(is.null(sample.idx)) stop("Cannot collapse GRanges without the index of your sample in elementMetadta()")
  
  # Parse out CN information for the sample
  each.sample <- elementMetadata(gr)[,sample.idx]
  sample.id <- colnames(elementMetadata(gr))[sample.idx]
  
  # Reduce continuous segments into a single segment based on rle
  rle.sample <- getRleIdx(each.sample)
  reduce.gr <- sapply(seq_along(rle.sample$start.idx), function(each_rle){
    if(!(na.rm && !(rle.sample$na.stat[each_rle]))){
      reduce(gr[rle.sample$start.idx[each_rle]:rle.sample$end.idx[each_rle],])
    } else if(!na.rm){
      reduce(gr[rle.sample$start.idx[each_rle]:rle.sample$end.idx[each_rle],])
    }
  })
  
  # Reassign copy-number values to the reduced GRanges object
  reduce.gr <-  Reduce(c, reduce.gr[!sapply(reduce.gr, is.null)])
  cn.values <- rle.sample$values
  if(na.rm) cn.values <- as.integer(na.omit(cn.values))
  elementMetadata(reduce.gr)$x <- cn.values
  colnames(elementMetadata(reduce.gr)) <- sample.id
  
  reduce.gr
}