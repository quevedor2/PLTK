# Generates the genomic bins of 10kb using QDNAseq algorithms to obtain mappability and GC content as well

library(QDNAseq)
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg19)
#### FUNCTIONS
createBins2 <- function (bsgenome, binSize, ignoreMitochondria = TRUE, 
                         excludeSeqnames = NULL, offset=NULL) {
  chrs <- GenomeInfoDb::seqnames(bsgenome)
  try({
    info <- GenomeInfoDb::genomeStyles(GenomeInfoDb::organism(bsgenome))
    style <- GenomeInfoDb::seqlevelsStyle(bsgenome)[1]
    chrs <- info[, style]
  }, silent = TRUE)
  if (!is.null(excludeSeqnames)) {
    chrs <- chrs[!chrs %in% excludeSeqnames]
  }
  if (ignoreMitochondria) {
    selectedMT <- grep("^(chr)?M(T)?$", chrs)
    if (length(selectedMT) != 0L) 
      chrs <- chrs[-selectedMT]
  }
  lengths <- GenomeInfoDb::seqlengths(bsgenome)[chrs]
  QDNAseq:::vmsg("Creating bins of ", binSize, " kbp for genome ", substitute(bsgenome))
  binWidth <- as.integer(binSize * 1000L)
  chrData <- QDNAseq:::flapply(chrs, function(chr) {
    QDNAseq:::vmsg("    Processing ", chr, " ...")
    chr.size <- lengths[chr]
    chr.starts <- seq(from = 1L, to = chr.size, by = binWidth)
    chr.ends <- chr.starts + binWidth - 1L
    if(!is.null(offset)){
      chr.starts <- chr.starts + (offset * 1000L)
      chr.ends <- chr.ends + (offset * 1000L)
    }
    chr.ends[length(chr.ends)] <- chr.size
    chr.seq <- BSgenome::getSeq(bsgenome, chr, as.character = TRUE)
    bin.seq <- substring(chr.seq, first = chr.starts, last = chr.ends)
    acgt <- gsub("[^ACGT]", "", bin.seq)
    cg <- gsub("[^CG]", "", acgt)
    chr.bases <- nchar(acgt)/(binWidth) * 100
    chr.gc <- nchar(cg)/nchar(acgt) * 100
    list(start = chr.starts, end = chr.ends, bases = chr.bases, 
         gc = chr.gc)
  })
  start <- unlist(lapply(chrData, "[[", "start"))
  end <- unlist(lapply(chrData, "[[", "end"))
  bases <- unlist(lapply(chrData, "[[", "bases"))
  gc <- unlist(lapply(chrData, "[[", "gc"))
  gc[is.nan(gc)] <- NA_real_
  bins <- data.frame(chromosome = rep(chrs, times = ceiling(lengths/binWidth)), 
                     start, end, bases, gc, stringsAsFactors = FALSE)
  bins$chromosome <- sub("^chr", "", bins$chromosome)
  rownames(bins) <- sprintf("%s:%i-%i", bins$chromosome, bins$start, 
                            bins$end)
  bins
}

#### MAIN
# set the bin size
binSize <- 10000
# create bins from the reference genome
setwd("/mnt/work1/users/pughlab/projects/shallowWGS_benchmark/reference/1000kb_stagger")
for(offset in 0){
  print(paste0("Generating an ", binSize, "kb bin with offset of ", offset))
  bins <- createBins2(bsgenome=BSgenome.Hsapiens.UCSC.hg19, binSize=binSize, offset=offset)
  
  map_file='/mnt/work1/users/pughlab/references/mappability/wgEncodeCrgMapabilityAlign50mer.bigWig'
  wig_over_bed='/mnt/work1/users/pughlab/references/mappability/bigWigAverageOverBed'
  # calculate mappabilites per bin from ENCODE mapability tracks
  mappability <- calculateMappability(bins,
                                      bigWigFile=map_file,
                                      bigWigAverageOverBed=wig_over_bed)
  if(length(mappability) != nrow(bins)){
    cat(paste0("\n\tSegmentation fault (core dumped) error while using UCSC bigWigAverageOverBed tool for ", offset, "\n"))
    next()
  }
  # mapbed <- '/mnt/work1/users/home2/quever/map800.bed'
  # map <- read.table(mapbed, sep = "\t", as.is = TRUE)
  # map <- map[order(map$V4), ]
  # mappability <- map$V5 * 100
  bins$mappability <- mappability
  
  encode_bl='/mnt/work1/users/pughlab/references/mappability/wgEncodeDacMapabilityConsensusExcludable.bed'
  encode_bl2='/mnt/work1/users/pughlab/references/mappability/wgEncodeDukeMapabilityRegionsExcludable.bed'
  # calculate overlap with ENCODE blacklisted regions
  bins$blacklist <- calculateBlacklist(bins,
                                       bedFiles=c(encode_bl, encode_bl2))
  
  # generic calculation of overlap with blacklisted regions
  # bins$blacklist <- calculateBlacklistByRegions(bins,
  #                                               cbind(chromosome, bpStart, bpEnd))
  bins$use <- bins$chromosome %in% as.character(1:22) & bins$bases > 0
  # convert to AnnotatedDataFrame and add metadata
  filt.idx <- which(bins$end < bins$start)
  if(length(filt.idx) > 0L){
    bins <- bins[-filt.idx,]
  }
  bins <- AnnotatedDataFrame(bins,
                             varMetadata=data.frame(labelDescription=c(
                               'Chromosome name',
                               'Base pair start position',
                               'Base pair end position',
                               'Percentage of non-N nucleotides (of full bin size)',
                               'Percentage of C and G nucleotides (of non-N nucleotides)',
                               'Average mappability of 50mers with a maximum of 2 mismatches',
                               'Percent overlap with ENCODE blacklisted regions',
                               'Whether the bin should be used in subsequent analysis steps'),
                               row.names=colnames(bins)))
  attr(bins, "QDNAseq") <- list(
    author="Ilari Scheinin",
    date=Sys.time(),
    organism="Hsapiens",
    build="hg19",
    version=packageVersion("QDNAseq"),
    md5=digest::digest(bins@data),
    sessionInfo=sessionInfo())
  bins <- makeGRangesFromDataFrame(bins@data, keep.extra.columns = TRUE)
  seqlevelsStyle(bins) <- 'UCSC'
  
  save(bins, file = paste0(binSize, "kb.", offset, ".Rda"))
  save(bins, file = file.path("annotated_df", paste0(binSize, "kb.", offset, ".Rda")))
}
