# gistic.file <- 'scores.gistic'
# out.file <- paste0("ACC_gistic.pdf")
gisticIdeogram <- function(gistic.file, out.file){
  require(Gviz)
  require(GenomicRanges)
  gistic <- read.table(gistic.file, header=T, stringsAsFactors = F, 
                       check.names = F, sep="\t")
  gr <- makeGRangesFromDataFrame(gistic, keep.extra.columns = T)
  grl <- split(gr, gr$Type)
  
  .combineGr <- function(gr){
    gr.chr <- split(gr, seqnames(gr))
    gr.c <- as(lapply(gr.chr, function(gr0){
      pos <- unique(sort(c(start(gr0), end(gr0))))
      GRanges(seqnames = rep(unique(as.character(seqnames(gr0))), length(pos)-1),
              IRanges(start=pos[-length(pos)],
                      end=pos[-1]))
    }), 'GRangesList')
    gr.c <- unlist(gr.c)
    gr.c
  }
  .populateGr <- function(grl, ref.gr, col.id){
    for(each.id in names(grl)){
      ov.idx <- findOverlaps(grl[[each.id]], ref.gr)
      mcols(ref.gr)[,each.id] <- NA
      mcols(ref.gr)[subjectHits(ov.idx),each.id] <- mcols(grl[[each.id]])[queryHits(ov.idx),col.id]
    }
    ref.gr
  }
  
  gr.c.raw <- .combineGr(gr)
  pdf(out.file, ...)
  lapply(c('frequency', 'average amplitude', 'G-score', '-log10(q-value)'), 
         function(metric){
    gr.c <- .populateGr(grl, gr.c.raw, metric)
    
    dT <- DataTrack(gr.c, name="ID")
    if(metric == 'average amplitude') y.max <- 1.7 else y.max <- 1
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(10,1)))
    for(chr in paste0("chr", c(1:10))){
      itrack <- IdeogramTrack(genome='hg19', chromosome=chr)
      idx <- as.numeric(gsub("^chr", "", chr))
      pushViewport(viewport(layout.pos.col=1, 
                            layout.pos.row=idx))
      plotTracks(list(itrack, dT), groups=rownames(values(dT)), type=c("a"), 
                 showSampleNames=T, chromosome=c(chr), ylim=c(0,y.max), showId=F, add=T, legend=F)
      popViewport(1)
    }
  })
  dev.off()
}


