maf.file <- "~/Desktop/tmp/sample_BRCA1_reversion.maf.txt"
maf <- read.table(maf.file, header = TRUE, sep="\t", quote='',
                  stringsAsFactors = FALSE, check.names = FALSE, 
                  fill=FALSE, comment.char = "")
muts <- mafLayer(maf)
genes <- sapply(muts, function(x) unique(x$Hugo_Symbol))

gene <- genes[1]
mut <- muts[[1]]
pd <- getProteinDat(genes[1], 
                    edb=EnsDb.Hsapiens.v86, 
                    bsg=BSgenome.Hsapiens.NCBI.GRCh38,
                    txid = names(genes)[1])

mut.mat <- mutSeq(mut = muts[[1]], pd = pd)

plotMuts(mut.mat, cdna.range=c(2990, 3015))

mafLayer <- function(maf, input.style='mutect'){
  switch(input.style,
         mutect=.input_mutect(maf))
}

.input_mutect <- function(maf){
  col.ids <- c('Hugo_Symbol', 'NCBI_Build', 'Chromosome', 
               'Start_Position', 'End_Position', 'Variant_Type',
               "HGVSc", "HGVSp_Short", "Transcript_ID", 
               "t_alt_count", "t_depth", "AF")
  maf$AF <- round(with(maf, t_alt_count/t_depth), 2)
  maf.red <- maf[,col.ids]
  
  # Validation check of data: ensure same BRCA1 gene and NCBI build
  
  split(maf.red[,c('Hugo_Symbol', 'Variant_Type', 'HGVSc', 'HGVSp_Short', 
                   'Transcript_ID', 'AF')], f=maf$Transcript_ID)
}



mutSeq <- function(mut, pd){
  require(Biostrings)
  seq <- strsplit(pd$cds, '')[[1]]
  seq.mat <- apply(mut, 1, function(m){
    switch(m['Variant_Type'],
           DEL=.hgsv_del(seq, m['HGVSc']))
  })
  
  seq.mat <- cbind(seq, seq.mat)
  seq.mat <- .aggregate_seq(seq.mat)
  
  
  seq.mat
}

.aggregate_seq <- function(seq.mat){
  # First column should be the reference sequence
  # Each column should be one mutation from cHGVS
  change.idx <- apply(seq.mat[,-1], 2, function(each.seq){
    each.seq != seq.mat[,1] 
  })
  
  # Goes through and modifies the reference seq one mut at a time
  seq.agg <- seq.mat[,1]
  for(cidx in c(1:ncol(change.idx))){
    cidx.pos <- which(change.idx[,cidx])
    seq.agg[cidx.pos] <- seq.mat[,cidx+1][cidx.pos]
  }
  cbind(seq.mat, seq.agg)
}

.hgsv_del <- function(seq, hgvs){
  pos <- gsub("^c.", "", hgvs) %>% gsub("[a-zA-Z]*$", "", .)
  pos <- strsplit(pos, split="_")[[1]]
  if(length(pos) > 1){
    seq[as.integer(pos[1]):as.integer(pos[2])] <- '-'
  } else {
    seq[as.integer(pos)] <- '-'
  }
  seq
}



cdna.range=c(2990, 3015)
cdna.range <- c(1, 9)
plotMuts <- function(mut.mat, cdna.range){
  aa.start <- floor((cdna.range[1] - 1)/3) 
  aa.end <- ceiling((cdna.range[2] - 1)/3)
  
  cds.start <- (aa.start * 3) + 1
  cds.end <- (aa.end * 3)
  
  cds.idx <- c(cds.start:cds.end)
  lapply(c(1:ncol(mut.mat)), function(e.mut){
    t(data.frame("aa"=.formatDNA(mut.mat[cds.idx,e.mut], TRUE),
               "cds"=.formatDNA(mut.mat[cds.idx,e.mut])))
  })
}

.formatDNA <- function(seq, translate.seq=FALSE){
  if(any(seq == '-')) seq <- seq[-which(seq=='-')]
  
  if(translate.seq){
    suppressWarnings(as.character(translate(DNAString(paste(seq, collapse="")))))
  } else {
    paste(seq, collapse="")
  }
}
