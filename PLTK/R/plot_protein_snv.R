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
           DEL=.hgsv_del(seq, m['HGVSc']),
           INS=.hgvs_ins(seq, m['HGVSc']),
           SNP=.hgvs_snv(seq, m['HGVSc']))
  })
  
  # check for failures
  null.idx <- sapply(seq.mat, is.null)
  if(any(null.idx)){
    error.hgvs <- paste(mut[which(null.idx), 'HGVSc'], collapse=", ")
    warning(paste0("Could not handle the following mutations: ", error.hgvs))
    seq.mat <- do.call(cbind, seq.mat[-which(null.idx)])
  }
  
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

.hgvs_del <- function(seq, hgvs){
  # hgvs <- "c.3008_3009del"
  pos <- gsub("^c.", "", hgvs) %>% gsub("[a-zA-Z]*$", "", .)
  pos <- strsplit(pos, split="_")[[1]]
  if(length(pos) > 1){
    seq[as.integer(pos[1]):as.integer(pos[2])] <- '-'
  } else {
    seq[as.integer(pos)] <- '-'
  }
  seq
}

.hgvs_ins <- function(seq, hgvs){
  # hgvs <- 'c.1982_1983insCT'
  # hgvs <- 'c.2560_2561dup'
  # hgvs <- 'c.2560dup'
  if(grepl("dup", hgvs)){
    pos <- gsub("^c.", "", hgvs) %>% gsub("dup$", "", .)  %>% strsplit(split="_")
    pos <- as.integer(pos[[1]])
    ins <- paste(seq[pos], collapse="")
    pos <- pos[2]
  } else if(grepl("ins", hgvs)) {
    pos <- gsub("^c.", "", hgvs) %>% gsub("ins.*$", "", .)  %>% strsplit(split="_")
    pos <- as.integer(pos[[1]])
    ins <- gsub("^c.[0-9_]*ins", "", hgvs) 
  }
  seq[pos[1]] <- paste0(seq[pos[1]], ins)
  
  seq
}

.hgvs_snv <- function(seq, hgvs){
  # hgvs <- "c.355G>A"
  # hgvs <- 'c.5137+1G>T'
  if(grepl("\\+", hgvs)){
    # Handle splice variants, not implemented yet
    pos <- 1
    ref <- 'X'
    alt <- 'Y'
  } else{
    pos <- as.integer(gsub("^c.", "", hgvs) %>% gsub("[^0-9]*$", "", .))
    ref <- gsub("^c.[0-9]*", "", hgvs) %>% gsub(">[a-zA-Z]*$", "", ., perl=TRUE)
    alt <- gsub("^c.[0-9]+[a-zA-Z]*>", "", hgvs) 
  }
  
  
  if(seq[pos] == toupper(ref)) {
    seq[pos] <- toupper(alt) 
  } else {
    warning(paste0(hgvs, ": Reference allele did not match reference sequence allele"))
    seq <- NULL
  }
  seq
}


plotMuts <- function(mut.mat, cdna.range, prot.coords){
  aa.start <- floor((cdna.range[1] - 1)/3) 
  aa.end <- ceiling((cdna.range[2] - 1)/3)
  cds.start <- (aa.start * 3) + 1
  cds.end <- (aa.end * 3)
  
  cds.idx <- c(cds.start:cds.end)
  cds.aa.mut <- lapply(c(1:ncol(mut.mat)), function(e.mut){
    list("aa"=.formatDNA(mut.mat[cds.idx,e.mut], TRUE),
         "cds"=.formatDNA(mut.mat[cds.idx,e.mut]),
         "cds.idx"=.idxDNA(mut.mat[cds.idx, c(1, e.mut)]))
  })

  scr.idx <- split.screen(c(length(cds.aa.mut),1))
  sapply(seq_along(cds.aa.mut), function(t.idx){
    screen(scr.idx[t.idx])
    par(mar=c(1, 4.1, 1, 2.1))
    cds.aa <- cds.aa.mut[[t.idx]]
    aa <- strsplit(cds.aa[['aa']], '')[[1]]
    cds <- substring(cds.aa[['cds']], 
                     seq(1, nchar(cds.aa[['cds']]), 3), 
                     seq(3, nchar(cds.aa[['cds']]), 3))
    
    plot(0, type='n', axes=FALSE, xlab='', ylab='', xaxt='n', yaxt='n',
         xlim=c(0, length(aa)), ylim=c(0.2,1))
    if(t.idx==1){
      par(xpd=TRUE)
      ratio <- .conv_coord(old.x=prot.coords[['xlim']], new.x=c(0, length(aa)))
      segments(x0 = prot.coords[['st']]*ratio,y0 = 1.17,x1 = 0,y1 = 1.03)
      segments(x0 = prot.coords[['end']]*ratio,y0 = 1.17,x1 = length(aa),y1 = 1.03)
      par(xpd=FALSE)
    }
    
    sapply(seq_along(aa), function(i) {
      rect(xleft = i-1, ybottom = 0.2, xright = i, ytop = 1, 
           col=if (i %% 2 ==1) 'lightblue' else 'aliceblue')
      text(x = i-0.5, y=0.6, labels=aa[i])
      text(x = i-0.5, y=0.9, labels=cds[i], cex=0.7)
      text(x = i-0.5, y=0.3, labels=c(aa.start:aa.end)[i], cex=0.7)
    })
  })
 
}

.idxDNA <- function(seq.mat){
  seq.mat <- mut.mat[cds.idx, c(1, e.mut)]
  mm.idx <- which(seq.mat[,2] != seq.mat[,1])
  aa.idx <- unique(ceiling(mm.idx/3))
  cds.idx <- mm.idx %% 3 %>% gsub("0", "3", .) %>% as.integer
}

.formatDNA <- function(seq, translate.seq=FALSE){
  if(any(seq == '-')) seq <- seq[-which(seq=='-')]
  
  if(translate.seq){
    suppressWarnings(as.character(translate(DNAString(paste(seq, collapse="")))))
  } else {
    paste(seq, collapse="")
  }
}

.conv_coord <- function(old.x, new.x){
  old.x <- old.x - min(old.x)
  new.x <- new.x - min(new.x)
  
  ratio <- max(new.x) / max(old.x)
  ratio
}

