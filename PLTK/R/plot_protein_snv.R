mafLayer <- function(maf, input.style='mutect'){
  maf.format <- switch(input.style,
                       mutect=.input_mutect(maf))
  
  # Reorder based on the CDS sequence
  lapply(maf.format, function(each.tx){
    pos <- gsub("\\+.*", "", each.tx$HGVSc) %>%
      gsub("[a-zA-Z\\.\\>]*", "", .) %>%
      strsplit(split="_")
    mins <- sapply(pos, function(i) min(as.integer(i)))
    each.tx[order(mins),]
  })
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
  mut.cds.idx <- apply(mut, 1, function(m){
    as.integer(strsplit(gsub("[a-zA-Z\\.]*", "", m['HGVSc']), "_")[[1]])
  })
  
  # check for failures
  null.idx <- sapply(seq.mat, is.null)
  if(any(null.idx)){
    error.hgvs <- paste(mut[which(null.idx), 'HGVSc'], collapse=", ")
    warning(paste0("Could not handle the following mutations: ", error.hgvs))
    seq.mat <- do.call(cbind, seq.mat[-which(null.idx)])
    mut.cds.idx <- mut.cds.idx[-which(null.idx)]
  }
  
  seq.mat <- cbind(seq, seq.mat)
  seq.mat <- .aggregate_seq(seq.mat)
  
  
  list('mat'=seq.mat, 'idx'=mut.cds.idx)
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


plotMuts <- function(mut.mat, mut.idx, cdna.range, prot.coords){
  aa.start <- floor((cdna.range[1] - 1)/3)
  aa.end <- ceiling((cdna.range[2] - 1)/3)
  cds.start <- (aa.start * 3) + 1
  cds.end <- (aa.end * 3)
  
  cds.idx <- c(cds.start:cds.end)
  mut.idx[[length(mut.idx)+1]] <- mut.idx
  cds.aa.mut <- lapply(c(1:ncol(mut.mat)), function(e.mut){
    mm.cds.idx <- if(e.mut > 1)  mut.idx[e.mut-1] else NULL
    list("aa"=.formatDNA(mut.mat[cds.idx,e.mut], TRUE),
         "cds"=.formatDNA(mut.mat[cds.idx,e.mut]),
         "cds.idx"=.idxDNA(mut.mat[cds.idx, c(1, e.mut)], 
                           mm.cds.idx, cds.idx))
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
    mm.stat <- cds.aa[['cds.idx']]
    
    plot(0, type='n', axes=FALSE, xlab='', ylab='', xaxt='n', yaxt='n',
         xlim=c(0, nchar(cds.aa.mut[[1]]$aa)), ylim=c(0.2,1))
    if(t.idx==1){
      par(xpd=TRUE)
      ratio <- .conv_coord(old.x=prot.coords[['xlim']], new.x=c(0, length(aa)))
      segments(x0 = prot.coords[['st']]*ratio,y0 = 1.17,x1 = 0,y1 = 1.03)
      segments(x0 = prot.coords[['end']]*ratio,y0 = 1.17,x1 = length(aa),y1 = 1.03)
      par(xpd=FALSE)
    }
    
    trigger.fs <- FALSE
    for(i in seq_along(aa)){
      ## Fill in un-altered stuff
      rect(xleft = i-1, ybottom = 0.2, xright = i, ytop = 1, 
           col=if (trigger.fs) "bisque" else if(i %% 2 ==1) 'lightblue' else 'aliceblue')
      text(x = i-0.5, y=0.6, labels=aa[i])
      text(x = i-0.5, y=0.9, labels=cds[i], cex=0.7)
      text(x = i-0.5, y=0.3, labels=c((aa.start+1):(aa.end+1))[i], cex=0.7)
      
      ## Fill in altered
      if(t.idx!=1) trigger.fs <- .fill_errors(i, mm.stat, trigger.fs,
                                              cds, aa)
    }
  })
 
}

.fill_errors <- function(i, mm.stat, trig=trigger.fs, cds, aa){
  # Checks if any of the AA residues have a mut in them
  if(any(i == sapply(mm.stat, function(i) i$aa))){
    #flags which mutation is being addressed
    mut.selection <- which(i == sapply(mm.stat, function(i) i$aa))
    t.shift <- sum(sapply(mm.stat[1:mut.selection], function(i) i$shift))
    mm.stat <- mm.stat[[mut.selection]]
    
    rect(xleft = i-1, ybottom = 0.2, xright = i, ytop = 1, 
         col='indianred1')
    text(x = i-0.5, y=0.9, labels=cds[i], cex=0.7)
    
    if(any(duplicated(mm.stat$cds))){
      idx <- unique((mm.stat$cds - 1) %% 3)
      lbl <- if(idx==1) "^ " else if(idx==2) " ^"
      text(x = i-0.5, y=0.9, labels=lbl, 
           cex=0.7, col="indianred4", pos = 1)
    }
    text(x = i-0.5, y=0.6, labels=aa[i], col='indianred4')
    text(x = i-0.5, y=0.3, labels=c((aa.start+1):(aa.end+1))[i], cex=0.7)
    
    trig <- if((t.shift %% 3) != 0) TRUE else FALSE
  }
  trig
}

.idxDNA <- function(seq.mat, mm.cds.idx, cds.idx){
  #seq.mat <- mut.mat[cds.idx, c(1, e.mut)]
  if(!is.null(mm.cds.idx)){
    if(class(mm.cds.idx[[1]]) == 'list') mm.cds.idx <- mm.cds.idx[[1]]
    lapply(mm.cds.idx, function(mm){
      mm.idx <- suppressWarnings(which(mm == cds.idx))
      
      # Identify frameshifts
      ref.seq <- paste(seq.mat[mm.idx,1], collapse="")
      mm.seq <- paste(seq.mat[mm.idx,2], collapse="") %>% gsub("-", "", .)
      mm.len <- nchar(mm.seq) - nchar(ref.seq)
      
      cds.idx <- c(min(mm.idx), min(mm.idx) + nchar(mm.seq)) # Which CDS bases are affected
      aa.idx <- unique(ceiling(mm.idx/3)) # Which AA residue is affected?
      fs.stat <- (mm.len %% 3) != 0 # Is there a frameshift?
      shift <- (mm.len %% 3) # How much of a shift
      
      list("cds"=cds.idx,
           "aa"=aa.idx,
           "fs"=fs.stat,
           "shift"=shift)
    })
  }
  
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

