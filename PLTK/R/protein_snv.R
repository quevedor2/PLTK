#' .getDomainBimap
#' @description Maps protein feature IDs to protein feature domain names
#'  (i.e.  PF02260 -> FATC)
#' @param domain.src 
#'
#' @return
#'
#' @examples .getDomainBimap('pfam')
.getDomainBimap <- function(domain.src){
  switch(domain.src,
         pfam={
           require(PFAM.db)
           db <- PFAM.db::PFAMID
           keys <- mappedkeys(db)
           as.list(db[keys])
         },
         test={
           stop("This does nothing")
          })
}



#' getProteinDat
#' @description Retrieves the protein structure and features from the ensembldb database.
#' 
#' @param gene Gene of interest in HUGO format (e.g., MTOR)
#' @param edb Ensembl DB R package, v86=GRCh38, v75=GRCh37 (Default: EnsDb.Hsapiens.v86)
#' @param txid Primary transcript ID analyzing. Returns first TXID if none provided (Default: NULL)
#' @param domain.src Protein feature databse, only pfam supported at this time (Default: pfam)
#' @param default.longest Selects the longest ENST txid as the canonical txid (Default: TRUE)
#'
#' @return
#' @export
#'
#' @examples getProteinDat("MTOR", EnsDb.Hsapiens.v86, BSgenome.Hsapiens.NCBI.GRCh38)
getProteinDat <- function(gene, edb, bsg, txid=NULL, 
                          domain.src='pfam', default.longest=TRUE){
  require(ensembldb)
  require(EnsDb.Hsapiens.v86)
  require(EnsDb.Hsapiens.v75)
  require(BSgenome.Hsapiens.NCBI.GRCh38)
  require(BSgenome.Hsapiens.UCSC.hg19)
  
  pd <- proteins(edb, filter = GeneNameFilter(gene),
                 columns = c("tx_id", "protein_id", "protein_domain_source", 
                             "protein_domain_id", "prot_dom_start", "prot_dom_end"),
                 return.type = "AAStringSet")
  transcripts(edb, filter = GeneNameFilter(gene))
  mpd <- mcols(pd) 
  spl.pd <- split(mpd, mpd$tx_id)
  
  # Gets the protein domains for each isoform
  splsrc.pd <- lapply(spl.pd, function(pd0, ...) {
    pd.src <- split(pd0, pd0$protein_domain_source)[[domain.src]]
    domain <-.getDomainBimap(...)[pd.src$protein_domain_id]
    pd.src$domain <- unlist(domain)
    return(pd.src)
  }, domain.src)
  
  # Retrieves the CDS sequence for each isoform
  CDS <- ensembldb::cdsBy(edb, by="tx", TxIdFilter(names(splsrc.pd)),
                          columns = c("tx_biotype", "gene_name"))
  CDS_seqs <- extractTranscriptSeqs(bsg, CDS)
  
  if(is.null(txid) & !default.longest){
    txid <- 1
  } else if (is.null(txid) & default.longest){
    CDS.len <- sapply(split(CDS, f=names(CDS)), function(i) as.integer(sum(width(i))))
    enst.idx <- grep("ENST", names(CDS.len))
    CDS.len <- CDS.len[enst.idx]
    txid <- names(which.max(CDS.len[enst.idx]))
  }
  
  
  
  pid <- unique(splsrc.pd[[txid]]$protein_id)
  aa <- as.vector(unique(pd[which(names(pd) == pid),]))
  cds <- as.vector(CDS_seqs[txid])
  
  list("aa"=aa,
       'cds'=as.vector(cds),
       "pdomains"=as.data.frame(splsrc.pd[[txid]]))
}

#' chooseTranscript
#' @description Select a single transcript using a gene name or mutation locus selection method
#' 
#' @param edb Ensembl DB R package, v86=GRCh38, v75=GRCh37 (Default: EnsDb.Hsapiens.v86)
#' @param gene Gene of interest in HUGO format (e.g., MTOR). If present, Locus infomation is not required.
#' @param Chromosome Chromosome
#' @param Start_Position Mutation start position
#' @param End_Position Mutation end position
#' @param default.longest Selects the longest ENST txid as the canonical txid (Default: TRUE)
#'
#' @return
#' @export
#'
#' @examples chooseTranscript(EnsDb.Hsapiens.v86,"MTOR")
chooseTranscript <- function(edb, gene, Chromosome, Start_Position, End_Position,
                             default.longest = TRUE){
  require(ensembldb)
  require(EnsDb.Hsapiens.v86)
  require(EnsDb.Hsapiens.v75)
  
  seqlevelsStyle(edb) <- "UCSC"
  
  if(!is.null(gene)){
    trans <- transcripts(edb, GeneNameFilter(gene),
                            columns = c("tx_biotype", "gene_name"))
  }else{
    
    Chromosome <- as.character(Chromosome)
    Start_Position <- as.numeric(as.character(Start_Position))
    End_Position <- as.numeric(as.character(End_Position))
    
    mutrange <-  GRanges(Chromosome, IRanges(start=Start_Position, end=End_Position))
    
    trans <- transcriptsByOverlaps(edb, mutrange, type = "any")
  }
  
    trans_vec <- trans@elementMetadata@listData[["tx_id"]]
    
    CDS.len <- sapply(trans@elementMetadata@listData[["tx_id"]], 
                      function(i) {
                        CDS <- ensembldb::cdsBy(edb, by="tx", TxIdFilter(t),
                                    columns = c("tx_biotype", "gene_name"))
                        sum(CDS@unlistData@ranges@width)})
    
    if(!default.longest){
      txid <- trans_vec[1]
    } else {

      CDS.len <- sapply(split(CDS, f=names(CDS)), function(i) as.integer(sum(width(i))))
      
      enst.idx <- grep("ENST", names(CDS.len))
      CDS.len <- CDS.len[enst.idx]
      txid <- names(which.max(CDS.len[enst.idx]))
    }
  }
  return(txid)
}

#' Title
#'
#' @param ids 
#' @param col.set 
#'
#' @return
#'
#' @examples 
.mapColors <- function(ids, col.set='Set1'){
  require(RColorBrewer)
  cols <- setNames(brewer.pal(length(ids), col.set), ids)
  cols
}

#' plotProteinStructure
#' @description Plots the top protein structure data track
#'
#' @param pdat Output list data structure from getProteinDat()
#' @param cex.val CEX value for the writing
#' @param st.range Start amino acid index for region of interest
#' @param end.range End amino acid index for region of interest
#'
#' @return
#' @export
#'
#' @examples
#' pdf("~/Desktop/test2.pdf", width=15)
#' split.screen(c(3,1))
#' screen(1)
#' plotProteinStructure(pdat)
#' close.screen(all.screens=TRUE)
plotProteinStructure <- function(pdat, st.range, end.range, cex.val=0.7){
  seqlen <- nchar(pdat$aa)

  cols <- .mapColors(unique(pdat$pdomains$domain))
  pdat$pdomains$cols <- cols[pdat$pdomains$domain]
  
  aa.idx <- seq(0, 5000, by=500)[-1]
  aa.idx <- c(1, aa.idx[aa.idx < seqlen], seqlen)
  
  # Plot the base structure of the protein
  par(xpd=FALSE)
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xaxt='n', yaxt='n',
       xlim=c(1, seqlen), ylim=c(0.2,1), main=unique(pdat$pdomains$gene_name))
  rect(xleft = 1, ybottom = 0, xright = seqlen, ytop = 1, col='grey')
  axis(side=1, at = aa.idx, cex.axis=cex.val)
  
  # Add rectangles for each pfam domain
  apply(pdat$pdomains, 1, function(dom){
    st <- as.integer(dom['prot_dom_start'])
    end <- as.integer(dom['prot_dom_end'])
    mid <- st + ((end - st)/2)
      
    rect(xleft=st, ybottom=0, xright=end, ytop=1, 
         col=dom['cols'], border=NA)
    axis(side = 3, at = mid, labels = dom['domain'], 
         line = NA, tick = FALSE, cex.axis=cex.val)
  })
  
  ## Add regions outside of margins
  par(xpd=TRUE)
  text(x=-1 * (seqlen*0.01), y=0.015, "Amino acid #", pos = 2, cex=cex.val)
  segments(x0 = st.range,y0 = 0.17,x1 = st.range,y1 = -100)
  segments(x0 = end.range,y0 = 0.17,x1 = end.range,y1 = -100)
  par(xpd=FALSE)
  
  list("xlim"=c(1, seqlen), "st"=st.range, "end"=end.range)
}


#' getCdsDat
#' @description For a given transcript, return the coding sequence
#' 
#' @param txid Transcript ID in Ensembl format
#' @param edb Ensembl DB R package, v86=GRCh38, v75=GRCh37 (Default: EnsDb.Hsapiens.v86)
#' @param gdb genome sequence database (Default: BSgenome.Hsapiens.UCSC.hg38)
#' 
#' @return
#' @export
#'
#' @examples getCdsDat("ENST00000471181", EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38)
getCdsDat <- function(txid, edb, gdb){
  require(ensembldb)
  require(EnsDb.Hsapiens.v86)
  require(EnsDb.Hsapiens.v75)
  require(BSgenome.Hsapiens.UCSC.hg38)
  require(BSgenome.Hsapiens.UCSC.hg19)

  
  CDS.frame <- as.data.frame(ensembldb::cdsBy(edb, by="tx", TxIdFilter(txid),
                            columns = c("tx_biotype", "gene_name")))
  

  CDS.frame$Reference_Sequence <- sapply(row.names(CDS.frame), function(i){
    gr <-  GRanges(as.character(CDS.frame[i,"seqnames"]), 
                   IRanges(start=as.numeric(CDS.frame[i,"start"]), 
                           end=as.numeric(CDS.frame[i,"end"])))
    as.character(getSeq(gdb, gr))})
  
  
  cds_start <- min(c(CDS.frame$start,CDS.frame$end))
  cds_end <- max(c(CDS.frame$start,CDS.frame$end))
  
  Indexed_Seq <- vector(mode = "character",
                        length = cds_end-cds_start+1)
  
  if(all(CDS.frame$strand=="-")){
  
    CDS.frame$Reference_Sequence <- sapply(CDS.frame$Reference_Sequence,
                                           function(i){
                                             as.character(reverseComplement(DNAString(i)))
                                           }
                                     )
  
  for (i in 1:nrow(CDS.frame)){
      Indexed_Seq[c(1+cds_end - as.numeric(CDS.frame[i,"end"]):
                    as.numeric(CDS.frame[i,"start"]))] <- 
      reverse(unlist(strsplit(CDS.frame[i,"Reference_Sequence"],"")))
      names(Indexed_Seq) <- cds_end:cds_start
  }
  }else{
    for (i in 1:nrow(CDS.frame)){
      Indexed_Seq[c(as.numeric(CDS.frame[i,"start"]):
                      as.numeric(CDS.frame[i,"end"])-cds_start+1)] <- 
        unlist(strsplit(CDS.frame[i,"Reference_Sequence"],""))
    }
    names(Indexed_Seq) <- cds_start:cds_end
  }
  
  AA_seq <- translate(DNAString(paste(CDS.frame$Reference_Sequence,collapse="")))
  
  return(list(cds = CDS.frame,
              Indexed_Seq = Indexed_Seq,
              AA_seq = AA_seq))
}
  
#' mutateCds
#' @description Given a CDS data frame and a mutation, return the Amino-Acid
#' sequence, aligned to the WT sequence
#' 
#' @param Indexed_Seq An Indexed and named vector of CDS sequences
#' Indeces are relative to transcription start, names are according to genomic
#' position.
#' @param Chromosome Chromosome
#' @param Start_Position Mutation start position
#' @param End_Position Mutation end position
#' @param Variant_Type Variant type, options: c("SNP","INS","DEL")
#' @param Alt_Allele Alternate allele to be added
#'
#' @return
#' @export
#'
#' @examples mutateCds("ENST00000471181", EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38)
mutateCds <- function(Indexed_Seq, Start_Position, End_Position, 
                      Variant_Type, Alt_Allele){
  
  mut_Indexed_Seq <- Indexed_Seq
  
  if(Variant_Type == "SNP")
  {
    mut_Indexed_Seq[as.character(Start_Position)] <- Alt_Allele
  }else if(Variant_Type == "DEL")
  {
    mut_Indexed_Seq[as.character(Start_Position:End_Position)] <- "-"
  }else if(Variant_Type == "INS")
  {
    ins_ind <- which(names(mut_Indexed_Seq)==as.character(Start_Position))
    mut_Indexed_Seq <- c(mut_Indexed_Seq[1:ins_ind], 
                     strsplit(Alt_Allele,""),
                     mut_Indexed_Seq[End_Position:length(mut_Indexed_Seq)])
  }
  
  AA.Seq <- translate(DNAString(paste0(mut_Indexed_Seq[mut_Indexed_Seq!="" & 
                                   mut_Indexed_Seq!="-"],collapse="")))
  
  return(list(Indexed_Seq = mut_Indexed_Seq,
              AA_seq = AA_seq))
}


