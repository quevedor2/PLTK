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
#' @examples getProteinDat("MTOR", EnsDb.Hsapiens.v86)
getProteinDat <- function(gene, edb, txid=NULL, 
                          domain.src='pfam', default.longest=TRUE){
  require(ensembldb)
  require(EnsDb.Hsapiens.v86)
  require(EnsDb.Hsapiens.v75)
  
  pd <- proteins(edb, filter = GeneNameFilter(gene),
                 columns = c("tx_id", "protein_id", "protein_domain_source", 
                             "protein_domain_id", "prot_dom_start", "prot_dom_end"),
                 return.type = "AAStringSet")
  transcripts(edb, filter = GeneNameFilter(gene))
  mpd <- mcols(pd) 
  spl.pd <- split(mpd, mpd$tx_id)
  
  splsrc.pd <- lapply(spl.pd, function(pd0, ...) {
    pd.src <- split(pd0, pd0$protein_domain_source)[[domain.src]]
    domain <-.getDomainBimap(...)[pd.src$protein_domain_id]
    pd.src$domain <- unlist(domain)
    return(pd.src)
  }, domain.src)
  
  if(is.null(txid) & !default.longest){
    txid <- 1
  } else if (default.longest){
    txlen <- lengthOf(edb, of = "tx", filter = TxIdFilter(names(splsrc.pd)))
    enst.idx <- grep("ENST", names(txlen))
    txlen <- txlen[enst.idx]
    txid <- names(which.max(txlen[enst.idx]))
  }
  
  pid <- unique(splsrc.pd[[txid]]$protein_id)
  seq <- as.vector(unique(pd[which(names(pd) == pid),]))
  
  list("seq"=seq,
       "pdomains"=as.data.frame(splsrc.pd[[txid]]))
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
  seqlen <- nchar(pdat$seq)

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
}




