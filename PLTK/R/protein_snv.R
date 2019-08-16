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
#'
#' @return
#' @export
#'
#' @examples getProteinDat("MTOR", EnsDb.Hsapiens.v86)
getProteinDat <- function(gene, edb, txid=NULL, domain.src='pfam'){
  require(ensembldb)
  require(EnsDb.Hsapiens.v86)
  require(EnsDb.Hsapiens.v75)
  
  pd <- proteins(edb, filter = GeneNameFilter(gene),
                 columns = c("tx_id", "protein_id", "protein_domain_source", 
                             "protein_domain_id", "prot_dom_start", "prot_dom_end"),
                 return.type = "AAStringSet")
  mpd <- mcols(pd) 
  spl.pd <- split(mpd, mpd$tx_id)
  
  splsrc.pd <- lapply(spl.pd, function(pd0, ...) {
    pd.src <- split(pd0, pd0$protein_domain_source)[[domain.src]]
    domain <-.getDomainBimap(...)[pd.src$protein_domain_id]
    pd.src$domain <- unlist(domain)
    return(pd.src)
  }, domain.src)
  
  if(is.null(txid)) txid <- 1
  pid <- unique(splsrc.pd[[txid]]$protein_id)
  seq <- as.vector(unique(pd[which(names(pd) == pid),]))
  
  list("seq"=seq,
       "pdomains"=as.data.frame(splsrc.pd[[txid]]))
}

getProteinDat("MTOR", EnsDb.Hsapiens.v86)



