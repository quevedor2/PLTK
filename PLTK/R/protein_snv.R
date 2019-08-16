if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("EnsDb.Hsapiens.v75")

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)

edb <- EnsDb.Hsapiens.v86
gene <- 'MTOR'

getProteinDat <- function(gene, edb){
  
}

pd <- proteins(edb, filter = GeneNameFilter(gene),
               columns = c("tx_id", "protein_domain_id", "prot_dom_start", "prot_dom_end"),
               return.type = "AAStringSet")
xpd <- split(mcols(pd), mcols(pd)$tx_id)
names(xpd)
xpd[[1]][order(xpd[[1]]$prot_dom_start),]
