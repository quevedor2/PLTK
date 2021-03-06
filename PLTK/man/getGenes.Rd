% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnTools.R
\name{getGenes}
\alias{getGenes}
\title{cnTools: get Genes in TxDb.Hsapiens.UCSC.hg19.knownGene}
\usage{
getGenes(genome.build = "hg19")
}
\value{
A Granges object containing strand-specific genes with EntrezIDs
}
\description{
Gets the genes from UCSC hg19 TxDb knownGene
}
\examples{
getGenes()
}
