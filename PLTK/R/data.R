#' 10Mb Genomic Bins (hg19)
#'
#' GRanges object of hg19 segmented into 10mb bins using QDNAseq
#'
#' @format A GRanges object in UCSC format; 322 ranges spanning 24 chromosomes; 5 columns of metadata
#' \describe{
#'   \item{bases}{QDNAseq: number of bases}
#'   \item{gc}{QDNAseq: GC fraction of bin}
#'   \item{mappability}{QDNAseq: Mappability and Alignability track of bin}
#'   \item{blacklist}{QDNAseq: Score related on whether to blacklist the region}
#'   ...
#' }
"bins"