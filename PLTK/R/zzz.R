.onLoad <- function(...){
  packageStartupMessage("The PLTK package heavily uses GRanges. It is best practice to use make use of the makeGRangesFromDataFrame() function from GenomicRanges\n")
}
