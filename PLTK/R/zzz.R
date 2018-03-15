.onLoad <- function(...){
  packageStartupMessage("\tPughLab ToolKit [PLTK]: This package relies heavily on GenomicRanges objects and it is best practice to use make use of the makeGRangesFromDataFrame() function from GenomicRanges before inputting segment files into copy-number functions\n")
}
