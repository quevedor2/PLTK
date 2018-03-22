#----------------------------------------------------------------------------------------
#' utils: genDemoData
#' @description Generates a random granges object to be used for demo purposes when trying out the cnSignature functions
#'
#' @param seed Set the random seed [default: 404]
#'
#' @return
#' @export
#' @import GenomicRanges
#' 
#' @examples 
#' demo <- genDemoData(404)
genDemoData <- function(seed=404){
  require(GenomicRanges)
  
  set.seed(seed=seed)
  #Generate random data
  all.segs <- c(runif(n = 10, min = 0, max = 1000),
                rnorm(10, mean=1200, sd=100),
                runif(n = 10, min = 0, max = 1000),
                runif(n = 3, min = 1001, max = 10000000),
                runif(n = 10, min = 10000001, max = 20000000),
                runif(n = 40, min = 20000001, max = 30000000))
  
  mkSegs <- function(data, chrid){
    ends <- sample(x = ceiling(data), size = length(data) * 0.7, replace = FALSE)
    ends <- sort(as.integer(ends))
    starts <- c(1, ends[-length(ends)]+1)
    intervals <- data.frame("chr"=rep(chrid, length(starts)),
                            "start"=starts,
                            "ends"=ends)
    intervals
  }
  intervals <- do.call("rbind", lapply(c("chr1", "chr2", "chr3"), function(x) mkSegs(all.segs, x)))
  
  #Dataframe to Granges object
  intervals.gr <- makeGRangesFromDataFrame(intervals, 
                                           ignore.strand = TRUE, 
                                           seqnames.field = "chr",
                                           start.field = "start", 
                                           end.field = "ends")
  
  intervals.gr
}


#----------------------------------------------------------------------------------------
#' utils: dataframeToGranges
#' @description OUTDATED: Uses a set of complex regex expression to automatically parse the chromosome, start_idx, and end_idx from a dataframe of segments to create a GRanges object. Also accepts standard GRanges arguments that you would normally pass in to define extra information.
#'
#' @param intdf An dataframe containing segment intervals
#' @param chr.regex Regex for the chromosome column
#' @param start.regex Regex for the start index column
#' @param end.regex Regex for the end index column
#' @param ... 
#'
#' @return
#' @import GenomicRanges
#' 
#' @examples
dataframeToGranges <- function(intdf, 
                               chr.regex="^chr(om|omosome)?$",
                               start.regex="^(loc|seg)?(.|_)?start(s|_idx)?$",
                               end.regex="^(loc|seg)?(.|_)?end(s|_idx)?$", ...){
  warning("Outdated function, use GenomicRanges::makeGRangesFromDataFrame() function")
  chr.idx <- grep(pattern=chr.regex, x = colnames(intdf), perl=TRUE, ignore.case = TRUE)
  start.idx <- grep(pattern=start.regex, x = colnames(intdf), perl=TRUE, ignore.case = TRUE)
  end.idx <- grep(pattern=end.regex, x = colnames(intdf), perl=TRUE, ignore.case = TRUE)
  
  if(all(sapply(list(chr.idx, start.idx, end.idx), length))){
    print(paste0("Chromosome column header: ", colnames(intdf)[chr.idx]))
    print(paste0("Seg start header: ", colnames(intdf)[start.idx]))
    print(paste0("Seg end header: ", colnames(intdf)[end.idx]))
    
    GRanges(seqnames=intdf[,chr.idx],
            ranges=IRanges(intdf[,start.idx], 
                           intdf[,end.idx], 
                           names=rownames(intdf)),
            strand=rep("+", nrow(intdf)), ...)
  } else {
    warning("Could not locate all chromosome, start, and end indices headers.  Specify your own chr.regex, start.regex, or end.regex where needed and try again")
    warning("Please raise an issue indicating what your headers were on the git repo so I can build it into the regex")
  }
}


#----------------------------------------------------------------------------------------
#' utils: getRleIdx
#' @description Similar to the rle() function, but gives start, end and NA indices on top of the existing "length" and "values". Easily viewed using str() on the return object().
#'
#' @param x A sequence to perform RLE on.  Can be either a vector or data frame
#' @param col.id If dataframe is provided as \code{x}, a column index is required
#' @param na.val Default value to set NAs too so RLE indexing doesn't break
#'
#' @return
#' @export
#'
#' @examples
#' getRleIdx(c(rep(1,5), rep(2,5), rep(3,5)))
getRleIdx <- function(x, col.id=NA, na.val=-100){
  if(is.vector(x)){
    reformat.na.x <- as.character(x)
    
  } else if(is.data.frame(x)){
    # Handles multiple columns
    if(length(col.id) > 1) {
      uniq.id <- apply(x, 1, function(y) paste(y[col.id], collapse="-"))
      x$uniq <- uniq.id
      col.id <- 'uniq'
    }
    reformat.na.x <- as.character(x[,col.id])
  }
  
  # Fill in NA values
  reformat.na.x[which(is.na(reformat.na.x))] <- na.val
  rle.x <- rle(reformat.na.x)
  
  #Get the array index for the start-to-end of each unique value/changepoint
  rle.x$start.idx <- c(1, (cumsum(rle.x$lengths) + 1)[-length(rle.x$lengths)])
  rle.x$end.idx <- rle.x$start.idx + (rle.x$lengths - 1)
  rle.x$values[which(rle.x$values == na.val)] <- NA
  rle.x$na.stat <- !is.na(rle.x$values) 
  
  return(rle.x)
}


#----------------------------------------------------------------------------------------
#' utils: Calculates the mode
#'
#' @param n [Vector]: A vector of numbers elements to calculate the mode
#'
#' @return 
#' @export
#'
#' @examples getMode(c(rep(1,5), 1:10))
getMode <- function(n, na.rm=FALSE){
  if(na.rm) table.n <- table(na.omit(n)) else table.n <- table(n)
  mode.table <- as.numeric(names(table.n[table.n == max(table.n)]))
  mode.table
}