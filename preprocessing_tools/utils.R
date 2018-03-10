# All-purpose function to generate seg/interval files into a GRanges object
dataframeToGranges <- function(intdf, 
                               chr.regex="^chr(om|omosome)?$",
                               start.regex="^(loc|seg)?(.|_)?start(s|_idx)?$",
                               end.regex="^(loc|seg)?(.|_)?end(s|_idx)?$", ...){
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


# RLE function with an added functionality to report start and end indices
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