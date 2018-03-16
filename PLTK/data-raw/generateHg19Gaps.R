library(GenomicRanges)

# File downloaded from UCSC Table Browser:  hg19, All Tables > gap
gap.file  <- '~/Desktop/ucsc.hg19.table_gap.centromere_telomere.tsv'
gap.table <- read.table(gap.file, header=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(gap.table) <- c("bin", "chrom", "chromStart", "chromEnd", "ix",  "n", "size", "type", "bridge")


hg19.centromeres <- makeGRangesFromDataFrame(gap.table[which(gap.table$type == 'centromere'),], keep.extra.columns = TRUE)
seqlevelsStyle(hg19.centromeres) <- 'UCSC'
hg19.centromeres <- sort(hg19.centromeres)
save(hg19.centromeres, file="~/git/PLTK/PLTK/data/hg19.centromeres.RData")

hg19.telomeres <- makeGRangesFromDataFrame(gap.table[which(gap.table$type == 'telomere'),], keep.extra.columns = TRUE)
seqlevelsStyle(hg19.telomeres) <- 'UCSC'
hg19.telomeres <- sort(hg19.telomeres)
save(hg19.telomeres, file="~/git/PLTK/PLTK/data/hg19.telomeres.RData")


# File downloaded from UCSC Table Browser:  hg19, All Tables > cytoBand
band.file  <- '~/Desktop/ucsc.hg19.table_cytoBand.chromosome_band.tsv'
band.table <- read.table(band.file, header=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(band.table) <- c("chrom", "chromStart", "chromEnd", "name",  "gieStain")

hg19.cytobands <- makeGRangesFromDataFrame(band.table, keep.extra.columns = TRUE)
seqlevelsStyle(hg19.cytobands) <- 'UCSC'
hg19.cytobands <- sort(hg19.cytobands)
save(hg19.cytobands, file="~/git/PLTK/PLTK/data/hg19.cytobands.RData")
