#install.packages("devtools")
library(devtools)

test.package=file.path("~/git", "PLTK/PLTK")
#devtools::create(test.package)
## Modify all package files here

devtools::document(test.package)
devtools::check(test.package)

#devtools::use_vignette("pltk-vignette", test.package)
#devtools::build_vignettes(test.package)
devtools::build(test.package)


install.packages("~/git/PLTK/PLTK_0.0.0.9000.tar.gz",
                 repos = NULL, type="source")
library(PLTK)


demo <- genDemoData()
sigClusterBreakpoints(demo, 50)
sigBinBreakpoints(demo, PLTK::bins)
sigGapDist(demo, gap.type = "telomeres", gap = PLTK::hg19.telomeres)
sigGapDist(demo, gap.type = "telomeres", gap = PLTK::hg19.telomeres, normalize=TRUE)
sigGapDist(demo, gap.type = "centromeres", gap = PLTK::hg19.centromeres)
sigGapDist(demo, gap.type = "centromeres", gap = PLTK::hg19.centromeres, normalize=TRUE)
sigSegSize(demo)
sigSegSize(demo, normalize=TRUE)
sigCnChangepoint(demo, collapse.segs = TRUE)
sigCnChangepoint(demo, collapse.segs = FALSE)


example.expr <- "~/Desktop/zscore_BRCA.txt"
expr.df <- read.table(example.expr, header=TRUE,
                      check.names = FALSE, stringsAsFactors = FALSE)
gr.expr <- sort(makeGRangesFromDataFrame(expr.df, keep.extra.columns = TRUE))


example <- "/mnt/work1/users/home2/quever/example.seg"
example <- "~/Desktop/example.seg"
seg <- read.table(example, header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)

#gr1 <- convertToGr(copyNumbersCalled)
gr2 <- convertToGr(seg, type='segfile')
#gr <- aggregateGr(list(gr1, gr2))
gr.cn <- gr2
sapply(c("gain", "loss", "all"), function(x) cnMetrics(analysis='wgii', gr=gr, cn.stat=x, copy.neutral=0))
sapply(c("gain", "loss", "all"), function(x) cnMetrics(analysis='gf', gr=gr, cn.stat=x, copy.neutral=0))
all.sigs <- runCnSignatures(gr=gr.cn, binsize=50000, bins=PLTK::bins, 
                            assign.amp.del = FALSE, cn.thresh=0.5, cn.scale=2,
                            numeric.return=TRUE)
summarizeSignatures(all.sigs, ids=colnames(elementMetadata(gr)))



PDIR='/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets'
tad <- file.path(PDIR, '/reference/TAD/TAD_IMR90.txt')
tad <- read.table(tad, sep="\t", header=FALSE, col.names = c("chr", "start", "end"))
tad.gr <- makeGRangesFromDataFrame(tad)
seqlevelsStyle(tad.gr) <- 'UCSC'
mapped.ref.gr <- mapGrToReference(gr, tad.gr, overlap='mode')

t.chr <- 'chr1'
split.screen(c(3, 1))
screen(3)
par(mar=c(5.1, 4.1, 0.5, 2.1))
plot.settings <- initializeGrPlot(PLTK::hg19.cytobands, plot.chrom=TRUE,
                                  plot.cband=TRUE, alpha.factor=3, label.side='bottom',
                                 target.chr=t.chr)

screen(2)
par(mar=c(0.5, 4.1, 0.5, 2.1))
suppressWarnings(plotGrMetadata(gr.cn, plot.settings, col.ids=c(1,2), data.type='cn', 
                                add.axis=FALSE, axis.mark=4, chr.lines=TRUE,
                                target.chr=t.chr))


screen(1)
par(mar=c(0.5, 4.1, 4.1, 2.1))
suppressWarnings(plotGrMetadata(gr.expr, plot.settings, data.type='expr', 
                                yrange=c(0,5), side=3,
                                add.axis=TRUE, axis.mark=4, add.y.axis=TRUE,
                                anno.track=0.25, add.annotations=TRUE,
                                target.chr=t.chr))
close.screen(all.screens=TRUE)




