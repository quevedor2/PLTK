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
sigSegSize(demo, normalize=FALSE)

example <- "/mnt/work1/users/home2/quever/example.seg"
example <- "~/Desktop/example.seg"
seg <- read.table(example, header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE)

gr1 <- convertToGr(copyNumbersCalled)
gr2 <- convertToGr(seg, type='segfile')
gr <- aggregateGr(list(gr1, gr2))
gr <- gr2
sapply(c("gain", "loss", "all"), function(x) cnMetrics(analysis='wgii', gr=gr, cn.stat=x, copy.neutral=0))
sapply(c("gain", "loss", "all"), function(x) cnMetrics(analysis='gf', gr=gr, cn.stat=x, copy.neutral=0))
all.sigs <- runCnSignatures(gr=gr, binsize=50000, bins=PLTK::bins)


