#install.packages("devtools")
library(devtools)

test.package=file.path("~/git", "PLTK/PLTK")
#devtools::create(test.package)
## Modify all package files here

devtools::document(test.package)
devtools::check(test.package)
devtools::build(test.package)

install.packages("~/git/PLTK/PLTK_0.0.0.9000.tar.gz",
                 repos = NULL, type="source")


library(PLTK)

demo <- genDemoData()
sigClusterBreakpoints(demo)