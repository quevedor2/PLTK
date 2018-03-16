# PLTK
Pughlab Toolkit (PLTK): Repository for all small analysis, plotting, preprocessing, and formatting scripts.  Visit the git wiki for a more detailed explanation of usage, purpose, and requirements for each tool/script

The package can be installed locally by forking or cloning the git repo to your local.  Then run the following installation commands:
```
install.packages("~/git/PLTK/PLTK_0.0.0.9000.tar.gz",
                 repos = NULL, type="source")
```

## visualization_tools
  * **plotScatterLine()**: an R function to combine violin plots, boxplots, and stripcharts between multiple groups and visually link elements that span all groups  (Cindy/Rene)
  * **plotLikelihoodRatio()**: an R function to plot the kernel density estimate of two datasets and calculate the likelihood ratio for a point estimate  (Rene)
## analysis_tools
  * **CNsignature**: A list of R functions meant to work with the GRanges object that will output CN-based signatures, similar to a somatic mutation signature
  * **CNtools**: A list of R functions meant to work with the GRanges object to easily import copy-number data and calculate a range of different metrics and analysis such as genomic-fraction altered and wGII scores.
  
## formatting_tools

## preprocessing_tools
  * **utils.R**: Contains a lot of little functions to help mediate tedious tasks or add extra functionalities to existing functions
  
## reference_files

