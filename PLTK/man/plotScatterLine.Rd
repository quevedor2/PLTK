% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plotScatterLine}
\alias{plotScatterLine}
\title{Visualization: plotScatterLine}
\usage{
plotScatterLine(groups, targ.pnt = NA, top.anno = 10, y.lim = NA,
  withViolin = F, withBox = F, vio.col = "lightblue", connect = T,
  show.all = T, ...)
}
\arguments{
\item{groups}{a labeled matrix, dataframe or list containing multiple groups (cols)}

\item{targ.pnt}{Individual targets based on rownames or labelled vectors in list}

\item{top.anno}{Annotate the top X differences between groups [default=10]}

\item{y.lim}{provide only if specific min/max required c(min,max)}

\item{withViolin}{boolean if violin contour to be drawn}

\item{withBox}{boolean if box contour to be drawn}

\item{vio.col}{vector containing custom colors for violin fill}

\item{connect}{boolean if dots should be connected}

\item{show.all}{Boolean to show all connecting segments if targ.pnt is provided}

\item{...}{}
}
\description{
Creates a stripchart, boxplot, or violin plot for each group and then offers the option to link changes between individual elements between groups. If rownames or a labelled vector in a list is provided, it will annotate the top X differences between groups based on taking the max difference for a given element across all groups.
}
\examples{
demo.plotScatterLine()
}
