% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotMethods.R
\name{PlotOverlapHeatmap}
\alias{PlotOverlapHeatmap}
\title{Plot of frequencies by cluster}
\usage{
PlotOverlapHeatmap(CDR3PairMat, ident, compute = "geomean")
}
\arguments{
\item{ident}{It can be cluster or celltype}

\item{compute}{geomean or logFC}
}
\description{
Plot of frequencies by cluster
}
\examples{
\dontrun{
sset_NKT = PlotTCR(sset, group.by = )
}
}
