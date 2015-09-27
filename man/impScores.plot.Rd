\name{impScores.plot}
\alias{impScores.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Importance score plot }
\description{
  Plots importance scores against genomic order.
}
\usage{
impScores.plot(impScores.table, yaxis.max = 10, lastchr = 22, color = "red")
}
\arguments{
  \item{impScores.table}{ Table containing region annotation information and importance scores as returned by the function \code{impScores.sc} or \code{impScores.hc}. }
  \item{yaxis.max}{ Maximum importance score to be plotted. }
  \item{lastchr}{ Integer corresponding to chromosome number. This is the last chromosome whose importance scores are included in the plot.}
  \item{color}{ Color of line connecting importance scores. }
}
\author{ Wessel N. van Wieringen: \email{wvanwie@few.vu.nl} }
\examples{
# generate object of class
data(WiltingCalled)

# make region data (soft and hard calls)
WiltingRegioned <- regioning(WiltingCalled)

# clustering with soft.calls
dendrogram <- WECCAsc(WiltingRegioned)

# generate a heatmap of the found clustering
WECCA.heatmap(WiltingRegioned, dendrogram)

# specify the number of clusters to be extracted from the dendrogram
nclusters <- 2
table.clusters.samples <- sample.cluster.table(WiltingRegioned, dendrogram, nclusters)

# calculate importance scores for each feature
impScores.table <- impScores.sc(WiltingRegioned, dendrogram, nclusters)

# plot importance scores
impScores.plot(impScores.table, yaxis.max=10, lastchr=22, color="red")
}
