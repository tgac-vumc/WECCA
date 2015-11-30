\name{sample.cluster.table}
\alias{sample.cluster.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sample-cluster table }
\description{
  Cuts a the \code{WECCAsc} or \code{WECCAhc} generated dendrogram at a desired number of clusters. 
  Then, produces a \code{data.frame} telling which samples belong to which cluster.
}
\usage{
sample.cluster.table(cghdata.regioned, dendrogram, nclusters)
}
\arguments{
  \item{cghdata.regioned}{ A \code{list}-object as returned by \code{regioning.sc}. }
  \item{dendrogram}{ Determines if and how the row dendrogram should be computed and reordered. Should be a dendrogram as returned by \code{WECCAsc} or \code{WECCAhc}. }
  \item{nclusters}{ An integer with the desired number of clusters. }
}
\value{
  A \code{data.frame} with the first two columns containing the sample number and name, and the third column the cluster information.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{\link{cutree}} }
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
}
