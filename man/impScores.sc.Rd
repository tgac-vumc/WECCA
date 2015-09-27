\name{impScores.sc}
\alias{impScores.sc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Importance scores }
\description{
Given the clusters, the importance score of each feature is calculated.
This is defined as the maximum pairwise symmetrized
Kullback Leibler divergence.
These scores reflect which features are most discriminative between clusters.
}
\usage{
impScores.sc(cghdata.regioned, dendrogram, nclusters)
}
\arguments{
  \item{cghdata.regioned}{ A \code{list}-object as returned by the function \code{regioning.sc}. }
  \item{dendrogram}{ Determines if and how the row dendrogram should be computed and reordered. Should be a dendrogram as returned by \code{WECCAsc}. }
  \item{nclusters}{ An integer with the desired number of clusters. }
}
\value{
  A \code{matrix} whose first five columns contain annotation information and the sixth the region's importance score.
}
\references{ Insert ref to article where importance scores are calculated. }
\author{ Wessel N. van Wieringen: \email{wvanwie@few.vu.nl} }
\seealso{ \code{KLdiv}, \code{regioning}, \code{WECCAsc} }
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
}
