\name{WECCA.heatmap}
\alias{WECCA.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Heatmap plotting with a chromosome indicator bar added. }
\description{
 A heatmap generation of the WECCA clustering with chromosome info color-coded.
}
\usage{
WECCA.heatmap(cghdata.regioned, dendrogram, ...)
}
\arguments{
  \item{cghdata.regioned}{ A \code{list}-object as returned by the function \code{regioning}. }
  \item{dendrogram}{ Determines if and how the row dendrogram should be computed and reordered. Should be a dendrogram as returned by the function \code{WECCAsc} or \code{WECCAhc}. }
}
\details{
  This function is merely a thin wrapper to the heatmap function. It adds a color bar with chromosome information to the heatmap. Note that hard calls are used for heatmap plotting. 
}
\value{
  See \code{heatmap}.
}
\author{ Wessel N. van Wieringen: \email{wvanwie@few.vu.nl} }
\seealso{ \code{heatmap}, \code{regioning}, \code{WECCAsc}, \code{WECCAhc}. }
\examples{
# generate object of class
data(WiltingCalled)

# make region data (soft and hard calls)
WiltingRegioned <- regioning(WiltingCalled)

# clustering with soft.calls
dendrogram <- WECCAsc(WiltingRegioned)

# generate a heatmap of the found clustering
WECCA.heatmap(WiltingRegioned, dendrogram)
}
