\name{WECCAhc}
\alias{WECCAhc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Weighted clustering of array CGH data. }
\description{
  Weighted distance matrix calculation, followed by dendrogram construction for regioned hard called array CGH data. 
}
\usage{
WECCAhc(cghdata.regioned, dist.measure = "agree", linkage = "ward", weight.type = "all.equal")
}
\arguments{
  \item{cghdata.regioned}{ A \code{list}-object as returned by the function \code{regioning}. }
  \item{dist.measure}{ The distance measure to be used. This is either \code{"agree"} (the agreement similarity) or \code{"conc"} (the concordance similarity). }
  \item{linkage}{ The linkage method to be used, like \code{"single"}, \code{"average"}, \code{"complete"}, \code{"ward"}, \code{"total"}. }
  \item{weight.type}{ Region weighting to be used in the calculation of the distance, either \code{"all.equal"} or \code{"heterogeneity"}. }
}
\details{
  The distance between the call profiles of all sample pairs are calculated using special distance measures suitable for this data type.
  The distance between the call signature of two regions is either defined via the agreement or concordance similarity.
  The distance between two call profiles is defined as a weighted average of the distances between individual regions.
  Region weighing currently implemented: no-weighing (all regions contribute equally to the overall distance), or heterogeneity weighing (weights are proportional to their Shannon's entropy).
  Once the distance matrix has been calculated, this is submitted to the \code{hclust} function, which clusters the samples hierarchically, returning a dendrogram.
}
\value{
  An object of class \code{hclust} which describes the tree produced by the clustering process. See \code{hclust} for a list of its components.
}
\references{ 
Van Wieringen, W.N., Van de Wiel, M.A., Ylstra, B. (2008), "Weighted clustering of called aCGH data", \emph{Biostatistics}, 9(3), 484-500.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{dist}, \code{hclust}, \code{KLdiv}, \code{regioning} }
\examples{
# generate object of cghCall-class
data(WiltingCalled)

# make region data (soft and hard calls)
WiltingRegioned <- regioning(WiltingCalled)

# clustering with soft.calls
dendrogram <- WECCAhc(WiltingRegioned)
}
