\name{regioning}
\alias{regioning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Regioning of called array CGH data }
\description{
Region construction in DNA copy number data using hard calls.
Given the resulting regions, the mediod call probability (soft calls) signature
of each region is calculated.}
\usage{
  regioning(cghdata.called, threshold = 1e-05)
}
\arguments{
  \item{cghdata.called}{ Object of class \code{cghCall}, containing annotion, calls and their probabilities. }
  \item{threshold}{ Threshold for the function \code{CGHregions}. It determines when a region is enforced. }
}
\details{
  CGHregions is applied to the calls of an object of class \code{cghCall} (as returned by CGHcall version 1.2.0 or higher). Using the regions demarcations the 
mediod call probability signatures are calculated. The matrix of these mediod signatures are considered 
the regioned soft call data (see warning below).
}
\value{
  A \code{list}-object with the following entries:
  \item{ann }{\code{data.frame} containing the annotation data of the regions. The first column contains the region's chromosome information, the second and third its start and end base pair, and the fourth column the number of features constituting the region.}
  \item{hardcalls }{The regioned hard calls \code{matrix}.}
  \item{softcalls }{The regioned soft calls \code{matrix}.}
}
\references{ 
Van de Wiel, M.A., Kim, K.I., Vosse, S.J., Van Wieringen, W.N., Wilting, S.M. , Ylstra, B. (2007), "CGHcall: an algorithm for calling aberrations for multiple array CGH tumor profiles", \emph{Bioinformatics}, 23, 892-894. 

Van de Wiel, M.A., Van Wieringen, W.N. (2007), "CGHregions: dimension reduction for array CGH data with minimal information loss", \emph{Cancer Informatics}, 2, 55-63.

Van Wieringen, W.N., Van de Wiel, M.A., Ylstra, B. (2007), "Normalized, segmented or called aCGH data?", \emph{Cancer Informatics}, 3, 331-337
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\section{Warning }{
  The current version of CGHregions is designed for hard calls only. This is likely to change in the near future, making this function obsolete.
}
\seealso{ \code{CGHregions}, \code{CGHcall} }
\examples{
# generate object of cghCall class object
data(WiltingCalled)

# make region data (soft and hard calls)
WiltingRegioned <- regioning(WiltingCalled)
}
