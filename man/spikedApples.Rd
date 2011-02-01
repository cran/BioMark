\name{spikedApples}
\alias{spikedApples}
\docType{data}
\title{Metabolomics data on spiked apples}
\description{An data set of LC-MS features, obtained from twenty
  apples. The last ten apples are spiked with known compounds. This set
  provides a test case for biomarker selection methods: the task is to
  retrieve the true biomarker variables. The raw LC-MS data have been
  converted to CDF format and processed with XCMS to obtain the features.
}
\usage{data(spikedApples)}
\format{
  The format is a list of four elements:
  \describe{
    \item{mz}{the m/z values of the features (rounded)}
    \item{rt}{the retention times of the features}
    \item{dataMatrix}{the intensities of the features in the individual samples}
    \item{biom}{the indices of the "true" biomarkers}
  }
}
\author{Pietro Franceschi}
\references{
R. Wehrens, P. Franceschi, U. Vrhovsek and F. Mattivi. Stability-based
biomarker selection. Analytica Chimica Acta (2011), in press.
}
\examples{
data(spikedApples)
## show features identified in all apples
plot(spikedApples$rt, spikedApples$mz)
}
\keyword{datasets}
