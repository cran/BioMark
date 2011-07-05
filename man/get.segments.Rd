\name{get.segments}
\alias{get.segments}
\title{Crossvalidation segments}
\description{Provides all possible combinations of samples from the two
  classes to be left out in crossvalidation, with a maximum given by
  parameter \code{max.seg}.}
\usage{
get.segments(i1, i2 = NULL, oob.size = 1, max.seg = 100)
}
\arguments{
  \item{i1}{either an index vector for objects in class 1, or a
    classification vector (factor, or numeric), from which the indices
    of both classes can be derived.}
  \item{i2}{if non-NULL, vector indexing objects in class 2.}
  \item{oob.size}{number of samples to be left out in every
    iteration. If one (the default), this corresponds to LOO
    crossvalidation.}
  \item{max.seg}{maximal number of segments to return. If null, all
    possible combinations are returned -- this option is only possible
    if oob.size equals 1. If oob.size is larger, max.seg must be defined
    since the number of possibilities becomes too large.}
}
\value{
  Returns a matrix where the columns contain the numbers of the samples to
  be left out in the respective crossvalidation iterations for both
  classes.
}
\author{
  Ron Wehrens
}
\seealso{\code{\link{get.biom}}}
\examples{
i1 <- seq(1, 10, by = 2)
i2 <- seq(2, 15, by = 2)
get.segments(i1, i2)
get.segments(i1, i2, max.seg = 10)
get.segments(i1, i2, oob.size = 2, max.seg = 10)

I <- rep(1:2, c(5,6))
get.segments(I)
get.segments(I, max.seg = 15)
}
\keyword{manip}
\keyword{iteration}
