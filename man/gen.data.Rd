\name{gen.data}
\alias{gen.data}
\title{Simulate data sets}
\description{
The function generates one or more data matrices where the first couple
of variables are changed in class 1. The aim is to provide an easy means
to evaluate the performance of biomarker identification methods.}
\usage{
gen.data(nobj1, nobj2 = nobj1, nvar, n.biom = 2, group.diff = 0.2,
         nrep = 100, means = rep(0, nvar), cormat = diag(nvar))
}
\arguments{
  \item{nobj1, nobj2}{Numbers of objects in the two classes. If only
    nobj1 is given, the two classes are assumed to be of equal size.}
  \item{nvar}{Number of variables.}
  \item{n.biom}{Number of biomarkers, i.e. the number of variables to be
  changed in class 1 compared to class 2. The variables that are changed
  are always the first variables in the data matrix.}
  \item{group.diff}{group difference; the average difference between
    values of the biomarkers in the two classes.}
  \item{nrep}{Number of data sets to simulate.}
  \item{means}{Mean values of all variables, a vector.}
  \item{cormat}{Correlation matrix to be used in the simulation. Default
    is the identity matrix.}
}
\value{A list with the following elements:
  \item{X}{An array of dimension \code{nobj1 + nobj2} times \code{nvar}
    times \code{nrep}.}
  \item{Y}{The class vector.}
  \item{n.biomarkers}{The number of biomarkers.}
  \item{cormat}{The correlation matrix}
  \item{means1, means2}{The vectors of means for the two classes}
  }
\author{Ron Wehrens}
\examples{
X <- gen.data(10, nvar = 200)
names(X)
dim(X$X)
}
\keyword{manip}
\keyword{multivariate}
