\name{get.biom}
\alias{get.biom}
\alias{pclda.coef}
\alias{pclda.biom}
\alias{plsda.coef}
\alias{vip.coef}
\alias{vip.biom}
\alias{plsda.biom}
\alias{plsda.biom}
\alias{shrinkt.biom}
\alias{shrinkt.coef}
\alias{studentt.fun}
\alias{studentt.biom}
\alias{studentt.coef}
\title{Get biomarkers discriminating between two classes}
\description{Biomarkers can be identified in two ways: the classical way
  is to look at those variables with large model coefficients. The other
  way is to assess the stability of these coefficients under subsampling
  of the data set. In this implementation, we use crossvalidation using
  every possible combination of single samples from the two classes as
  out-of-bag set.
}
\usage{
pclda.coef(X, Y, ncomp = 2, scale.p = NULL, \dots)
pclda.biom(X, Y, ncomp = 2, scale.p = NULL,
           segments = NULL, variables = NULL, \dots)
%
plsda.coef(X, Y, ncomp = 2, scale.p = NULL, \dots)
plsda.biom(X, Y, ncomp = 2, scale.p = NULL,
           segments = NULL, variables = NULL, \dots)
%
vip.coef(X, Y, ncomp = 2, scale.p = NULL, \dots)
vip.biom(X, Y, ncomp = 2, scale.p = NULL,
         segments = NULL, variables = NULL, \dots)
%
shrinkt.coef(X, Y, scale.p = NULL, \dots)
shrinkt.biom(X, Y, scale.p = NULL,
           segments = NULL, variables = NULL, \dots)
%
studentt.fun(Y)
studentt.coef(X, Y, scale.p = NULL, \dots)
studentt.biom(X, Y, scale.p = NULL,
           segments = NULL, variables = NULL, \dots)
%
get.biom(X, Y, fmethod, type, ncomp, max.seg, oob.size, oob.fraction,
         variable.fraction, ntop, min.present, scale.p, \dots)
}
\arguments{
  \item{X}{Data matrix. Usually the number of columns (variables) is
    (much) larger than the number of rows (samples).}
  \item{Y}{Class indication. Either a factor, or a numeric vector.}
  \item{ncomp}{Number of latent variables to use in PCLDA and PLSDA
    modelling. In function \code{get.biom} this may be a vector; in all
    other functions it should be one number.}
  \item{scale.p}{Scaling. This is performed individually in every
    crossvalidation iteration, and can have a profound effect on the
    results.}
  \item{segments}{In stability-based biomarker identification, the
    segments constituting the out-of-box samples in each crossvalidation
    iteration. If not specified, all possible combinations between
    individual samples from both classes will be used.}
  \item{variables}{indices of variables to include, or NULL to include all.}
  \item{fmethod}{Modelling method(s) employed. The default is to use
    \code{"all"}, which will test \code{"pclda"}, \code{"plsda"},
    \code{"vip"}, and \code{"studentt"}. Also available is
    \code{"shrinktt"}.}
  \item{type}{Whether to use coefficient size as a criterion
    (\code{"coef"}), or \code{"stability"}.}
  \item{ntop}{The number of "top" coefficients taken into account in
    stability-based biomarker identification. If a variable appears
    consistently among the \code{ntop} biggest coefficients, it is said
    tobe stable.}
  \item{min.present}{The minimal fraction of times a variable should be in
    the top list to be considered as a potential biomarker
    (stability-based identification). Setting this
    argument to 0 will lead to a list containing all coefficients that
    were present in the top list at least once - a value of 1 only
    returns those variables that are selected in every
    iteration. Default: .1.}
  \item{oob.size, oob.fraction}{Size of the out-of-bag fraction, either
    given as an absolute number (oob.size) or as a fraction. Default is
    to leave out ten percent. If oob.size is given explicitly, it takes
    precedence over oob.fraction.}
  \item{max.seg}{Maximal number of jackknife iterations.}
  \item{variable.fraction}{Use 1 to always include all variables - use a
    smaller fraction to have a different random subset of all variables
    in each iteration (stability-based identification).}
  \item{\dots}{Further arguments for modelling functions. Often used to
    catch unused arguments.}
}
\value{When \code{type == "coef"}, the function returns a list
  with elements \code{biom.indices} and \code{coef.size}: these contain
  the values of the coefficients, and the indices of the variables
  sorted on absolute value of the coefficients. Note that the
  coefficients themselves are in the original order. When \code{type ==
    "stability"}, the second list element is 
  \code{times.selected} rather than \code{coef.size}: this indicates in
  what fraction of the crossvalidation iterations a particular variable
  has been selected. The more times it has been selected, the more
  stable it is as a biomarker.}
\author{Ron Wehrens}
\seealso{\code{\link{get.segments}}, \code{\link{scalefun}}}
\examples{
## Real apple data
data(spikedApples)
apple.coef <- get.biom(X = spikedApples$dataMatrix,
                       Y = rep(1:2, each = 10))
TPs <- sapply(apple.coef, function(x, y) sum(x$biom.indices \%in\% y),
              spikedApples$biom) ## all methods find 4/5
FPs <- sapply(apple.coef, function(x, y) sum(!(x$biom.indices \%in\% y)),
              spikedApples$biom) ## VIP finds one fewer FP

## Simulated data
set.seed(7)
huhn <- gen.data(10, nvar = 200, nrep = 1)

huhn.coef <- get.biom(huhn$X[,,1], huhn$Y, fmethod = "all", type = "coef",
                      ncomp = 3, scale.p = "auto")

huhn.pclda <- get.biom(huhn$X[,,1], huhn$Y, fmethod = "pclda",
                       type = "stability",
                       ncomp = 2:4, scale.p = "pareto")
}
\keyword{multivariate}
\keyword{models}
