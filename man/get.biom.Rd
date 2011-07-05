\name{get.biom}
\alias{get.biom}
\alias{print.BMark}
\alias{summary.BMark}
\alias{coef.BMark}
\title{Get biomarkers discriminating between two classes}
\description{Biomarkers can be identified in several ways: the classical way
  is to look at those variables with large model coefficients or large t
  statistics. One other is based on the higher criticism approach (HC),
  and the third possibility assesses the stability of these coefficients
  under subsampling of the data set.
}
\usage{
get.biom(X, Y, fmethod = "all", type = c("stab", "HC", "coef"),
         ncomp = 2, biom.opt = biom.options(), scale.p = "auto",
         \dots)
\method{coef}{BMark}(object, \dots)
\method{print}{BMark}(x, \dots)
\method{summary}{BMark}(object, \dots)
}
\arguments{
  \item{X}{Data matrix. Usually the number of columns (variables) is
    (much) larger than the number of rows (samples).}
  \item{Y}{Class indication. Either a factor, or a numeric vector.}
  \item{fmethod}{Modelling method(s) employed. The default is to use
    \code{"all"}, which will test all methods in the current
    \code{biom.options$fmethods} list.}
  \item{type}{Whether to use coefficient size as a criterion
    (\code{"coef"}), \code{"stab"} or \code{"HC"}.}
  \item{ncomp}{Number of latent variables to use in PCLDA and PLSDA
    modelling. In function \code{get.biom} this may be a vector; in all
    other functions it should be one number. Default: 2.}
  \item{biom.opt}{Options for the biomarker selection - a list with
    several named elements. See \code{\link{biom.options}}.}
  \item{scale.p}{Scaling. This is performed individually in every
    crossvalidation iteration, and can have a profound effect on the
    results. Default: "none". Other possible choices: "auto" for
    autoscaling, "pareto" for pareto scaling, "log" and "sqrt" for log
    and square root scaling, respectively.}
  \item{object, x}{A BMark object.}
  \item{\dots}{Further arguments for modelling functions. Often used to
    catch unused arguments.}
}
\value{Returns an object of class "BMark", a list always containing an
  element \code{info}. When \code{type == "coef"} and
  \code{type == "HC"}, the object also contains elements
  \code{biom.indices} and \code{coef.size}: these contain the indices of
  the variables sorted on absolute value of the coefficients, and the
  values of the coefficients, respectively. Note that the 
  coefficients themselves are in the original order.
  
  When \code{type == "stab"}, \code{times.selected} is returned
  rather than \code{coef.size}: this indicates in what fraction of the
  crossvalidation iterations a particular variable has been
  selected. The more times it has been selected, the more stable it is
  as a biomarker.}
\author{Ron Wehrens}
\seealso{\code{\link{biom.options}}, \code{\link{get.segments}},
  \code{\link{scalefun}}}
\examples{
## Real apple data (small set)
data(spikedApples)
apple.coef <- get.biom(X = spikedApples$dataMatrix,
                       Y = rep(1:2, each = 10),
                       ncomp = 2:3, type = "coef")

## stability-based selection
set.seed(7)
apple.stab <- get.biom(X = spikedApples$dataMatrix,
                       Y = rep(1:2, each = 10),
                       ncomp = 2:3, type = "stab")
sapply(apple.stab[names(apple.stab) != "info"],
       function(x, y) sum(x$biom.indices \%in\% y),
       spikedApples$biom)
## TPs (stab): all find 5/5
sapply(apple.stab[names(apple.stab) != "info"],
       function(x, y) sum(!(x$biom.indices \%in\% y)),
       spikedApples$biom)
## FPs (stab): all methods find many FPs, PCLDA even more
## than the others

\dontrun{
## Same for HC-based selection
## Warning: takes a long time!
apple.HC <- get.biom(X = spikedApples$dataMatrix,
                     Y = rep(1:2, each = 10),
                     ncomp = 2:3, type = "HC")
sapply(apple.HC[names(apple.HC) != "info"],
       function(x, y) sum(x$biom.indices \%in\% y),
       spikedApples$biom)
sapply(apple.HC[names(apple.HC) != "info"],
       function(x, y) sum(!(x$biom.indices \%in\% y)),
       spikedApples$biom)
}
}
\keyword{multivariate}
\keyword{models}
\keyword{htest}
