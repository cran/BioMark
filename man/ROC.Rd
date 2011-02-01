\name{ROC}
\alias{ROC}
\alias{print.ROC}
\alias{plot.ROC}
\alias{lines.ROC}
\alias{points.ROC}
\alias{identify.ROC}
\alias{ROC.default}
\alias{roc.value}
\alias{AUC}
\title{ROC curves}
\description{Functions for making, plotting and analysing ROC curves.}
\usage{
ROC(TestResult, \dots)
\method{ROC}{default}(TestResult, D, take.abs = TRUE, \dots)
\method{plot}{ROC}(x, type = "b", null.line = TRUE,
xlab = "False Pos. Rate", ylab = "True Pos. Rate",
xlim = c(0, 1), ylim = c(0, 1), main = "ROC", \dots) 
\method{points}{ROC}(x, \dots)
\method{lines}{ROC}(x, \dots)
\method{identify}{ROC}(x, labels = NULL, \dots, digits = 1)
\method{print}{ROC}(x, \dots)
roc.value(found, true, totalN)
AUC(x, max.mspec)
}
\arguments{
  \item{TestResult}{Typically regression coefficients or t
    statistics. Note that when p values are used directly, the least
    significant values would be selected first. In this case one should
    use 1/p.}
  \item{D}{True, known, differences, either expressed as a vector of 0
    and 1 of the same length as \code{TestResult} or as a vector of
    indices.}
  \item{take.abs}{Logical, indicating whether to take absolute values of
    the test statistic.}
  \item{x}{An object of class ROC.}
  \item{type, xlab, ylab, xlim, ylim, main, labels, digits}{Standard
    arguments to functions like \code{plot} and \code{identify}.}
  \item{null.line}{Logical, whether to draw the line y = x,
    corresponding to random guessing.}
  \item{max.mspec}{Maximal value of the True Positive Rate to consider
    in AUC calculations. Setting this to a value smaller than one (which
    is the default) leads to a partial AUC value, which may in many
    cases be more useful.}
  \item{found}{The indices of the coefficients identified with a
    biomarker identification method.}
  \item{true}{The indices of the true biomarkers.}
  \item{totalN}{The total number of variables to choose from.}
  \item{\dots}{Further arguments, especially useful in the plotting
    functions.}
}
\value{
Function \code{ROC} returns a list with elements:
  \item{sens}{Sensitivity, or True Positive Rate (TPR).}
  \item{mspec}{1 - Specificity, or False Positive Rate (FPR).}
  \item{test}{levels of the test statistic.}
  \item{call}{Function call.}

Function \code{roc.value} returns a list with elements \code{sens} and
\code{mspec}, i.e., one point on a ROC curve.

Function \code{AUC} returns the area under the curve, measured up to the
value of \code{max.mspec} - if the latter is smaller than 1, it is a
partial AUC curve.
}

\references{
T. Lumley: ROC curves - in Programme's Niche, R News 4/1, June 2004.
}
\author{Ron Wehrens}
\examples{
set.seed(7)
huhn <- gen.data(10, nvar = 200, nrep = 20)

methods <- c("shrinkt", "studentt", "pclda", "plsda")
setnr <- 14 # check other sets too if you feel like it...
huhn.coef <- get.biom(huhn$X[,,setnr], huhn$Y,
                      fmethod = methods, type = "coef",
                      ncomp = 3, scale.p = "auto")
plot(ROC(abs(huhn.coef[[1]]$coef.size), 1:2), type = "l",
     main = paste("Set", setnr))
for (i in 2:4) 
  lines(ROC(abs(huhn.coef[[i]]$coef.size), 1:2), col = i, lty = i)
legend("bottomright", col = 1:4, lty = 1:4,
       legend = names(huhn.coef))

## Plot average ROC curves for simulated data
rock <- list(sens = matrix(0, dim(huhn$X)[[3]], dim(huhn$X)[[2]]),
             mspec = matrix(0, dim(huhn$X)[[3]], dim(huhn$X)[[2]]))
rocs <- rep(list(rock), length(methods))

for (i in 1:dim(huhn$X)[3]) {
  huhn.coef <- get.biom(huhn$X[,,i], huhn$Y,
                        fmethod = methods, type = "coef",
                        ncomp = 3, scale.p = "auto")
  for (j in seq(along = methods)) {
    woppa <- ROC(abs(huhn.coef[[j]]$coef.size), 1:huhn$n.biomarkers)
    rocs[[j]][[1]][i,] = c(woppa$sens)
    rocs[[j]][[2]][i,] = c(woppa$mspec)
  }
}

par(mfrow = c(2,2))
for (j in seq(along = methods)) {
  matplot(t(rocs[[j]][[2]]), t(rocs[[j]][[1]]), type = "l", col = "gray",
          xlim = c(0, 1), ylim = c(0, 1), main = methods[j],
          xlab = "TPR", ylab = "FPR")
  lines(colMeans(rocs[[j]][[2]]), colMeans(rocs[[j]][[1]]), col = 2)
}

## Show distributions of AUC values for a simulated scenario
all.aucs <- matrix(0, dim(huhn$X)[3], length(methods))
for (i in 1:dim(huhn$X)[3]) {
  huhn.coef <- get.biom(huhn$X[,,i], huhn$Y,
                        fmethod = methods, type = "coef",
                        ncomp = 3, scale.p = "auto")
  all.aucs[i,] <- sapply(huhn.coef,
                         function(x) AUC(ROC(abs(x$coef.size),
                                             1:huhn$n.biomarkers)))
}
boxplot(all.aucs ~ col(all.aucs), notch = TRUE, main = "AUC",
        col = 1:length(methods), names = methods)
}
\keyword{classif}
