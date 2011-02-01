require(MASS)
gen.data <- function(nobj1, nobj2 = nobj1, nvar, n.biom = 2,
                     group.diff = .2, nrep = 100,
                     means = rep(0, nvar), cormat = diag(nvar))
{
  nobj <- nobj1 + nobj2

  X <- array(0, c(nobj, nvar, nrep))
  diffvec <- .5*rep(c(group.diff, 0), c(n.biom, nvar-n.biom))
  means1 <- means + diffvec
  means2 <- means - diffvec

  for (i in 1:nrep)
    X[,,i] <- rbind(mvrnorm(nobj1, means1, cormat),
                    mvrnorm(nobj2, means2, cormat))
  
  list(X = X, Y = rep(c(0, 1), each = c(nobj1, nobj2)),
       n.biomarkers = n.biom,
       means1 = means1, means2 = means2, cormat = cormat)
}

