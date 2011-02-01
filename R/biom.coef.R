pclda.coef <- function(X, Y, ncomp = 2, scale.p = NULL, ...)
{
  if (is.factor(Y)) Y <- as.numeric(Y)
  if (length(table(Y)) != 2) stop("only two-class discrimination implemented")
  
  FUN <- scalefun(scale.p)
  matrix(svdpc.fit(FUN(X), Y, ncomp = max(ncomp),
                   stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}

plsda.coef <- function(X, Y, ncomp = 2, scale.p = NULL, ...)
{
  if (is.factor(Y)) Y <- as.numeric(Y)
  if (length(table(Y)) != 2) stop("only two-class discrimination implemented")
  
  FUN <- scalefun(scale.p)
  matrix(kernelpls.fit(FUN(X), Y, ncomp = max(ncomp),
                       stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}

vip.coef <- function(X, Y, ncomp = 2, scale.p = NULL, ...)
{
  FUN <- scalefun(scale.p)

  plsmod <- plsr(Y ~ FUN(X), ncomp = max(ncomp))
  ww <- loading.weights(plsmod)

  result <- matrix(NA, ncol(X), length(ncomp))
  for (i in 1:length(ncomp)) {
    var.exp <- diff(c(0, R2(plsmod, estimate = "train",
                            ncomp = 1:ncomp[i], intercept = FALSE)$val))

    result[,i] <- sqrt(ncol(X) * ww[,1:ncomp[i]]^2 %*% var.exp / sum(var.exp))
  }

  result
}

## direct copy from the st package, just to get rid of all the output...
studentt.fun <- function (Y)
{
  function(X) {
    tmp <- centroids(X, Y, var.pooled = TRUE, var.groups = FALSE,
                     shrink = FALSE, verbose = FALSE)
    diff <- tmp$means[, 1] - tmp$means[, 2]
    n1 <- tmp$samples[1]
    n2 <- tmp$samples[2]
    v <- tmp$var.pooled
    sd <- sqrt((1/n1 + 1/n2) * v)
    
    diff/sd
  }
}

studentt.coef <- function(X, Y, scale.p = NULL, ...)
{
  if (is.factor(Y)) Y <- as.numeric(Y)
  if (length(table(Y)) != 2) stop("only two-class discrimination implemented")
  
  FUN <- scalefun(scale.p)
  TFUN <- studentt.fun(Y)
  
  TFUN(FUN(X))
}

shrinkt.coef <- function(X, Y, scale.p = NULL, ...)
{
  if (is.factor(Y)) Y <- as.numeric(Y)
  if (length(table(Y)) != 2) stop("only two-class discrimination implemented")
  
  FUN <- scalefun(scale.p)
  TFUN <- shrinkt.fun(L =  Y, FALSE, FALSE)
  
  TFUN(FUN(X))
}

