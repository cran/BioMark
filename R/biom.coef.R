pclda.coef <- function(X, Y, ncomp = 2, scale.p = NULL, ...)
{
  if (is.factor(Y)) Y <- as.numeric(Y)
  if (length(table(Y)) != 2) stop("only two-class discrimination implemented")
  
  FUN <- scalefun(scale.p)
  matrix(svdpc.fit(FUN(X), Y, ncomp = max(ncomp),
                   stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}


## Changed to widekernelpls.fit because this probably is the most
## relevant situation  
plsda.coef <- function(X, Y, ncomp = 2, scale.p = NULL, ...)
{
  if (is.factor(Y)) Y <- as.numeric(Y)
  if (length(table(Y)) != 2) stop("only two-class discrimination implemented")
  
  FUN <- scalefun(scale.p)
  matrix(widekernelpls.fit(FUN(X), Y, ncomp = max(ncomp),
                           stripped = TRUE)$coefficients[, 1, ncomp],
         ncol(X), length(ncomp))
}

vip.coef <- function(X, Y, ncomp = 2, scale.p = NULL, ...)
{ ## careful with the next line! VIPs change depending on the scale of
  ## Y, so it matters whether you use 0/1, -1/1, or other codings.
  if (is.factor(Y)) Y <- as.numeric(Y)
  FUN <- scalefun(scale.p)

  plsmod <- plsr(Y ~ FUN(X), ncomp = max(ncomp), method = "widekernelpls")
  ww <- loading.weights(plsmod)

  result <- matrix(NA, ncol(X), length(ncomp))
  for (i in 1:length(ncomp)) {
    var.exp <- diff(c(0, R2(plsmod, estimate = "train",
                            ncomp = 1:ncomp[i], intercept = FALSE)$val))

    result[,i] <- sqrt(ncol(X) * ww[,1:ncomp[i],drop = FALSE]^2 %*%
                       var.exp / sum(var.exp))
  }

  result
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

