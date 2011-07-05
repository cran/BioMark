pclda.stab <- function(X, Y, ncomp = 2, scale.p = NULL,
                       segments = NULL, variables = NULL, ...)
{
  if (is.null(segments)) 
    segments <- get.segments(Y)
  if (is.null(variables)) ## default: use all variables
    variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))

  x.coef <- array(NA, c(ncol(segments), ncol(X), length(ncomp)))
  for (i in 1:ncol(segments))
    x.coef[i,variables[,i],] <-
      pclda.coef(X[-segments[,i], variables[,i]], Y[-segments[,i]],
                 ncomp = ncomp, scale.p = scale.p, ...)

  x.coef
}

plsda.stab <- function(X, Y, ncomp = 2, scale.p = NULL,
                       segments = NULL, variables = NULL, ...)
{
  if (is.null(segments))
    segments <- get.segments(Y)
  if (is.null(variables))
    variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))
  
  x.coef <- array(NA, c(ncol(segments), ncol(X), length(ncomp)))
  for (i in 1:ncol(segments))
    x.coef[i,variables[,i],] <-
      plsda.coef(X[-segments[,i], variables[,i]], Y[-segments[,i]],
                 ncomp = ncomp, scale.p = scale.p, ...)
  
  x.coef
}

vip.stab <- function(X, Y, ncomp = 2, scale.p = NULL,
                     segments = NULL, variables = NULL, ...)
{
  if (is.null(segments))
    segments <- get.segments(Y)
  if (is.null(variables))
    variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))
  
  x.coef <- array(NA, c(ncol(segments), ncol(X), length(ncomp)))
  for (i in 1:ncol(segments)) 
    x.coef[i,variables[,i],] <-
      vip.coef(X[-segments[,i], variables[,i]], Y[-segments[,i]],
               ncomp = ncomp, scale.p = scale.p, ...)
  
  x.coef
}

### the dots in the shrinkt.stab and studentt.stab functions are
### necessary to catch extra arguments to other functions.
shrinkt.stab <- function(X, Y, scale.p = NULL,
                         segments = NULL, variables = NULL, ...)
{
  if (is.null(segments))
    segments <- get.segments(Y)
  if (is.null(variables))
    variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))
  
  x.coef <- matrix(NA, ncol(segments), ncol(X))
  for (i in 1:ncol(segments))
    x.coef[i,variables[,i]] <- shrinkt.coef(X[-segments[,i], variables[,i]], 
                                            Y[-segments[,i]],
                                            scale.p = scale.p, ...)
  
  array(x.coef, c(ncol(segments), ncol(X), 1))
}

studentt.stab <- function(X, Y, scale.p = NULL,
                          segments = NULL, variables = NULL, ...)
{
  if (is.null(segments))
    segments <- get.segments(Y)
  if (is.null(variables))
    variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))
  
  x.coef <- matrix(NA, ncol(segments), ncol(X))
  for (i in 1:ncol(segments))
    x.coef[i,variables[,i]] <- studentt.coef(X[-segments[,i], variables[,i]], 
                                             Y[-segments[,i]],
                                             scale.p = scale.p, ...)
  
  array(x.coef, c(ncol(segments), ncol(X), 1))
}

