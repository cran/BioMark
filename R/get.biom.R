get.biom <- function(X, Y,
                     fmethod = c("all", "pclda", "plsda",
                       "vip", "studentt"),
                     type = c("stability", "coef"),
                     ncomp = 2, max.seg = 100,
                     oob.size = NULL, oob.fraction = .1,
                     variable.fraction = 1,
                     ntop = 10, min.present = .1,
                     scale.p = "none", ...)
{
  type <- match.arg(type)

  univ.methods <- c("shrinkt", "studentt")
  
  fmethod <- match.arg(fmethod, c("all",
                                  "shrinkt", "studentt",
                                  "pclda", "plsda", "vip"),
                       several.ok = TRUE)
  if ("all" %in% fmethod)
    fmethod <- c("studentt", "pclda", "plsda",  "vip")

  multiv <- fmethod[!(fmethod  %in% univ.methods)]
  nmultiv <- length(multiv)
  univ <- fmethod[(fmethod  %in% univ.methods)]
  nuniv <- length(univ)
  fmethod <- c(univ, multiv) # do univariate methods first
  nncomp <- rep(c(1, length(ncomp)), c(nuniv, nmultiv))
  
  if (length(ncomp) > 1 & nmultiv > 0) {
    result <- vector(nuniv + nmultiv*length(ncomp), mode = "list")
    names(result) <- c(univ, paste(rep(multiv, each = length(ncomp)),
                                   " (",
                                   rep(ncomp, length(multiv)), ")",
                                   sep = ""))
  } else {
    result <- vector(length(fmethod), mode = "list")
    names(result) <- fmethod
  }
  
  if (type == "stability") {
    if (is.null(oob.size))
      oob.size <- round(.5 * oob.fraction * length(Y))
    ## .5 is necessary because there are two classes and oob.size is
    ## the number of oob samples in one class. For unequal class sizes
    ## this may be a problem.
    segments <- get.segments(Y, oob.size = oob.size, max.seg = max.seg)

    if (variable.fraction < 1) { # use different subsets of variables
      nvar <- round(variable.fraction * ncol(X))
      variables <- sapply(1:ncol(segments),
                          function(i) sample(ncol(X), nvar))
      nvars <- table(variables)
      if (length(nvars) < ncol(X))
        stop("Too few variables in resampling scheme:\ntry with a larger variable.fraction or use more segments.")
    } else {
      variables <- matrix(1:ncol(X), nrow = ncol(X), ncol = ncol(segments))
      nvars <- rep(max.seg, ncol(X))
    }
    
    fname <- paste(fmethod, "biom", sep = ".")
  } else {
    fname <- paste(fmethod, "coef", sep = ".")
    variables <- NULL
  }
  
  coef.order <- function(xx) {
    huhn.order <- order(abs(xx), decreasing = TRUE)
    list(biom.indices = huhn.order, coef.size = xx)
  }
  stab.order <- function(xx){
    huhn.order <- t(apply(abs(xx), 1, order, decreasing = TRUE))
    
    npicked <- table(huhn.order[,1:ntop])
    which.picked <- as.numeric(names(npicked))
    fraction.picked <- npicked / nvars[which.picked]
    
    selection <- sort(fraction.picked, decreasing = TRUE)
    selection <- selection[selection >= min.present]
    list(biom.indices = as.numeric(names(selection)),
         fraction.selected = selection)
  }

  counter <- 1
  for (m in seq(along = fmethod)) {
    ## for coef-based selection this is always a vector or a matrix
    ## for stability-based selection a vector or an array (if nncomp > 1)
    huhn.models <- do.call(fname[m], 
                           list(X = X, Y = Y,
                                segments = segments,
                                ncomp = ncomp,
##                                ntop = ntop,
##                                min.present = min.present,
                                scale.p = scale.p,
                                variables = variables, ...))
    
    if (type == "coef") {
      orderfun <- coef.order
      huhn.models <- array(huhn.models, c(1, ncol(X), nncomp[m]))
    } else {
      orderfun <- stab.order
      huhn.models <-
        array(huhn.models, c(nrow(huhn.models), ncol(huhn.models), nncomp[m]))
    }

    woppa <- lapply(1:(dim(huhn.models)[3]),
                    function(i, x) orderfun(x[,,i]),
                    huhn.models)
    for (mm in 1:nncomp[m]) {
      result[[counter]] <-  woppa[[mm]]
      counter <- counter + 1
    }
  }
  
  result
}

