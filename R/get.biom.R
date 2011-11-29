get.biom <- function(X, Y, fmethod = "all",
                     type = c("stab", "HC", "coef"),
                     ncomp = 2, biom.opt = biom.options(),
                     scale.p = "auto", ...)
{
  fmethod <- match.arg(fmethod, c("all", biom.opt$fmethods),
                       several.ok = TRUE)
  if ("all" %in% fmethod)
    fmethod <- biom.opt$fmethods

  multiv <- fmethod[!(fmethod  %in% biom.opt$univ.methods)]
  nmultiv <- length(multiv)
  univ <- fmethod[(fmethod  %in% biom.opt$univ.methods)]
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
  
  type <- match.arg(type)
  mod.method <- ifelse(type == "stab", "stab", "coef")
  fname <- paste(fmethod, mod.method, sep = ".")

  ## Get settings from the biom.opt argument, mostly for stability-based BS
  if (type == "stab") {
    oob.size <- biom.opt$oob.size
    oob.fraction <- biom.opt$oob.fraction
    smallest.class.fraction <- min(table(Y) / length(Y))
    ## for equal class sizes this is .5
    if (is.null(oob.size))
      oob.size <- round(smallest.class.fraction * oob.fraction * length(Y))
    max.seg <- biom.opt$max.seg
    segments <- get.segments(Y, oob.size = oob.size, max.seg = max.seg)
    
    variable.fraction <- biom.opt$variable.fraction
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
  } else {
    variables <- NULL
    
    if (type == "HC") {
      nset <- biom.opt$nset
      HCalpha <- biom.opt$HCalpha
    }
  }

  ## Compared to earlier versions: treat HC separately because of the
  ## expensive evaluation of null distributions
  ## Temporary solution - not pretty though - is to do HC in the same
  ## way only for the univariate methods and to treat the multivariate
  ## methods separately. Take care that with future versions this may
  ## have to be revised.
  counter <- 1
  for (m in seq(along = fmethod)) {
    ## Here the real work is done: call the modelling functions
    huhn.models <- do.call(fname[m], 
                           list(X = X, Y = Y,
                                segments = segments,
                                ncomp = ncomp,
                                scale.p = scale.p,
                                variables = variables, ...))
    
    ## for coef-based selection this is always a vector or a matrix,
    ## for stability-based selection a vector or an array (if nncomp > 1)
    switch(type,
           coef = {
             orderfun <- function(xx) {
               huhn.order <- order(abs(xx), decreasing = TRUE)
               list(biom.indices = huhn.order, coef.size = xx)
             }
             huhn.models <- array(huhn.models, c(1, ncol(X), nncomp[m]))
             woppa <- lapply(1:dim(huhn.models)[3],
                             function(i, x) orderfun(x[,,i]),
                             huhn.models)
           },
           stab = {
             orderfun <- function(xx){
               huhn.order <- t(apply(abs(xx), 1, order, decreasing = TRUE))

               npicked <- 1:ncol(huhn.order)
               names(npicked) <- 1:ncol(huhn.order)
               npicked <- table(huhn.order[,1:biom.opt$ntop])
               which.picked <- as.numeric(names(npicked))
               fraction.picked <- npicked / nvars[which.picked]
               
               selection <- sort(fraction.picked, decreasing = TRUE)
               selection <- selection[selection >= biom.opt$min.present]
               list(biom.indices = as.numeric(names(selection)),
                    fraction.selected = selection)
             }
             
             huhn.models <-
               array(huhn.models, c(nrow(huhn.models), ncol(huhn.models),
                                    nncomp[m]))
             woppa <- lapply(1:dim(huhn.models)[3],
                             function(i, x) orderfun(x[,,i]),
                             huhn.models)
           },
           HC = {
             if (m <= nuniv) {
               if (is.vector(huhn.models))
                 huhn.models <- matrix(huhn.models, ncol = 1)
               huhn.pvals <-
                 apply(huhn.models, 2,
                       function(x) 2*(1 - pt(abs(x), nrow(X) - 2)))
               woppa <-
                 lapply(1:ncol(huhn.models),
                        function(i)
                        list(biom.indices =
                             HCthresh(huhn.pvals[,i], alpha = HCalpha,
                                      plotit = FALSE),
                             coef.size = huhn.models[,i]))
             } else { # just return something, real calcs later
               woppa <- lapply(1:ncol(huhn.models),
                               function(i)
                               list(biom.indices = NULL,
                                    coef.size = huhn.models[,i]))
             }
           })
    
    for (mm in 1:nncomp[m]) {
      result[[counter]] <-  woppa[[mm]]
      counter <- counter + 1
    }
  }

  if (type == "HC" & nmultiv > 0) {
    ## Possible PCLDA, PLSDA and VIP calculations for HC are done here
    which.pclda <- which(substr(names(result), 1, 5) == "pclda")
    if (length(which.pclda) > 0) {
      ## pval.pclda returns pvalues for each number of components, i.e.,
      ## a matrix, possibly with one column
      huhn.models <- pval.pclda(X, Y, ncomp, scale.p, nset)
      for (mm in seq(along = ncomp)) {
        result[[ which.pclda[mm] ]]$biom.indices <-
          HCthresh(huhn.models[,mm], alpha = HCalpha)
      }
    }
    which.plsda <- which(substr(names(result), 1, 5) == "plsda")
    which.vip <- which(substr(names(result), 1, 3) == "vip")
    
    if (length(which.plsda) > 0 | length(which.vip > 0)) {
      ## pval.pclda returns pvalues for each number of components, i.e.,
      ## a matrix, possibly with one column
      if (length(which.plsda) > 0) {
        if (length(which.vip > 0)) {
          smethod <- "both"
        } else {
          smethod <- "plsda"
        }
      } else {
        smethod <- "vip"
      }
      huhn.models <- pval.plsdavip(X, Y, ncomp, scale.p, nset, smethod)
      if (length(which.plsda) > 0) {
        for (mm in seq(along = ncomp)) {
          result[[ which.plsda[mm] ]]$biom.indices <-
            HCthresh(huhn.models[,mm,"plsda"], alpha = HCalpha)
        }
      }
      if (length(which.vip) > 0) {
        for (mm in seq(along = ncomp)) {
          result[[ which.vip[mm] ]]$biom.indices <-
            HCthresh(huhn.models[,mm, "vip"], alpha = HCalpha)
        }
      }
    }
  }
  
  result2 <- c(result, list(info = list(call = match.call(),
                              type = type, fmethod = fmethod,
                              nvar = ncol(X),
                              biom.opt = biom.opt)))
  class(result2) <- "BMark"

  result2
}


coef.BMark <- function(object, ...) {
  if (object$info$type != "coef")
    stop("Coefficients only available from coefficient-based biomarker selection")

  sapply(object[names(object) != "info"], function(xx) xx$coef.size)
}

print.BMark <- function(x, ...) {
  type <- x$info$type
  switch(type,
         coef = cat("Result of coefficient-based biomarker selection using ",
           length(x)-1, " modelling method",
           ifelse(length(x) > 2, "s", ""), ".\n", sep = ""),
         HC = cat("Result of HC-based biomarker selection using ",
           length(x)-1, " modelling method",
           ifelse(length(x) > 2, "s", ""), ".\n", sep = ""),
         cat("Result of stability-based biomarker selection using ",
             length(x)-1, " modelling method",
             ifelse(length(x) > 2, "s", ""), ".\n", sep = ""))
}

summary.BMark <- function(object, ...) {
  type <- object$info$type
  nslots <- length(object)
  infoslot <- which(names(object) == "info")
  switch(type,
         coef = {
           cat("Result of coefficient-based biomarker selection using ",
               nslots-1, " modelling method",
           ifelse(length(object) > 2, "s", ""), ":\n", sep = "")
           cat(names(object)[-infoslot], "\n")
           cat("Total number of variables in X matrix:", 
               object[infoslot]$nvar, "\n")
         },
         {
           typestr <- ifelse(type == "HC", "HC-based", "stability-based")
           cat("Result of ", typestr, " biomarker selection using ",
               nslots-1, " modelling method",
           ifelse(length(object) > 2, "s", ""), ":\n", sep = "")
           cat(names(object)[-infoslot], "\n")
           cat("Total number of variables:",
               object[infoslot]$nvar, "\n")
           cat("Number of variables selected:\n")
           nsel <- sapply(object[-infoslot],
                          function(xx) length(xx$biom.indices))
           names(nsel) <- names(object)[-infoslot]
           print(nsel)
         }
         )
}
