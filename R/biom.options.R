### biom.options.R: modeled after pls.options.R

## The list of initial options, for the moment all pertaining to
## stability-based biomarkers selection

.biom.options <-
  list(max.seg = 100, oob.size = NULL, oob.fraction = .3,
       variable.fraction = .5, ntop = 10, min.present = .1,
       fmethods = c("studentt", "shrinkt", "pclda", "plsda", "vip"),
       univ.methods = c("studentt", "shrinkt"),
       nset = 10000, HCalpha = .1)

biom.options <- function(...) {
    if (nargs() == 0) return(.biom.options)
    current <- .biom.options
    temp <- list(...)
    if (length(temp) == 1 && is.null(names(temp))) {
        arg <- temp[[1]]
        switch(mode(arg),
               list = temp <- arg,
               character = return(.biom.options[arg]),
               stop("invalid argument: ", sQuote(arg)))
    }
    if (length(temp) == 0) return(current)
    n <- names(temp)
    if (is.null(n)) stop("options must be given by name")
    changed <- current[n]
    current[n] <- temp
    ## This assigns .biom.options in the global environment.  That way one
    ## can get back to the `factory defaults' by removing the variable from
    ## the global environment.  It also means that options are remembered
    ## between sessions (if the environment is saved).  Except for renaming
    ## .pls.Options to .biom.options, this is the only modification of the
    ## function:
    assign(".biom.options", current, pos = .GlobalEnv)
    invisible(current)
}
