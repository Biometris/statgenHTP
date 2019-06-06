#' Fit spatial models per time point
#'
#' Fit spatial models per time point in an object of class TP.
#'
#' @param TP An object of class TP.
#' @param trait A character string indicating the trait used as response
#' variable in the model.
#' @param covariates A character vector indicating the variables used as
#' fixed effects in the model.
#' @param geno.decomp A character vector indicating the variables used to
#' group the genotypes in the model.
#' @param useCheck Should check genotypes be used as an extra factor in the
#' model?
#'
#' @return An object of class fitMod, a list of fitted spatial models.
#'
#' @export
fitModels <- function(TP,
                      trait,
                      covariates = NULL,
                      geno.decomp = NULL,
                      useCheck = FALSE,
                      engine = c("SpATS", "asreml")) {
  ## Checks.
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  if (!all(sapply(X = TP, FUN = hasName, name = trait))) {
    stop(trait, " should be a column in TP for all timePoints.\n")
  }
  if (!is.null(covariates)) {
    for (covar in covariates) {
      if (!all(sapply(X = TP, FUN = hasName, name = covar))) {
        stop(covar, " should be a column in TP for all timePoints.\n")
      }
    }
  }
  if (!is.null(geno.decomp)) {
    for (gd in geno.decomp) {
      if (!all(sapply(X = TP, FUN = hasName, name = gd))) {
        stop(gd, " should be a column for all timePoints.\n")
      }
    }
  }
  if (useCheck) {
    if (!all(sapply(X = TP, FUN = hasName, name = "check"))) {
      stop("check should be a column in TP for all timePoints.\n")
    }
  }
  engine <- match.arg(engine)
  ## If geno.decomp is used genotype and covariates have to be replaced by
  ## an interaction of genotype and covariates with the geno.decomp variables.
  ## Construct an interaction of all variables in geno.decomp.
  if (length(geno.decomp) > 1) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      timePoint[["geno.decomp"]] <- interaction(timePoint[geno.decomp],
                                                sep = "_")
      return(timePoint)
    })
    ## Set geno.decomp to newly constructed variable.
    geno.decomp <- "geno.decomp"
  }
  ## Replace genotype and covariates by their interaction with geno.decomp.
  if (!is.null(geno.decomp)) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      timePoint[["genotype"]] <- interaction(timePoint[[geno.decomp]],
                                             timePoint[["genotype"]], sep = "_")
      for (covar in covariates) {
        timePoint[[covar]] <- interaction(timePoint[[geno.decomp]],
                                          timePoint[[covar]], sep = "_")
      }
      return(timePoint)
    })
  }
  ## All covariates should be factors. If not convert them to factor.
  for (covar in covariates) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      if (!is.factor(timePoint[[covar]])) {
        timePoint[[covar]] <- as.factor(timePoint[[covar]])
      }
      return(timePoint)
    })
  }
  ## Fixed part consists of covariates, geno.decomp and check.
  if (!is.null(c(covariates, geno.decomp))) {
    fixedForm <- formula(paste("~", paste(c(covariates, geno.decomp),
                                          collapse = "+"),
                               if (useCheck) "+ check"))
  } else {
    fixedForm <- if (useCheck) formula("~ check") else NULL
  }
  ## Get column containing genotype.
  genoCol <- if (useCheck) "genoCheck" else "genotype"
  if (engine == "SpATS") {
    ## Loop on timepoint to run SpATS.
    fitMods <- lapply(X = TP, function(timePoint) {
      message(timePoint[["timePoint"]][1])
      ## Only keep columns needed for analysis.
      modCols <- c("timePoint", "plotId", "genotype", "genoCheck", "check",
                   "colId", "rowId", "colNum", "rowNum", covariates,
                   geno.decomp, trait)
      modDat <- timePoint[colnames(timePoint) %in% modCols]
      modDat <- droplevels(modDat)
      ## number of segments for SpATS.
      nseg = c(nlevels(modDat[["colId"]]), nlevels(modDat[["rowId"]]))
      ## Fit and return the model.
      SpATS::SpATS(response = trait, fixed = fixedForm,
                   random = ~ colId + rowId,
                   spatial = ~ SpATS::PSANOVA(colNum, rowNum, nseg = nseg,
                                              nest.div = c(2, 2)),
                   genotype = genoCol, genotype.as.random = TRUE,
                   geno.decomp = geno.decomp, data = modDat,
                   control = list(maxit = 50, tolerance = 1e-03,
                                  monitoring = 0))
    })
  } else if (engine == "asreml") {
    ## fixed in asreml needs response variable on lhs of formula.
    fixedForm <- update(fixedForm, paste(trait, "~ ."))
    ## Construct formula for random part of the model.
    randForm <- formula(paste("~ colId + ", if (is.null(geno.decomp)) genoCol else
      paste0("at(", geno.decomp, "):genotype")))
    ## Loop on timepoint to run asreml.
    fitMods <- lapply(X = TP, function(timePoint) {
      message(timePoint[["timePoint"]][1])
      ## Only keep columns needed for analysis.
      modDat <- droplevels(timePoint)
      ## Run model.
      asrFit <- asreml::asreml(fixed = update(fixedForm, paste(trait, "~ .")),
                               random = randForm, data = modDat, trace = FALSE,
                               na.action = asreml::na.method(x = "include"),
                               maxiter = 200)
      ## evaluate call terms in mr and mfTrait so predict can be run.
      asrFit$call$fixed <- eval(asrFit$call$fixed)
      asrFit$call$random <- eval(asrFit$call$random)
      asrFit$call$data <- substitute(modDat)
      return(asrFit)
    })
  }
  return(createFitMod(fitMods))
}



