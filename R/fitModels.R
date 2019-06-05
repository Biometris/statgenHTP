#' @export
fitModels <- function(TP,
                      trait,
                      covariates = NULL,
                      geno.decomp = NULL,
                      useCheck = FALSE) {
  ## Checks
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  if (!all(sapply(X = TP, FUN = hasName, name = trait))) {
    stop(trait, " should be a column for all timePoints.\n")
  }
  if (!is.null(covariates)) {
    for (covar in covariates) {
      if (!all(sapply(X = TP, FUN = hasName, name = covar))) {
        stop(covar, " should be a column for all timePoints.\n")
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
  ## Construct formula for fixed part.
  ## Fixed part consists of covariates, geno.decomp and check.
  if (!is.null(c(covariates, geno.decomp))) {
    fixedForm <- formula(paste("~", paste(c(covariates, geno.decomp),
                                          collapse = "+"),
                               if (useCheck) "+ check"))
  } else {
    fixedForm <- if (useCheck) formula("~ check") else NULL
  }
  ## Loop on timepoint to run SpATS.
  fitMods <- lapply(X = TP, function(timePoint) {
    message(timePoint[["timePoint"]][1])
    ## Only keep columns needed for analysis.
    modCols <- c("timePoint", "plotId", "genotype", "genoCheck", "check",
                 "colId", "rowId", "colNum", "rowNum", covariates, geno.decomp,
                 trait)
    modDat <- timePoint[colnames(timePoint) %in% modCols]
    modDat <- droplevels(modDat)
    ## number of segments for SpATS.
    nseg = c(nlevels(modDat[["colId"]]), nlevels(modDat[["rowId"]]))
    ## Fit and return the model.
    SpATS::SpATS(response = trait, fixed = fixedForm,
                 random = ~ colId + rowId,
                 spatial = ~ SpATS::PSANOVA(colNum, rowNum, nseg = nseg,
                                            nest.div = c(2, 2)),
                 genotype = if (useCheck) "genoCheck" else "genotype",
                 genotype.as.random = TRUE, geno.decomp = geno.decomp,
                 data = modDat,
                 control = list(maxit = 50, tolerance = 1e-03, monitoring = 0))
  })
  return(createFitMod(fitMods))
}
