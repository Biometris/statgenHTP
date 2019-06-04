#' @export
fitModels <- function(TD,
                      trait,
                      covariates = NULL,
                      geno.decomp = NULL,
                      useCheck = FALSE) {
  # In these analyses, the genotype is always included as random
  # I have to check with Coté if we can also include it as fixed
  genotype.as.random <- TRUE
  ## All covariates should be factors.
  for (covar in covariates) {
    TD <- lapply(X = TD, FUN = function(timePoint) {
      if (!is.factor(timePoint[[covar]])) {
        timePoint[[covar]] <- as.factor(timePoint[[covar]])
      }
      return(timePoint)
    })
  }
  ## If geno.decomp is used genotype and covariates have to be replaced by
  ## an interaction of genotype and covariates with the geno.decomp variable.
  ## Add these new variables to the data.
  if (!is.null(geno.decomp)) {
    TD <- lapply(X = TD, FUN = function(timePoint) {
      timePoint[["genotype"]] <- interaction(timePoint[[geno.decomp]],
                                         timePoint[["genotype"]], sep = "_")
      for (covar in covariates) {
        timePoint[[covar]] <- interaction(timePoint[[geno.decomp]],
                                          timePoint[[covar]], sep = "_")
      }
      return(timePoint)
    })
  }
  ## Construct formula for fixed part.
  if (!is.null(covariates)) {
    fixedForm <- formula(paste("~", paste(covariates, collapse = "+"),
                         if (useCheck) "+ check"))
    if (!is.null(geno.decomp)) {
      geno.decomp <- covariates
    }
  } else {
    fixedForm <- if (useCheck) formula("~ check") else NULL
  }
  ## Loop on timepoint to run SpATS.
  fitMods <- lapply(X = TD, function(timePoint) {
    message(timePoint[["timePoint"]][1])
    ## Only keep columns needed for analysis.
    modCols <- c("timePoint", "plotId", "genotype", "genoCheck", "check",
                 "colId", "rowId", "colNum", "rowNum", covariates, trait)
    modDat <- timePoint[colnames(timePoint) %in% modCols]
    modDat <- droplevels(modDat)
    ## number of segments for SpATS.
    nseg = c(nlevels(modDat[["colId"]]), nlevels(modDat[["rowId"]]))
    ## Fit the model.
    SpATS::SpATS(response = trait, fixed = fixedForm,
                 random = ~ colId + rowId,
                 spatial = ~ SpATS::PSANOVA(colNum, rowNum, nseg = nseg,
                                            nest.div = c(2, 2)),
                 genotype = if (useCheck) "genoCheck" else "genotype",
                 genotype.as.random = genotype.as.random,
                 geno.decomp = geno.decomp, data = modDat,
                 control = list(maxit = 50, tolerance = 1e-03, monitoring = 0))
  })
  return(createFitMod(fitMods))
}
