#' @export
fitModels <- function(TD,
                      trait,
                      covariates = NULL,
                      geno.decomp = NULL,
                      useCheck = FALSE) {
  # In these analyses, the genotype is always included as random
  # I have to check with Coté if we can also include it as fixed
  genotype.as.random <- TRUE
  ## Construct formula for fixed part.
  if (!is.null(covariates)) {
    fixedForm <- formula(paste("~", paste(covariates, collapse = "+"),
                         if (useCheck) "+ check"))
  } else {
    fixedForm <- if (useCheck) formula("~check") else NULL
  }
  ##############################################################################
  ######## Loop on timepoint to run SpATS
  fitMods <- lapply(X = TD, function(timePoint) {
    message(timePoint[["timePoint"]][1])
    modDat <- droplevels(timePoint)
    ## number of segments for SpATS:
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
