fitModels <- function(TD,
                      trait,
                      covariates = NULL,
                      geno.decomp = NULL) {
  # In these analyses, the genotype is always included as random
  # I have to check with Coté if we can also include it as fixed
  genotype.as.random <- TRUE
  ## Construct formula for fixed part.
  if (!is.null(covariates)) {
  fixedForm <- formula(paste("~", paste(covariates, collapse = "+")))
  } else {
    fixedForm <- NULL
  }
  ##############################################################################
  ######## Loop on timepoint to run SpATS

  fitMods <- lapply(X = seq_along(TD), function(i) {
    message(names(TD)[i])
    modDat <- droplevels(TD[[i]])
    ### number of segments for SpATS:
    nseg = c(nlevels(modDat[["colId"]]), nlevels(modDat[["rowId"]]))

    # Fit the model using check
    # fit.SpATS2 <- SpATS::SpATS(response = trait,
    #                           fixed = ~ Sowing_Block + Image_pos + Check ,
    #                           random = ~ Col + Row ,
    #                           spatial = ~ SpATS::PSANOVA(colNum, rowNum,
    #                                                      nseg = nseg,
    #                                                      nest.div=c(2,2)),
    #                           genotype = "Genobis",
    #                           genotype.as.random = genotype.as.random,
    #                           data = dat.ti,
    #                           control = list(maxit = 50,
    #                                          tolerance = 1e-03,
    #                                          monitoring = 0))

    # Fit the model without the check effect
    SpATS::SpATS(response = trait, fixed = fixedForm,
                 random = ~ colId + rowId,
                 spatial = ~ SpATS::PSANOVA(colNum, rowNum, nseg = nseg,
                                            nest.div = c(2, 2)),
                 genotype = "genotype", genotype.as.random = genotype.as.random,
                 geno.decomp = geno.decomp, data = modDat,
                 control = list(maxit = 50, tolerance = 1e-03, monitoring = 0))

  })

  return(fitMods)
}
