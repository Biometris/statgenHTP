#' Approach 2: correct for spatial effects and other unnecesary factors
#' @keywords internal
correctSpatial <- function(fitMod) {
  ## Get engine from fitted model.
  engine <- class(fitMod)
  ## All steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (engine == "SpATS") {
    pred <- correctSpatialSpATS(fitMod)
  } else {
    pred <- correctSpatialAsreml(fitMod)
  }
  ## return results.
  return(pred)
}

correctSpatialSpATS <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
  ## Get trait from fitted model.
  trait <- fitMod$model$response
  ## Get geno.decomp from fitted models.
  geno.decomp <- fitMod$model$geno$geno.decomp
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing.
  if (!is.null(fitMod$model$fixed)) {
    fixVars <- attr(terms(fitMod$model$fixed), "term.labels")
  } else {
    fixVars <- NULL
  }
  predVars <- setdiff(c(fixVars, "colNum", "rowNum", "colId", "rowId"),
                      c(geno.decomp, "check"))
  pred <- predict(fitMod, which = predVars)
  ## Merge genotype and timepoint to data
  pred <- merge(pred, fitMod$data[c("rowNum", "colNum", "genotype",
                                    "plotId", "timePoint", trait,
                                    geno.decomp)],
                by = c("rowNum", "colNum"))
  ## Temporary fix for difference between SpATS and asreml predictions.
  ## asreml predicts marginal means whereas SpATS predicts conditional means.
  ## By adding the means of the fixed effects to the conditional means the
  ## marginal means are calculated.
  ## Note that this means the standard errors are no longer correct.
  corVars <- setdiff(fixVars, c(geno.decomp, "check"))
  intercept <- fitMod$coeff["Intercept"]
  if (length(corVars) > 0) {
    ## Order in descreasing order so variables that are substrings of other
    ## variables are treated correctly.
    corVars <- corVars[order(nchar(corVars), decreasing = TRUE)]
    ## Get coefficients for fixed variables.
    coeffs <- fitMod$coeff[!attr(fitMod$coeff, "random")]
    ## Loop over corVars and adjust predicted value by mean of fixed effects
    ## for corVar. Then remove it from coeff so it isn't used again by a
    ## shorter variable, i.e. repId1 and repId
    for (corVar in corVars) {
      corMean <- mean(c(0, coeffs[grepl(corVar, names(coeffs))]))
      intercept <- intercept + corMean
      coeffs <- coeffs[!grepl(corVar, names(coeffs))]
    }
  }
  ## Obtain the corrected trait.
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.values"]] +
    intercept
  ## Select the variables needed for subsequent analyses.
  if (!useCheck) {
    pred[["genotype"]] <- pred[["genotype.y"]]
  }
  if (!is.null(geno.decomp) && !hasName(pred , "geno.decomp")) {
    pred[[geno.decomp]] <- pred[[paste0(geno.decomp, ".y")]]
  }
  pred <- pred[c("newTrait", "genotype", geno.decomp, predVars, "plotId",
                 "timePoint")]
}

correctSpatialAsreml <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- "check" %in% all.vars(update(fitMod$formulae$fixed, 0~.))
  ## Get trait from fitted model.
  trait <- all.vars(update(fitMod$formulae$fixed, .~0))
  ## Get geno.decomp from fitted models.
  if ("geno.decomp" %in% all.vars(fitMod$formulae$random)) {
    geno.decomp <- "geno.decomp"
  } else {
    geno.decomp <- NULL
  }
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing.
  if (!is.null(fitMod$formulae$fixed)) {
    fixVars <- all.vars(update(fitMod$formulae$fixed, 0~.))
  } else {
    fixVars <- NULL
  }
  randVars <- all.vars(fitMod$formulae$random)
  predVars <- setdiff(c(fixVars, randVars), c("genotype", "genoCheck",
                                              if (useCheck) "check",
                                              geno.decomp))
  if (length(predVars) > 0) {
    pred <- predict(fitMod, classify = paste(predVars, collapse = "+"),
                    present = c(predVars, geno.decomp))$pvals
    ## Merge genotype and timepoint to data
    pred <- merge(pred, fitMod$call$data[union(c("genotype",
                                                 if (useCheck) "check",
                                                 "plotId", "timePoint", trait,
                                                 geno.decomp), predVars)],
                  by = predVars)
    if (!is.null(geno.decomp)) {
      predGD <- predict(fitMod, classify = "geno.decomp")$pvals
      pred <- merge(pred, predGD, by = geno.decomp)
    } else {
      predInt <- predict(fitMod, classify = "(Intercept)",
                         present = predVars)$pvals
      pred[["predicted.value.x"]] <- pred[["predicted.value"]]
      pred[["predicted.value.y"]] <- predInt$predicted.value
    }
    ## Obtain the corrected trait.
    pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.value.x"]] +
      pred[["predicted.value.y"]]
  } else {
    ## Nothing to correct for. Return raw values.
    pred <- fitMod$call$data[c("genotype", if (useCheck) "check", "plotId",
                               "timePoint", trait, geno.decomp)]
    pred[["newTrait"]] <- pred[[trait]]
    warning("No spatial or fixed effects to correct for. Returning raw data.\n")
  }
  ## Select the variables needed for subsequent analyses.
  pred <- pred[c("newTrait", "genotype", geno.decomp, predVars, "plotId",
                 "timePoint")]
}
