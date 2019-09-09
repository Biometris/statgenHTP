#' Correct for spatial effects and other unnecesary factors.
#' @keywords internal
correctSpatial <- function(fitMod) {
  ## Get engine from fitted model.
  engine <- class(fitMod)
  ## All steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (engine == "SpATS") {
    pred <- correctSpatialSpATS(fitMod)
  } else if (engine == "asreml") {
    pred <- correctSpatialAsreml(fitMod)
  }
  ## return results.
  return(pred)
}

#' @keywords internal
correctSpatialSpATS <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
  ## Get trait from fitted model.
  trait <- fitMod$model$response
  ## Set name for new trait.
  newTrait <- paste0(trait, "_corr")
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
  if (!is.null(geno.decomp)) {
    if (!hasName(x = pred, name = "geno.decomp.y")) {
      pred[["geno.decomp.y"]] <- pred[["geno.decomp"]]
    }
    ## Genotype was converted to an interaction term of genotype and
    ## geno.decomp in the proces of fitting the model. That needs to be
    ## undone to get the genotype back in the output again.
    genoStart <- nchar(as.character(pred[["geno.decomp.y"]])) + 2
    pred[["genotype.y"]] <- as.factor(substring(pred[["genotype.y"]],
                                                first = genoStart))

  }
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
  pred[[newTrait]] <- pred[[trait]] - pred[["predicted.values"]] +
    intercept
  ## Select the variables needed for subsequent analyses.
  if (!useCheck) {
    pred[["genotype"]] <- pred[["genotype.y"]]
  }
  if (!is.null(geno.decomp) && !hasName(pred , "geno.decomp")) {
    pred[[geno.decomp]] <- pred[[paste0(geno.decomp, ".y")]]
  }
  pred <- pred[c(newTrait, "genotype", geno.decomp,
                 setdiff(predVars, c("rowNum", "colNum")), "plotId",
                 "timePoint")]
}

#' @keywords internal
correctSpatialAsreml <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- "check" %in% all.vars(update(fitMod$formulae$fixed, 0~.))
  ## Get trait from fitted model.
  trait <- all.vars(update(fitMod$formulae$fixed, .~0))
  ## Set name for new trait.
  newTrait <- paste0(trait, "_corr")
  ## Get geno.decomp from fitted models.
  if ("geno.decomp" %in% all.vars(fitMod$formulae$random) ||
      "geno.decomp" %in% all.vars(fitMod$formulae$fixed)) {
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
  fixVars <- setdiff(fixVars, c("genotype", "genoCheck", geno.decomp,
                                if (useCheck) "check"))
  randVars <- all.vars(fitMod$formulae$random)
  randVars <- setdiff(randVars, c("genotype", "genoCheck", geno.decomp))
  predVars <- c(fixVars, randVars)
  pred <- fitMod$call$data[union(c("genotype", if (useCheck) "check",
                                   "plotId", "timePoint", trait,
                                   geno.decomp), c(randVars, fixVars))]
  pred[[newTrait]] <- pred[[trait]]
  ## Predict fixed effects.
  if (length(fixVars) > 0) {
    predFix <- predictAsreml(fitMod, classify = paste0(fixVars, collapse = ":"),
                             vcov = FALSE, present = fixVars)$pvals
    pred <- merge(pred, predFix[c(fixVars, "predicted.value")])
    pred[[newTrait]] <- pred[[newTrait]] - pred[["predicted.value"]]
    pred[["predicted.value"]] <- NULL
  }
  ## Predict random effects.
  if (length(randVars) > 0) {
    predRand <- predictAsreml(fitMod,
                              classify = paste0(randVars, collapse = ":"),
                              vcov = FALSE, only = randVars)$pvals
    pred <- merge(pred, predRand[c(randVars, "predicted.value")])
    pred[[newTrait]] <- pred[[newTrait]] - pred[["predicted.value"]]
    pred[["predicted.value"]] <- NULL
  }
  ## Predict intercept.
  if (length(fixVars) > 0) {
    if (!is.null(geno.decomp)) {
      predGD <- predictAsreml(fitMod, classify = "geno.decomp",
                              vcov = FALSE)$pvals
      pred <- merge(pred, predGD[c("geno.decomp", "predicted.value")])
      pred[[newTrait]] <- pred[[newTrait]] + pred[["predicted.value"]]
      pred[["predicted.value"]] <- NULL
    } else {
      predInt <- predictAsreml(fitMod, classify = "(Intercept)",
                               vcov = FALSE,
                               present = predVars)$pvals$predicted.value
      pred[[newTrait]] <- pred[[newTrait]] + predInt
    }
  }
  if (length(predVars) == 0) {
    warning("No spatial or fixed effects to correct for. Returning raw data.\n")
  }
  ## Remove row/col combinations added when fitting models.
  pred <- pred[!is.na(pred[["plotId"]]), ]
  ## Select the variables needed for subsequent analyses.
  pred <- pred[c(newTrait, "genotype", geno.decomp, fixVars, randVars,
                 "plotId", "timePoint")]
}
