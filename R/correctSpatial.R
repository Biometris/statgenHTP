#' Correct for spatial effects and other unnecesary factors.
#' @keywords internal
correctSpatial <- function(fitMod) {
  ## All steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (inherits(fitMod, "SpATS")) {
    pred <- correctSpatialSpATS(fitMod)
  } else if (inherits(fitMod, "asreml")) {
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
  pred <- predict(fitMod, which = predVars, predFixed = "marginal")
  ## Merge genotype and timepoint to data
  pred <- merge(pred, fitMod$data[c("rowNum", "colNum", "genotype",
                                    "plotId", "timePoint", trait,
                                    geno.decomp)],
                by = c("rowNum", "colNum"))
  if (!is.null(geno.decomp) && !useCheck) {
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
  ## Predict intercept.
  if (!is.null(geno.decomp) && !useCheck) {
    predGD <- predict(fitMod, which = "geno.decomp")
    intercept <- mean(predGD[["predicted.values"]])
  } else {
    intercept <- fitMod$coeff["Intercept"]
  }
  ## Obtain the corrected trait.
  pred[[newTrait]] <- pred[[trait]] - pred[["predicted.values"]] + intercept
  ## Select the variables needed for subsequent analyses.
  if (!useCheck) {
    pred[["genotype"]] <- pred[["genotype.y"]]
  }
  if (!is.null(geno.decomp) && !hasName(pred , "geno.decomp")) {
    pred[[geno.decomp]] <- pred[[paste0(geno.decomp, ".y")]]
  }
  pred <- pred[c(newTrait, trait, "genotype", geno.decomp,
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
  fixVars <- setdiff(fixVars, c("genotype", "genoCheck", geno.decomp))
  randVars <- all.vars(fitMod$formulae$random)
  randVars <- setdiff(randVars, c("genotype", "genoCheck", geno.decomp))
  predVars <- c(fixVars, randVars)
  pred <- fitMod$call$data[union(c("genotype", if (useCheck) "check",
                                   "plotId", "timePoint", trait,
                                   geno.decomp), predVars)]
  pred[[newTrait]] <- pred[[trait]]
  ## Predict fixed + random effects.
  if (length(predVars) > 0) {
    predFix <- predictAsreml(fitMod,
                             classify = paste0(predVars, collapse = ":"),
                             vcov = FALSE, present = predVars)$pvals
    pred <- merge(pred, predFix[c(predVars, "predicted.value")])
    pred[[newTrait]] <- pred[[newTrait]] - pred[["predicted.value"]]
    pred[["predicted.value"]] <- NULL
  }
  ## Predict intercept.
  if (length(predVars) > 0) {
    if (!is.null(geno.decomp)) {
      predGD <- predictAsreml(fitMod, classify = "geno.decomp",
                              present = c("geno.decomp", predVars),
                              vcov = FALSE)$pvals
      intercept <- mean(predGD[["predicted.value"]])
    } else {
      ## No predict here. Just get the intercept from the coefficients.
      intercept <- fitMod$coefficients$fixed["(Intercept)", ]
    }
    pred[[newTrait]] <- pred[[newTrait]] + intercept
  }
  if (length(predVars) == 0) {
    warning("No spatial or fixed effects to correct for. Returning raw data.\n")
  }
  ## Remove row/col combinations added when fitting models.
  pred <- pred[!is.na(pred[["plotId"]]), ]
  ## Select the variables needed for subsequent analyses.
  pred <- pred[c(newTrait, trait, "genotype", geno.decomp, fixVars, randVars,
                 "plotId", "timePoint")]
}
