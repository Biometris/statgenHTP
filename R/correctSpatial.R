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
  ## Obtain the corrected trait.
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.values"]] +
    fitMod$coeff["Intercept"]
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
  pred <- predict(fitMod, classify = paste(predVars, collapse = "+"),
                  present = predVars)$pvals
  ## Merge genotype and timepoint to data
  pred <- merge(pred, fitMod$call$data[union(c("genotype", "check",
                                               "plotId", "timePoint", trait,
                                               geno.decomp), predVars)],
                by = predVars)
  if (!is.null(geno.decomp)) {
    predGD <- predict(fitMod, classify = "geno.decomp")$pvals
    pred <- merge(pred, predGD, by = geno.decomp)

    #predGD <- predict(fitMod, classify = "geno.decomp + check")$pvals
    #pred <- merge(pred, predGD, by = c("geno.decomp", "check"))

  } else {
    predInt <- predict(fitMod, classify = "(Intercept)",
                       present = fixVars)$pvals
    pred[["predicted.value.x"]] <- pred[["predicted.value"]]
    pred[["predicted.value.y"]] <- predInt$predicted.value
  }
  ## Obtain the corrected trait.
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.value.x"]] +
    pred[["predicted.value.y"]]
  ## Select the variables needed for subsequent analyses.
  pred <- pred[c("newTrait", "genotype", geno.decomp, predVars, "plotId",
                 "timePoint")]
}
