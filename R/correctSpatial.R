#' Approach 2: correct for spatial effects and other unnecesary factors
#' @keywords internal
correctSpatial <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
  ## Get trait from fitted model.
  trait <- fitMod$model$response
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing.
  if (!is.null(fitMod$model$fixed)) {
    fixVars <- attr(terms(fitMod$model$fixed), "term.labels")
  } else {
    fixVars <- NULL
  }
  genoDec <- fitMod$model$geno$geno.decomp
  predVars <- setdiff(c(fixVars, "colNum", "rowNum", "colId", "rowId"),
                      genoDec)
  pred <- predict(fitMod, which = predVars)
  ## Merge genotype, pos, time and timepoint to data
  pred <- merge(pred, fitMod$data[c("rowNum", "colNum", "genotype",
                                       "plotId", "timePoint", trait,
                                       genoDec)],
                by = c("rowNum", "colNum"))
  ## Obtain the corrected trait.
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.values"]] +
    fitMod$coeff["Intercept"]
  ## Select the variables needed for subsequent analyses.
  if (!useCheck) {
    pred[["genotype"]] <- pred[["genotype.y"]]
  }
  if (!is.null(genoDec)) pred[[genoDec]] <- pred[[paste0(genoDec, ".y")]]
  pred <- pred[c("newTrait", "genotype", genoDec, predVars, "plotId",
                 "timePoint")]
  ## return results.
  return(pred)
}
