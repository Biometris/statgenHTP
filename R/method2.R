method2 <- function(fit.SpATS) {
  ##############################################################################
  # Approach 2: correct for spatial effects and other unnecesary factors
  ##############################################################################

  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fit.SpATS$model$fixed))
  ## Get trait from fitted model.
  trait <- fit.SpATS$model$response
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing.
  if (!is.null(fit.SpATS$model$fixed)) {
    fixVars <- attr(terms(fit.SpATS$model$fixed), "term.labels")
  } else {
    fixVars <- NULL
  }
  genoDec <- fit.SpATS$model$geno$geno.decomp
  predVars <- setdiff(c(fixVars, "colNum", "rowNum", "colId", "rowId"),
                      genoDec)
  pred <- predict(fit.SpATS, which = predVars)
  ## Merge genotype, pos, time and timepoint to data
  pred <- merge(pred, fit.SpATS$data[c("rowNum", "colNum", "genotype",
                                       "plotId", "time", "timePoint", trait,
                                       genoDec)],
                by = c("rowNum", "colNum"))
  ## Obtain the corrected trait.
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.values"]] +
    fit.SpATS$coeff["Intercept"]
  ## Select the variables needed for subsequent analyses.
  if (!useCheck) {
    pred[["genotype"]] <- pred[["genotype.y"]]
  }
  if (!is.null(genoDec)) pred[[genoDec]] <- pred[[paste0(genoDec, ".y")]]
  pred[["time"]] <- as.numeric(pred[["time"]])
  pred <- pred[c("newTrait", "genotype", "time", genoDec, predVars, "plotId",
                 "timePoint")]
  ## return results.
  return(pred)
}
