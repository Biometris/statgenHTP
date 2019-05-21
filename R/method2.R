method2 <- function(fit.SpATS) {
  ##############################################################################
  # Approach 2: correct for spatial effects and other unnecesary factors
  ##############################################################################

  ## Get trait from fitted model.
  trait <- fit.SpATS$model$response
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing
  pred <- predict(fit.SpATS, which = c("colNum", "rowNum", "colId", "rowId",
                                       "Sowing_Block", "Image_pos"))
  ## Merge genotype, pos, time and timepoint to data
  pred <- merge(pred, fit.SpATS$data[c("rowNum", "colNum", "genotype",
                                       "pos", "time", "timePoint", trait)],
                by = c("rowNum", "colNum"))
  # Obtain the corrected trait
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.values"]] +
    fit.SpATS$coeff["Intercept"]
  # Select the needed variables for subsequent analyses
  pred[["genotype"]] <- pred[["genotype.y"]]
  pred[["time"]] <- as.numeric(pred[["time"]])
  pred <- pred[c("newTrait", "genotype", "time", "Sowing_Block", "Image_pos",
                 "colNum", "rowNum", "colId", "rowId", "pos", "timePoint")]
  ## return results
  return(pred)
}
