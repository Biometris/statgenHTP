method2 <- function(fit.SpATS) {
  ##############################################################################
  # Approach 2: correct for spatial effects and other unnecesary factors
  ##############################################################################

  ## Get trait from fitted model.
  trait <- fit.SpATS$model$response
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing
  pred <- predict(fit.SpATS, which = c("Colnum", "Rownum", "Col", "Row",
                                       "Sowing_Block", "Image_pos"))
  ## Merge genotype, pos and time to data
  pred <- merge(pred, fit.SpATS$data[c("Rownum", "Colnum", "Geno",
                                       "pos", "Time", trait)],
                by = c("Rownum", "Colnum"))
  # Obtain the corrected trait
  pred[["newTrait"]] <- pred[[trait]] - pred[["predicted.values"]] +
    fit.SpATS$coeff["Intercept"]
  # Select the needed variables for subsequent analyses
  pred[["Genotype"]] <- pred[["Geno.y"]]
  pred <- pred[c("newTrait", "Genotype", "Time", "Sowing_Block", "Image_pos",
                 "Colnum", "Rownum", "Col", "Row", "pos")]
  ## return results
  return(pred)
}
