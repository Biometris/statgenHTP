#' Fit Splines
#'
#' Function for fitting splines.
#'
#' @param corrDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param n Just here for testing
#' @param knots The number of knots to use when fitting the spline.
#'
#' @export
fitSpline <- function(corrDat,
                      trait,
                      n = nlevels(corrDat[["plotId"]]),
                      knots = 50) {
  ## Create data.frame with plants and genotypes for adding genotype to results.
  plantGeno <- unique(corrDat[c("plotId", "genotype")])
  ## determine minimum number of time points requiered.
  minTP <- 0.8 * length(unique(corrDat[["timeNumber"]]))
  ## Construct formula for fitting model.
  modForm <- as.formula(paste(trait, "~s(timeNumber, bs = 'ps', k = ",
                              knots, ")"))
  ## Fit splines.
  fitSp <- lapply(X = levels(corrDat[["plotId"]])[1:n], FUN = function(plant) {
    ## Restrict data to current plant.
    dat <- corrDat[corrDat[["plotId"]] == plant, c("timeNumber", trait)]
    dat <- droplevels(dat)
    ## Manually select minimum number of time points.
    if (length(unique(dat[["timeNumber"]])) >= minTP) {
      ## Fit the P-spline using gam() function in mgcv.
      ## Manually set the number of knots
      ## Depends on the nb of time points and shape of the curve)
      obj <- mgcv::gam(modForm, data = dat, method = "REML")
      ## Extract the spline coefficients.
      coeff <- data.frame(obj$coefficients, plotId = plant)
      ## Predictions on a dense grid.
      ## check how to set the grid.
      x0 = seq(from = min(dat[["timeNumber"]]), to = max(dat[["timeNumber"]]),
               by = 0.1)
      dfX0 = data.frame(timeNumber = x0)
      yPred = predict(obj, newdata = dfX0)
      predDat <- data.frame(timeNumber = x0, pred.value = yPred, plotId = plant)
      return(list(coeff, predDat))
    }
  })
  ## Bind all coefficients into one data.frame.
  coefTot <- do.call(rbind, lapply(fitSp, `[[`, 1))
  ## Remove brackets in coeffient names.
  coefTot[["type"]] <- gsub(pattern = "[()]", replacement = "",
                            x = row.names(coefTot))
  ## Merge genotype.
  coefTot <- merge(coefTot, plantGeno)
  ## Bind all predictions into one data.frame.
  predTot <- do.call(rbind, lapply(fitSp, `[[`, 2))
  ## Merge genotype.
  predTot <- merge(predTot, plantGeno)
  return(list(coefTot, predTot))
}
