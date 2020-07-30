#' Fit Splines
#'
#' Function for fitting splines.
#'
#' @param corrDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param knots The number of knots to use when fitting the spline.
#' @param perMinTP The percentage of minimum number of time points to use.
#' @param timeNumStep Time step of the time scale for the prediction.
#'
#' @export
fitSpline <- function(corrDat,
                      trait,
                      knots = 50,
                      perMinTP = 0.8,
                      timeNumStep = 0.1) {
  ## Create data.frame with plants and genotypes for adding genotype to results.
  plantGeno <- unique(corrDat[c("plotId", "genotype")])
  ## Determine minimum number of time points required.
  minTP <- perMinTP * length(unique(corrDat[["timePoint"]])) # 0.8
  ## Construct formula for fitting model.
  modForm <- as.formula(paste(trait, "~s(timeNumber, bs = 'ps', k = ",
                              knots, ")"))
  ## Create time range for prediction.
  ## Get range for time number and time point from data.
  timeNumRange <- range(corrDat[["timeNumber"]])
  timePointRange <- range(corrDat[["timePoint"]])
  ## Create data.frame with time number and time point on prediction scale.
  timeRange <- data.frame(timeNumber = seq(from = timeNumRange[1],
                                           to = timeNumRange[2],
                                           by = timeNumStep),
                          timePoint = seq(from = timePointRange[1],
                                          to = timePointRange[2],
                                          by = timeNumStep * diff(timePointRange) /
                                            diff(timeNumRange)))
  ## Fit splines.
  fitSp <- lapply(X = levels(plantGeno[["plotId"]]), FUN = function(plant) {
    ## Restrict data to current plant.
    dat <- corrDat[corrDat[["plotId"]] == plant, c("timeNumber", trait)]
    ## Manually select minimum number of time points.
    if (length(unique(dat[["timeNumber"]])) >= minTP) {

      ## Fit the P-spline using gam() function in mgcv.
      ## Manually set the number of knots.
      ## Depends on the number of time points and shape of the curve.
      obj <- mgcv::gam(modForm, data = dat, method = "REML")
      ## Extract the spline coefficients.
      coeff <- data.frame(obj$coefficients, plotId = plant)
      coeff$type <- row.names(coeff)
      ## Predictions on a dense grid.
      dfX0 <- timeRange[, "timeNumber", drop = FALSE]
      yPred <- predict(obj, newdata = dfX0)
      ## Merge time, predictions and plotId.
      predDat <- data.frame(timeRange, pred.value = yPred, plotId = plant)
      return(list(coeff, predDat))
    } else {
      return(list(coeff = NULL, predDat = NULL))
    }
  })
  ## Bind all coefficients into one data.frame.
  coefTot <- do.call(rbind, lapply(fitSp, `[[`, 1))
  ## Remove brackets in coeffient names.
  coefTot[["type"]] <- gsub(pattern = "[()]", replacement = "",
                            x =  coefTot[["type"]])
  ## Add genotype.
  coefTot[["genotype"]] <- plantGeno[match(coefTot[["plotId"]],
                                           plantGeno[["plotId"]]), "genotype"]
  ## Bind all predictions into one data.frame.
  predTot <- do.call(rbind, lapply(fitSp, `[[`, 2))
  ## Add genotype.
  predTot[["genotype"]] <- plantGeno[match(predTot[["plotId"]],
                                           plantGeno[["plotId"]]), "genotype"]
  return(list(coefDat = coefTot, predDat = predTot))
}
