#' Fit Splines
#'
#' Function for fitting splines.
#'
#' @param corrDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param knots The number of knots to use when fitting the spline.
#' @param perMinTP The percentage of minimum number of time points to use.
#'
#' @export
fitSpline <- function(corrDat,
                      trait,
                      knots = 50,
                      perMinTP = 0.8) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(corrDat, "data.frame")) {
    stop("corrDat should be a data.frame.\n")
  }
  corrCols <- c("plotId", "genotype", trait)
  if (!all(hasName(x = corrDat, name = corrCols))) {
    stop("corrDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!is.numeric(knots) || length(knots) > 1 || knots < 0) {
    stop("knots should be a positive numerical value.\n")
  }
  if (!is.numeric(perMinTP) || length(perMinTP) > 1 || perMinTP < 0 ||
      perMinTP > 1) {
    stop("perMinTP should be a numerical value between 0 and 1.\n")
  }
  if (!hasName(x = corrDat, name = "timeNumber")) {
    if (hasName(x = corrDat, name = "timePoint")) {
      corrDat[["timeNumber"]] <- as.numeric(corrDat[["timePoint"]])
    } else {
      stop("corrDat should contain at least one of the colums timeNumber and ",
           "timePoint.\n")
    }
  }
  ## Create data.frame with plants and genotypes for adding genotype to results.
  plantGeno <- unique(corrDat[c("plotId", "genotype")])
  ## Determine minimum number of time points required.
  minTP <- perMinTP * length(unique(corrDat[["timeNumber"]]))
  ## Construct formula for fitting model.
  modForm <- as.formula(paste(trait, "~s(timeNumber, bs = 'ps', k = ",
                              knots, ")"))
  ## Compute step size for prediction grid.
  ## Use smallest time gap between two points and divide that in 10 segments.
  minStep <- min(diff(sort(unique(corrDat[["timeNumber"]]))))
  timeNumStep <- minStep / 9
  ## Create time range for prediction.
  ## Get range for time number and time point from data.
  timeNumRange <- range(corrDat[["timeNumber"]])
  ## Create data.frame with time number and, if present,
  ## time point on prediction scale.
  timeRange <- data.frame(timeNumber = seq(from = timeNumRange[1],
                                           to = timeNumRange[2],
                                           by = timeNumStep))
  if (hasName(x = corrDat, name = "timePoint")) {
    timePointRange <- range(corrDat[["timePoint"]])
    timeRange[["timePoint"]] <- seq(from = timePointRange[1],
                                    to = timePointRange[2],
                                    by = timeNumStep * diff(timePointRange) /
                                      diff(timeNumRange))
  }
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
      yPred <- predict(obj, newdata = timeRange)
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
