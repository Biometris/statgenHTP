#' Fit Splines
#'
#' Function for fitting splines.
#'
#' @param corrDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param genotypes A character vector indicating the genotypes for which
#' splines are fitted. If \code{NULL}, splines will be fitted for all genotypes.
#' @param plotIds A character vector indicating the plotIds for which splines
#' are fitted. If \code{NULL}, splines will be fitted for all plotIds.
#' @param knots The number of knots to use when fitting the spline.
#' @param useTimeNumber Should the timeNumber be used instead of the timePoint?
#' @param timeNumber If \code{useTimeNumber = TRUE}, a character vector indicating the column
#' containing the numerical time to use.
#' @param perMinTP The percentage of minimum number of time points to use.
#'
#' @export
fitSpline <- function(corrDat,
                      trait,
                      genotypes = NULL,
                      plotIds = NULL,
                      knots = 50,
                      useTimeNumber = FALSE,
                      timeNumber = NULL,
                      perMinTP = 0.8) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(corrDat, "data.frame")) {
    stop("corrDat should be a data.frame.\n")
  }
  if (useTimeNumber && (is.null(timeNumber) || !is.character(timeNumber) ||
                        length(timeNumber) > 1)) {
    stop("timeNumber should be a character string of length 1.\n")
  }
  fitLevel <- if (hasName(x = corrDat, name = "plotId")) "plotId" else "genotype"
  corrCols <- c("genotype", trait, if (fitLevel == "plotId") "plotId",
                if (useTimeNumber) timeNumber else "timePoint")
  if (!all(hasName(x = corrDat, name = corrCols))) {
    stop("corrDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!is.null(genotypes) &&
      (!is.character(genotypes) && !all(genotypes %in% corrDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in corrDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) && !all(plotIds %in% corrDat[["genotype"]]))) {
    stop("plotIds should be a character vector of plotIds in corrDat.\n")
  }
  if (!is.numeric(knots) || length(knots) > 1 || knots < 0) {
    stop("knots should be a positive numerical value.\n")
  }
  if (!is.numeric(perMinTP) || length(perMinTP) > 1 || perMinTP < 0 ||
      perMinTP > 1) {
    stop("perMinTP should be a numerical value between 0 and 1.\n")
  }
  if (!useTimeNumber) {
    corrDat[["timeNumber"]] <- as.numeric(corrDat[["timePoint"]])
  } else {
    if (!is.numeric(corrDat[[timeNumber]])) {
      stop("timeNumber should be a numerical column.\n")
    }
    corrDat[["timeNumber"]] <- corrDat[[timeNumber]]
  }
  ## Restrict corrDat to selected genotypes and plotIds.
  if (!is.null(genotypes)) {
    corrDat <- corrDat[corrDat[["genotype"]] %in% genotypes, ]
  }
  if (!is.null(plotIds)) {
    corrDat <- corrDat[corrDat[["plotId"]] %in% plotIds, ]
  }
  if (nrow(corrDat) == 0) {
    stop("At least one valid combination of genotype and plotId should be ",
         "selected.\n")
  }
  if (fitLevel == "genotype") {
    ## If fitting at a genotype level set plotId to genotype.
    ## This way all further code can remain intact.
    corrDat[["plotId"]] <- corrDat[["genotype"]]
  }
  corrDat <- droplevels(corrDat)
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
      ## Restrict dense grid to points within observation range.
      timeRangePl <- timeRange[timeRange[["timeNumber"]] >= min(dat[["timeNumber"]]) &
                                 timeRange[["timeNumber"]] <= max(dat[["timeNumber"]]),
                               , drop = FALSE]
      ## Predictions on a dense grid.
      yPred <- predict(obj, newdata = timeRangePl)
      yDeriv <- gratia::derivatives(obj, newdata = timeRangePl,
                                    eps = 1e-7 * mean(timeNumRange[1:2]))
      ## Merge time, predictions and plotId.
      predDat <- data.frame(timeRangePl, pred.value = yPred,
                            deriv = yDeriv[["derivative"]],
                            plotId = plant)
      return(list(coeff, predDat))
    } else {
      return(list(coeff = NULL, predDat = NULL))
    }
  })
  ## Bind all coefficients into one data.frame.
  coefTot <- do.call(rbind, lapply(fitSp, `[[`, 1))
  ## Remove brackets in coefficient names.
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
  if (fitLevel == "genotype") {
    ## Remove plotId (duplicated genotype) from output.
    coefTot[["plotId"]] <- NULL
    predTot[["plotId"]] <- NULL
  }
  ## Create output.
  res <- structure(list(coefDat = coefTot, predDat = predTot),
                   modDat = corrDat,
                   trait = trait,
                   useTimeNumber = useTimeNumber,
                   fitLevel = fitLevel,
                   class = c("HTPSpline", "list"))
  return(res)
}

#' Plot the results of a fitted spline.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class HTPSpline.
#' @param genotypes A character vector indicating the genotypes for which
#' splines should be plotted.
#' @param plotIds A character vector indicating the plotIds for which splines
#' should be plotted.
#'
#' @export
plot.HTPSpline <- function(x,
                           ...,
                           plotType = c("predictions", "derivatives"),
                           genotypes = NULL,
                           plotIds = NULL,
                           output = TRUE) {
  plotType <- match.arg(plotType)
  plotVar <- if (plotType == "predictions") "pred.value" else "deriv"
  modDat <- attr(x, which = "modDat")
  trait <- attr(x, which = "trait")
  fitLevel <- attr(x, which = "fitLevel")
  useTimeNumber <- attr(x, which = "useTimeNumber")
  predDat <- x$predDat
  if (!is.null(genotypes) &&
      (!is.character(genotypes) && !all(genotypes %in% predDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in predDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) && !all(plotIds %in% predDat[["genotype"]]))) {
    stop("plotIds should be a character vector of plotIds in predDat.\n")
  }
  ## Restrict predDat and modDat to selected genotypes and plotIds.
  if (!is.null(genotypes)) {
    predDat <- predDat[predDat[["genotype"]] %in% genotypes, ]
    modDat <- modDat[modDat[["genotype"]] %in% genotypes, ]
  }
  if (!is.null(plotIds)) {
    predDat <- predDat[predDat[["plotId"]] %in% plotIds, ]
    modDat <- modDat[modDat[["plotId"]] %in% plotIds, ]
  }
  ## Remove plotIds with only NA from data.
  ## This can be caused by removing outliers.
  fitLevels <- unique(modDat[[fitLevel]])
  allNA <- sapply(X = fitLevels, FUN = function(x) {
    all(is.na(modDat[modDat[[fitLevel]] == x, trait]))
  })
  modDat <- modDat[!modDat[[fitLevel]] %in% fitLevels[allNA], ]
  predDat <- predDat[!predDat[[fitLevel]] %in% fitLevels[allNA], ]
  modDat <- droplevels(modDat)
  predDat <- droplevels(predDat)
  if (nrow(predDat) == 0) {
    stop("At least one valid combination of genotype and plotId should be ",
         "selected.\n")
  }
  timeVar <- if (useTimeNumber) "timeNumber" else "timePoint"
  p <- ggplot(modDat, aes_string(x = timeVar, y = trait)) +
    geom_line(data = predDat,
              aes_string(x = timeVar, y = plotVar), col = "blue", na.rm = TRUE) +
    labs(y = trait, x = timeVar)
  if (plotType == "predictions") {
    p <- p + geom_point(na.rm = TRUE) +
      ggtitle("Corrected data and Pspline prediction")
  } else {
    p <- p + ggtitle("Pspline first derivatives")
  }
  if (!useTimeNumber) {
    ## Compute the number of breaks for the time scale.
    ## If there are less than 3 time points use the number of time points.
    ## Otherwise use 3.
    nBr <- min(length(unique(modDat[["timePoint"]])), 3)
    ## Format the time scale to Month + day.
    p <- p + ggplot2::scale_x_datetime(breaks = prettier(n = nBr),
                                       labels = scales::date_format("%B %d"))
  }
  ## Calculate the total number of plots.
  nPlots <- length(unique(modDat[[fitLevel]]))
  ## 25 Plots per page.
  nPag <- ceiling(nPlots / 25)
  if (nPlots >= 25) {
    ## More than 25 plots.
    ## For identical layout on all pages use 5 x 5 plots throughout.
    rowPag <- colPag <- rep(x = 5, times = nPag)

    # 28-7-2020. ggforce has a bug that prevents this identical layout
    # https://github.com/thomasp85/ggforce/issues/201
    # When fixed the code above can be reactivated and the three lines below
    # removed.
    # plotsLastPag <- nPlots %% 25
    # rowPag <- c(rep(x = 5, times = nPag - 1), min(plotsLastPag %/% 5 + 1, 5))
    # colPag <- c(rep(x = 5, times = nPag - 1),
    #             ifelse(plotsLastPag >= 5, 5, plotsLastPag))
  } else {
    ## Less than 25 plots.
    ## Fill page by row of 5 plots.
    plotsPag <- nPlots %% 25
    rowPag <- min(plotsPag %/% 5 + 1, 5)
    colPag <- ifelse(plotsPag >= 5, 5, plotsPag)
  }
  ## Build pages of plots.
  pPag <- vector(mode = "list", length = nPag)
  for (i in 1:nPag) {
    pPag[[i]] <- p +
      ggforce::facet_wrap_paginate(facets = fitLevel, nrow = rowPag[i],
                                   ncol = colPag[i], page = i)
    if (output) {
      suppressMessages(plot(pPag[[i]]))
    }
  }
  invisible(pPag)
}

#' Extract estimates from fitted splines.
#'
#' Function for extracting parameter estimates from fitted splines on a
#' specified interval.
#'
#' @param HTPSpline An object of class HTPSpline, the output of the
#' \code{\link{fitSpline}} function.
#' @param estimate The type of estimate that should be extracted,
#' the predictions or the first derivatives.
#' @param what The type of estimate that should be extracted.
#' @param timeMin The lower bound of the time interval from which the
#' estimates should be extracted. If \code{NULL} the smallest time value for
#' which the splines were fitted is used.
#' @param timeMax The upper bound of the time interval from which the
#' estimates should be extracted. If \code{NULL} the largest time value for
#' which the splines were fitted is used.
#' @param genotypes A character vector indicating the genotypes for which
#' splines are fitted. If \code{NULL}, splines will be fitted for all genotypes.
#' @param plotIds A character vector indicating the plotIds for which splines
#' are fitted. If \code{NULL}, splines will be fitted for all plotIds.
#'
#' @export
estimateSplineParameters <- function(HTPSpline,
                                     estimate = c("predictions", "derivatives"),
                                     what = c("min", "max", "mean"),
                                     timeMin = NULL,
                                     timeMax = NULL,
                                     genotypes = NULL,
                                     plotIds = NULL) {
  estimate <- match.arg(estimate)
  what <- match.arg(what)
  estVar <- if (estimate == "predictions") "pred.value" else "deriv"
  useTimeNumber <- attr(HTPSpline, which = "useTimeNumber")
  fitLevel <- attr(HTPSpline, which = "fitLevel")
  predDat <- HTPSpline$predDat
  if (!is.null(genotypes) &&
      (!is.character(genotypes) && !all(genotypes %in% predDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in predDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) && !all(plotIds %in% predDat[["genotype"]]))) {
    stop("plotIds should be a character vector of plotIds in predDat.\n")
  }
  ## Restrict predDat to selected genotypes and plotIds.
  if (!is.null(genotypes)) {
    predDat <- predDat[predDat[["genotype"]] %in% genotypes, ]
  }
  if (!is.null(plotIds)) {
    predDat <- predDat[predDat[["plotId"]] %in% plotIds, ]
  }
  if (nrow(predDat) == 0) {
    stop("At least one valid combination of genotype and plotId should be ",
         "selected.\n")
  }
  timeVar <- if (useTimeNumber) "timeNumber" else "timePoint"
  if (is.null(timeMin)) {
    timeMin <- min(predDat[[timeVar]])
  } else {
    if (timeMin < min(predDat[[timeVar]]) || timeMin > max(predDat[[timeVar]])) {
      stop("timeMin should be within the time interval in the data.\n")
    }
  }
  if (is.null(timeMax)) {
    timeMax <- max(predDat[[timeVar]])
  } else {
    if (timeMax < min(predDat[[timeVar]]) || timeMax > max(predDat[[timeVar]])) {
      stop("timeMax should be within the time interval in the data.\n")
    }
  }
  if (timeMin >= timeMax) {
    stop("timeMax should be larger than timeMin.\n")
  }
  ## Restrict predDat to time interval.
  predDat <- predDat[predDat[[timeVar]] >= timeMin &
                       predDat[[timeVar]] <= timeMax, ]
  ## Get estimates.
  res <- aggregate(x = predDat[[estVar]],
                   by = predDat[c("genotype",
                                  if (fitLevel == "plotId") "plotId")],
                   FUN = what)
  return(res)
}
