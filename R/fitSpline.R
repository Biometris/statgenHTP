#' Fit Splines
#'
#' Fit P-Splines on corrected or raw data. The number of
#' knots is chosen by the user. The function outputs are predicted P-Spline
#' values and their first and second derivatives on a dense grid. The
#' outputs can then be used for outlier detection for time series
#' (see \code{\link{detectSerieOut}}) and to estimate relevant parameters from
#' the curve for further analysis (see \code{\link{estimateSplineParameters}}).
#'
#' @param inDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param genotypes A character vector indicating the genotypes for which
#' splines should be fitted. If \code{NULL}, splines will be fitted for all
#' genotypes.
#' @param plotIds A character vector indicating the plotIds for which splines
#' should be fitted. If \code{NULL}, splines will be fitted for all plotIds.
#' @param knots The number of knots to use when fitting the spline.
#' @param useTimeNumber Should the timeNumber be used instead of the timePoint?
#' @param timeNumber If \code{useTimeNumber = TRUE}, a character vector
#' indicating the column containing the numerical time to use.
#' @param minNoTP The minimum number of time points for which data should be
#' available for a plant. Defaults to 80% of all time points present in the
#' TP object. No splines are fitted for plants with less than the minimum number
#' of timepoints.
#'
#' @return An object of class \code{HTPSpline}, a list with two
#' \code{data.frames}, \code{predDat} with predicted values and \code{coefDat}
#' with P-Spline coefficients on a dense grid.
#'
#' @examples
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## Fit P-Splines on a subset of genotypes
#' subGeno <- c("G070", "G160")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGeno,
#'                         knots = 50)
#'
#' ## Extract the data.frames with predicted values and P-Spline coefficients.
#' predDat <- fit.spline$predDat
#' head(predDat)
#'
#' coefDat <- fit.spline$coefDat
#' head(coefDat)
#'
#' ## Visualize the P-Spline predictions for one genotype.
#' plot(fit.spline, genotypes = "G160")
#'
#' ## Visualize the P-Spline predictions and first derivatives for one plant.
#' plot(fit.spline, plotIds = "c10r29", plotType =  "predictions")
#' plot(fit.spline, plotIds = "c10r29", plotType =  "derivatives")
#'
#' @family functions for fitting splines
#'
#' @export
fitSpline <- function(inDat,
                      trait,
                      genotypes = NULL,
                      plotIds = NULL,
                      knots = 50,
                      useTimeNumber = FALSE,
                      timeNumber = NULL,
                      minNoTP = NULL) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(inDat, "data.frame")) {
    stop("inDat should be a data.frame.\n")
  }
  if (isTRUE(useTimeNumber) &&
      (is.null(timeNumber) || !is.character(timeNumber) ||
       length(timeNumber) > 1)) {
    stop("timeNumber should be a character string of length 1.\n")
  }
  fitLevel <- if (hasName(x = inDat, name = "plotId")) "plotId" else
    "genotype"
  corrCols <- c("genotype", trait, if (fitLevel == "plotId") "plotId",
                if (useTimeNumber) timeNumber else "timePoint")
  if (!all(hasName(x = inDat, name = corrCols))) {
    stop("inDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!inherits(inDat[[fitLevel]], "factor")) {
    stop(fitLevel, " should be a factor column in inDat.\n")
  }
  if (!is.null(genotypes) &&
      (!is.character(genotypes) ||
       !all(genotypes %in% inDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in inDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) || !all(plotIds %in% inDat[["plotId"]]))) {
    stop("plotIds should be a character vector of plotIds in inDat.\n")
  }
  if (!is.numeric(knots) || length(knots) > 1 || knots < 0) {
    stop("knots should be a positive numerical value.\n")
  }
  if (knots < 4) {
    stop("Number of knots should be at least 4 for proper spline fitting.\n")
  }
  if (!is.null(minNoTP) && (!is.numeric(minNoTP) || length(minNoTP) > 1)) {
    stop("minNoTP should be a numerical value.\n")
  }
  if (!useTimeNumber) {
    if (!inherits(inDat[["timePoint"]], "POSIXct")) {
      stop("Column timePoint should be of class POSIXct.\n")
    }
    ## Convert time point to time number with the first time point as 0.
    minTime <- min(inDat[["timePoint"]], na.rm = TRUE)
    inDat[["timeNumber"]] <- as.numeric(inDat[["timePoint"]] - minTime)
  } else {
    if (!is.numeric(inDat[[timeNumber]])) {
      stop("timeNumber should be a numerical column.\n")
    }
    inDat[["timeNumber"]] <- inDat[[timeNumber]]
  }
  ## Restrict inDat to selected genotypes and plotIds.
  if (!is.null(genotypes)) {
    inDat <- inDat[inDat[["genotype"]] %in% genotypes, ]
  }
  if (!is.null(plotIds)) {
    inDat <- inDat[inDat[["plotId"]] %in% plotIds, ]
  }
  if (nrow(inDat) == 0) {
    stop("At least one valid combination of genotype and plotId should be ",
         "selected.\n")
  }
  ## Check if geno.decomp in present in inDat.
  useGenoDecomp <- hasName(x = inDat, name = "geno.decomp")
  if (fitLevel == "genotype") {
    ## If fitting at a genotype level set plotId to genotype.
    ## This way all further code can remain intact.
    if (useGenoDecomp) {
      inDat[["plotId"]] <- interaction(inDat[["genotype"]],
                                       inDat[["geno.decomp"]], drop = TRUE)
    } else {
      inDat[["plotId"]] <- inDat[["genotype"]]
    }
  }
  inDat <- droplevels(inDat)
  ## Create data.frame with plants and genotypes for adding genotype to results.
  if (useGenoDecomp) {
    plantGeno <- unique(inDat[c("plotId", "genotype", "geno.decomp")])
  } else {
    plantGeno <- unique(inDat[c("plotId", "genotype")])
  }
  ## Determine minimum number of time points required.
  nTimeNumber <- length(unique(inDat[["timeNumber"]]))
  if (is.null(minNoTP)) {
    minTP <- 0.8 * nTimeNumber
  } else if (minNoTP < 0 || minNoTP > nTimeNumber) {
    stop("minNoTP should be a number bewtween 0 and ", nTimeNumber, ".\n")
  } else {
    minTP <- minNoTP
  }
  ## Get number of non NA observations per plot for determining minimum
  ## number of knots.
  plotObs <- aggregate(x = inDat[[trait]], by = list(inDat[["plotId"]]),
                       FUN = function(plant) {
                         sum(!is.na(plant))
                       })
  ## Compute step size for prediction grid.
  ## Use smallest time gap between two points and divide that in 10 segments.
  minStep <- min(diff(sort(unique(inDat[["timeNumber"]]))))
  timeNumStep <- minStep / 9
  ## Create time range for prediction.
  ## Get range for time number and time point from data.
  timeNumRange <- range(inDat[["timeNumber"]])
  ## Create data.frame with time number and, if present,
  ## time point on prediction scale.
  timeRange <- data.frame(timeNumber = seq(from = timeNumRange[1],
                                           to = timeNumRange[2],
                                           by = timeNumStep))
  if (hasName(x = inDat, name = "timePoint")) {
    timePointRange <- range(inDat[["timePoint"]])
    timeRange[["timePoint"]] <- seq(from = timePointRange[1],
                                    to = timePointRange[2],
                                    by = timeNumStep * diff(timePointRange) /
                                      diff(timeNumRange))
  }
  ## Check for plotIds that have a limited amount of observations.
  plotTab <- table(inDat[!is.na(inDat[[trait]]), "plotId"])
  plotLimObs <- names(plotTab[plotTab < minTP])
  if (length(plotLimObs) > 5) {
    warning("More than 5 ", fitLevel, "s have observations for less than the ",
            "minimum number of time points, which is ", round(minTP), ". The  ",
            "first 5 are printed, to see them all run attr(..., 'plotLimObs') ",
            "on the output\n",
            paste(plotLimObs[1:5], collapse = ", "), "\n", call. = FALSE)
  } else if (length(plotLimObs) > 0) {
    warning("The following ", fitLevel, "s have observations for less than ",
            "the minimum number of time points, which is ", round(minTP), ":\n",
            paste(plotLimObs, collapse = ", "), "\n", call. = FALSE)
  }
  ## Fit splines.
  fitSp <- lapply(X = levels(plantGeno[["plotId"]]), FUN = function(plant) {
    ## Restrict data to current plant.
    dat <- inDat[inDat[["plotId"]] == plant, c("timeNumber", trait)]
    ## Manually select minimum number of time points.
    if (length(unique(dat[!is.na(dat[[trait]]), "timeNumber"])) >= minTP) {
      ## Manually set the number of knots.
      xmin <- min(dat[!is.na(dat[[trait]]), "timeNumber"]) - 1e-10
      xmax <- max(dat[!is.na(dat[[trait]]), "timeNumber"]) + 1e-10
      ## Construct vector of knots.
      knotsVec <- PsplinesKnots(xmin = xmin, xmax = xmax, degree = 3,
                                nseg = knots)
      ## Fit the P-spline using PsplinesREML.
      obj <- PsplinesREML(x = dat[["timeNumber"]], y = dat[[trait]],
                          knots = knotsVec)
      ## Extract the spline coefficients.
      coeff <- data.frame(obj.coefficients = obj$splineCoeffs, plotId = plant)
      coeff[["type"]] <- paste0("timeNumber", seq_len(nrow(coeff)))
      ## Restrict dense grid to points within observation range.
      timeRangePl <- timeRange[timeRange[["timeNumber"]] >=
                                 min(dat[!is.na(dat[[trait]]), "timeNumber"]) &
                                 timeRange[["timeNumber"]] <=
                                 max(dat[!is.na(dat[[trait]]), "timeNumber"]),
                               , drop = FALSE]
      ## Predictions on a dense grid.
      yPred <- predict(obj, x = timeRangePl$timeNumber)
      yDeriv <- predict(obj, x = timeRangePl$timeNumber, deriv = TRUE)
      yDeriv2 <-  predict(obj, x = timeRangePl$timeNumber, deriv2 = TRUE)
      ## Merge time, predictions and plotId.
      predDat <- data.frame(timeRangePl, pred.value = yPred,
                            deriv = yDeriv, deriv2 = yDeriv2, plotId = plant)
      return(list(coeff = coeff, predDat = predDat))
    } else {
      return(list(coeff = NULL, predDat = NULL))
    }
  })
  ## Bind all coefficients into one data.frame.
  coefTot <- do.call(rbind, lapply(fitSp, `[[`, 1))
  ## Add genotype and optionally geno.decomp.
  addCols <- colnames(plantGeno)[colnames(plantGeno) != "plotId"]
  coefTot[addCols] <- plantGeno[match(coefTot[["plotId"]],
                                      plantGeno[["plotId"]]), addCols]
  ## Bind all predictions into one data.frame.
  predTot <- do.call(rbind, lapply(fitSp, `[[`, 2))
  ## Add genotype.
  predTot[addCols] <- plantGeno[match(predTot[["plotId"]],
                                      plantGeno[["plotId"]]), addCols]
  if (fitLevel == "genotype") {
    ## Remove plotId (duplicated genotype) from output.
    coefTot[["plotId"]] <- NULL
    predTot[["plotId"]] <- NULL
  }
  ## Create output.
  res <- structure(list(coefDat = coefTot, predDat = predTot),
                   modDat = inDat,
                   trait = trait,
                   useTimeNumber = useTimeNumber,
                   fitLevel = fitLevel,
                   useGenoDecomp = useGenoDecomp,
                   plotLimObs = plotLimObs,
                   class = c("HTPSpline", "list"))
  return(res)
}

#' Plot the results of a fitted spline.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class \code{HTPSpline}.
#' @param ... Ignored.
#' @param plotType A character string indicating which spline component
#' should be plotted, either predictions, derivatives or second derivatives
#' ("derivatives2").
#' @param genotypes A character vector indicating the genotypes for which
#' spline components should be plotted.
#' @param plotIds A character vector indicating the plotIds for which spline
#' components should be plotted.
#'
#' @return A list of object of class ggplot is invisibly returned.
#'
#' @examples
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## Fit P-Splines on a subset of genotypes
#' subGeno <- c("G070", "G160")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGeno,
#'                         knots = 50)
#'
#' ## Visualize the P-Spline predictions for one genotype.
#' plot(fit.spline, genotypes = "G160")
#'
#' ## Visualize the first and second derivatives of the predictions for one plant.
#' plot(fit.spline, plotIds = "c10r29", plotType =  "derivatives")
#' plot(fit.spline, plotIds = "c10r29", plotType =  "derivatives2")
#'
#' @family functions for fitting splines
#'
#' @export
plot.HTPSpline <- function(x,
                           ...,
                           plotType = c("predictions", "derivatives",
                                        "derivatives2"),
                           genotypes = NULL,
                           plotIds = NULL,
                           title = NULL,
                           output = TRUE,
                           outFile = NULL,
                           outFileOpts = NULL) {
  plotType <- match.arg(plotType)
  plotVar <- if (plotType == "predictions") "pred.value" else if
  (plotType == "derivatives") "deriv" else "deriv2"
  modDat <- attr(x, which = "modDat")
  trait <- attr(x, which = "trait")
  fitLevel <- attr(x, which = "fitLevel")
  useGenoDecomp <- attr(x, which = "useGenoDecomp")
  useTimeNumber <- attr(x, which = "useTimeNumber")
  predDat <- x$predDat
  if (!is.null(genotypes) &&
      (!is.character(genotypes) ||
       !all(genotypes %in% predDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in predDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) || !all(plotIds %in% predDat[["plotId"]]))) {
    stop("plotIds should be a character vector of plotIds in predDat.\n")
  }
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    on.exit(dev.off(), add = TRUE)
    do.call(pdf, args = outFileOpts)
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
  if (nrow(predDat) > 0) {
    if (fitLevel == "genotype" && useGenoDecomp) {
      modDat[["fitLevInt"]] <-
        interaction(modDat[["genotype"]], modDat[["geno.decomp"]], drop = TRUE)
      predDat[["fitLevInt"]] <-
        interaction(predDat[["genotype"]], predDat[["geno.decomp"]], drop = TRUE)
      fitLevel <- "fitLevInt"
      fitLevels <- levels(modDat[["fitLevInt"]])
      plotLevel <- c("genotype", "geno.decomp")
    } else {
      fitLevels <- unique(modDat[[fitLevel]])
      plotLevel <- fitLevel
    }
    allNA <- sapply(X = fitLevels, FUN = function(x) {
      all(is.na(modDat[modDat[[fitLevel]] == x, trait]))
    })
    modDat <- modDat[!modDat[[fitLevel]] %in% fitLevels[allNA], ]
    predDat <- predDat[!predDat[[fitLevel]] %in% fitLevels[allNA], ]
  }
  modDat <- droplevels(modDat)
  predDat <- droplevels(predDat)
  if (nrow(predDat) == 0) {
    stop("At least one valid combination of genotype and plotId should be ",
         "selected.\n")
  }
  ## Construct plot title.
  if (is.null(title)) {
    if (plotType == "predictions") {
      title <- "Corrected data and P-spline prediction"
    } else if (plotType == "derivatives") {
      title <- "P-spline first derivatives"
    } else {
      title <- "P-spline second derivatives"
    }
  }
  timeVar <- if (useTimeNumber) "timeNumber" else "timePoint"
  p <- ggplot2::ggplot(modDat, ggplot2::aes_string(x = timeVar, y = trait)) +
    ggplot2::geom_line(data = predDat,
                       ggplot2::aes_string(x = timeVar, y = plotVar),
                       col = "blue", na.rm = TRUE) +
    ggplot2::labs(title = title, y = trait, x = timeVar)
  if (plotType == "predictions") {
    p <- p + ggplot2::geom_point(na.rm = TRUE)
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
    rowPag <- min(ceiling(plotsPag / 5), 5)
    colPag <- ifelse(plotsPag >= 5, 5, plotsPag)
  }
  ## Build pages of plots.
  pPag <- vector(mode = "list", length = nPag)
  for (i in 1:nPag) {
    pPag[[i]] <- p +
      ggforce::facet_wrap_paginate(facets = plotLevel, nrow = rowPag[i],
                                   ncol = colPag[i],
                                   labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
                                   page = i)
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
#' @param estimate The P-Spline component for which the estimate should be
#' extracted, the predictions, the first derivatives or the second derivatives
#' ("derivatives2")
#' @param what The types of estimate that should be extracted. Either minimum
#' ("min"), maximum ("max"), mean, area under the curve ("AUC") or a percentile.
#' Percentiles should be given as p + percentile. E.g. for the 10th percentile
#' specify what = "p10". Multiple types of estimate can be extracted at once.
#' @param AUCScale The area under the curve is dependent on the scale used on
#' the x-axis. By default the area is computed assuming a scale in minutes. This
#' can be changed to either hours or days.
#' @param timeMin The lower bound of the time interval from which the
#' estimates should be extracted. If \code{NULL} the smallest time value for
#' which the splines were fitted is used.
#' @param timeMax The upper bound of the time interval from which the
#' estimates should be extracted. If \code{NULL} the largest time value for
#' which the splines were fitted is used.
#' @param genotypes A character vector indicating the genotypes for which
#' estimates should be extracted. If \code{NULL}, estimates will be extracted
#' for all genotypes for which splines where fitted.
#' @param plotIds A character vector indicating the plotIds for which
#' estimates should be extracted. If \code{NULL}, estimates will be extracted
#' for all plotIds for which splines where fitted.
#'
#' @return A data.frame containing the estimated parameters.
#'
#' @examples
#' ## Run the function to fit P-splines on a subset of genotypes.
#' subGeno <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGeno,
#'                         knots = 50)
#'
#' ## Estimate the maximum value of the predictions at the beginning of the time course.
#' paramVator <- estimateSplineParameters(HTPSpline = fit.spline,
#'                                        estimate = "predictions",
#'                                        what = "max",
#'                                        timeMin = 1527784620,
#'                                        timeMax = 1528500000,
#'                                        genotypes = subGeno)
#' head(paramVator)
#'
#' ## Estimate the minimum and maximum value of the predictions.
#' paramVator2 <- estimateSplineParameters(HTPSpline = fit.spline,
#'                                         estimate = "predictions",
#'                                         what = c("min", "max"),
#'                                         genotypes = subGeno)
#' head(paramVator2)
#'
#' @family functions for fitting splines
#'
#' @export
estimateSplineParameters <- function(HTPSpline,
                                     estimate = c("predictions", "derivatives",
                                                  "derivatives2"),
                                     what = c("min", "max", "mean", "AUC", "p"),
                                     AUCScale = c("min", "hour", "day"),
                                     timeMin = NULL,
                                     timeMax = NULL,
                                     genotypes = NULL,
                                     plotIds = NULL) {
  estimate <- match.arg(estimate)
  estVar <- if (estimate == "predictions") "pred.value" else if
  (estimate == "derivatives") "deriv" else "deriv2"
  useTimeNumber <- attr(HTPSpline, which = "useTimeNumber")
  useGenoDecomp <- attr(HTPSpline, which = "useGenoDecomp")
  fitLevel <- attr(HTPSpline, which = "fitLevel")
  predDat <- HTPSpline$predDat
  if (!is.null(genotypes) &&
      (!is.character(genotypes) ||
       !all(genotypes %in% predDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in predDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) || !all(plotIds %in% predDat[["plotId"]]))) {
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
    if (timeMin < min(predDat[[timeVar]]) ||
        timeMin > max(predDat[[timeVar]])) {
      stop("timeMin should be within the time interval in the data.\n")
    }
  }
  if (is.null(timeMax)) {
    timeMax <- max(predDat[[timeVar]])
  } else {
    if (timeMax < min(predDat[[timeVar]]) ||
        timeMax > max(predDat[[timeVar]])) {
      stop("timeMax should be within the time interval in the data.\n")
    }
  }
  if (timeMin >= timeMax) {
    stop("timeMax should be larger than timeMin.\n")
  }
  ## Restrict predDat to time interval.
  predDat <- predDat[predDat[[timeVar]] >= timeMin &
                       predDat[[timeVar]] <= timeMax, ]
  ## Loop over what
  whatRes <- lapply(X = what, FUN = function(w) {
    if (substr(w, 1, 1) == "p") {
      percentile <- suppressWarnings(as.numeric(substring(what, 2))) / 100
      if (is.na(percentile) || percentile < 0 || percentile > 1) {
        stop("A percentile should be give as pN, with N between 0 and 100.\n")
      }
      ## Percentiles are calculated using quantile
      estFun <- "quantile"
    } else {
      estFun <- match.arg(arg = w, choices = c("min", "max", "mean", "AUC"))
    }
    ## Area under the curve corresponds to sum.
    if (w == "AUC") {
      intWidth <- diff(predDat[1:2, timeVar])
      if (timeVar == "timePoint") {
        ## x-axis scale for time variables as computed by diff is in minutes.
        ## For conversino to hours/days divide by appropriate factor.
        AUCScale <- match.arg(AUCScale)
        if (AUCScale == "hour") {
          intWidth <- intWidth / 60
        } else if (AUCScale == "day") {
          intWidth <- intWidth / (24 * 60)
        }
      }
      estFun <- function(x, ...) {
        ## All intervals have the same (small) width.
        ## Just summing and multiplying by this width gives a good
        ## approximation of the area under the curve.
        return(as.numeric(sum(x) * intWidth))
      }
    }
    ## Get estimates.
    resW <- aggregate(x = predDat[[estVar]],
                      by = predDat[c("genotype", if (useGenoDecomp) "geno.decomp",
                                     if (fitLevel == "plotId") "plotId")],
                      FUN = estFun,
                      probs = if (is.character(estFun) &&
                                  estFun == "quantile") percentile)
    colnames(resW)[colnames(resW) == "x"] <- paste0(w, "_", estimate)
    ## For min and max get corresponding time point.
    if (w %in% c("min", "max")) {
      resW <- merge(resW, predDat, by.x = colnames(resW),
                    by.y = c(colnames(resW)[-ncol(resW)], estVar))
      resW <- resW[, 1:(4 + (fitLevel == "plotId"))]
      colnames(resW)[colnames(resW) %in% c("timeNumber", "timePoint")] <-
        paste0(w, "_", colnames(resW)[colnames(resW) %in% c("timeNumber", "timePoint")])
    }
    return(resW)
  })
  if (length(whatRes) > 1) {
    res <- do.call(what = merge, args = whatRes)
  } else {
    res <- whatRes[[1]]
  }
  res <- structure(res,
                   what = what,
                   class = c("HTPSplineEst", "data.frame"))
  return(res)
}

#### Helper function for fitting splines.

#' Spectral decomposition of D'D
#'
#' Spectral decomposition of D'D, returns a q x (q - ord) matrix.
#' @noRd
#' @keywords internal
calcUsc <- function(q,
                    ord) {
  D <- diff(diag(q), diff = ord)
  DtD <- crossprod(D)
  decomp <- eigen(DtD, symmetric = TRUE)
  U <- decomp$vectors[, 1:(q - ord)]
  d <- decomp$values[1:(q - ord)]
  return(U %*% diag(1 / sqrt(d)))
}

#' Equally placed knots.
#'
#' @noRd
#' @keywords internal
PsplinesKnots <- function(xmin,
                          xmax,
                          degree,
                          nseg) {
  dx <- (xmax - xmin) / nseg
  knots <- seq(xmin - degree * dx, xmax + degree * dx, by = dx)
  attr(knots, "degree") <- degree
  return(knots)
}

#' Fit P-Spline
#'
#' @noRd
#' @keywords internal
PsplinesREML <- function(x,
                         y,
                         knots,
                         lambda = 1,
                         optimize = TRUE) {
  ## Remove missing values in y from x and y.
  x <- x[!is.na(y)]
  y <- y[!is.na(y)]
  ## Set defaults.
  pord <- 2
  degree <- 3
  ## Construct B-Spline base.
  B <- splines::splineDesign(knots = knots, x = x, derivs = rep(0, length(x)),
                             ord = degree + 1)
  q <- ncol(B)
  Usc <- calcUsc(q = q, ord = pord)
  ## Calculate the linear/fixed parts.
  UX1 <- cbind(1, scale(1:q))
  UX1[, 2] <- UX1[, 2] / norm(UX1[, 2], type = "2")
  X <- B %*% UX1
  Z <- B %*% Usc
  U <- cbind(X, Z)
  UtU <- crossprod(U)
  ## max ED for random part.
  maxED <- length(unique(x)) - pord
  UtY <- t(U) %*% y
  P <- diag(c(0, 0, rep(1, q - 2)))
  phi <- 1
  psi <- lambda
  n <- length(x)
  p <- 2
  for (it in 1:100) {
    C <- phi * UtU + psi * P
    ## calculate effective dimensions.
    Cinv <- solve(C)
    EDfSpline <- p
    EDrSpline <- ncol(Z) - psi * sum(diag(Cinv %*% P))
    EDRes <- n - p - EDrSpline
    ## calculate coefficients and residuals.
    coeffs <- phi * Cinv %*% UtY
    resid <- y - U %*% coeffs
    if (optimize) {
      ## Optimization steps.
      ## Compute new values for psi and phi
      psiNew <- EDrSpline / (sum(coeffs * (P %*% coeffs)) + 1.0e-20)
      phiNew <- EDRes / (sum(resid ^ 2) + 1.0e-20)
      pOld <- log(c(phi, psi))
      pNew <- log(c(phiNew, psiNew))
      diff <- max(abs(pOld - pNew))
      phi <- phiNew
      psi <- psiNew
      if (diff < 1.0e-12) break
    }
  }
  ## Calculate spline coefficients.
  U <- cbind(UX1, Usc)
  splineCoeffs <- U %*% coeffs
  ## Construct output.
  res <- structure(
    list(maxED = maxED, ED = EDrSpline + p, coeffs = coeffs,
         splineCoeffs = splineCoeffs, knots = knots, nObs = n,
         x = x, y = y, optimize = optimize, UX1 = UX1, Usc = Usc),
    class = c("PsplinesREML", "list"))
  return(res)
}

#' Prediction from P-Spline.
#'
#' @noRd
#' @keywords internal
predict.PsplinesREML <- function(obj,
                                 x,
                                 deriv = FALSE,
                                 deriv2 = FALSE) {
  d <- if (deriv) 1 else if (deriv2) 2 else 0
  Bgrid = splines::splineDesign(obj$knots, x, derivs = d, ord = 4)
  U <- cbind(obj$UX1, obj$Usc)
  pred <- as.vector(Bgrid %*% U %*% obj$coeffs)
  return(pred)
}

