#' Fit Splines
#'
#' Function for fitting splines.
#'
#' @param inDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param genotypes A character vector indicating the genotypes for which
#' splines are fitted. If \code{NULL}, splines will be fitted for all genotypes.
#' @param plotIds A character vector indicating the plotIds for which splines
#' are fitted. If \code{NULL}, splines will be fitted for all plotIds.
#' @param knots The number of knots to use when fitting the spline.
#' @param useTimeNumber Should the timeNumber be used instead of the timePoint?
#' @param timeNumber If \code{useTimeNumber = TRUE}, a character vector
#' indicating the column containing the numerical time to use.
#' @param perMinTP The percentage of minimum number of time points to use.
#'
#' @examples ## The data from the Phenovator platform have been corrected for
#' # spatial trends and time points outliers have been removed
#' # Format the timepoint
#' spatCorrVator$timePoint <- lubridate::as_datetime(spatCorrVator$timePoint)
#' # Run the function to fit P-spline using the mgcv package on a subset of
#' # genotypes
#' subGeno <- c("G70","G160","G151","G179","G175","G4","G55")
#' fit.spline <- fitSpline(inDat = spatCorrVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGeno,
#'                         knots = 50,
#'                         perMinTP = 0.8)
#' # Extracting the tables of predicted values and P-spline coeficients
#' pred.Dat <- fit.spline$predDat
#' coef.Dat <- fit.spline$coefDat
#'
#' ## We can then visualise the P-spline predictions and first derivatives
#' # for a subset of genotypes...
#' plot(fit.spline, genotypes = "G160")
#' # ...or for a subset of plots.
#' plot(fit.spline, plotIds = "c10r29", plotType =  "predictions")
#' plot(fit.spline, plotIds = "c10r29", plotType =  "derivatives")
#'
#' @family Fit splines
#'
#' @export
fitSpline <- function(inDat,
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
  if (!inherits(inDat, "data.frame")) {
    stop("inDat should be a data.frame.\n")
  }
  if (useTimeNumber && (is.null(timeNumber) || !is.character(timeNumber) ||
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
  if (!is.null(genotypes) &&
      (!is.character(genotypes) &&
       !all(genotypes %in% inDat[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes in inDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) && !all(plotIds %in% inDat[["genotype"]]))) {
    stop("plotIds should be a character vector of plotIds in inDat.\n")
  }
  if (!is.numeric(knots) || length(knots) > 1 || knots < 0) {
    stop("knots should be a positive numerical value.\n")
  }
  if (knots < 4) {
    stop("Number of knots should be at least 4 for proper spline fitting.\n")
  }
  if (!is.numeric(perMinTP) || length(perMinTP) > 1 || perMinTP < 0 ||
      perMinTP > 1) {
    stop("perMinTP should be a numerical value between 0 and 1.\n")
  }
  if (!useTimeNumber) {
    inDat[["timeNumber"]] <- as.numeric(inDat[["timePoint"]])
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
  minTP <- perMinTP * length(unique(inDat[["timeNumber"]]))
  ## Get number of non NA observations per plot for determining minimum
  ## number of knots.
  plotObs <- aggregate(x = inDat[[trait]], by = list(inDat[["plotId"]]),
                       FUN = function(plant) {
                         sum(!is.na(plant))
                       })
  ## Number of knots cannot be larger than the minimum number of observations
  ## determined over all plots that have at least the number of required
  ## observations.
  minKnots <- min(plotObs[plotObs[["x"]] > minTP, "x"])
  if (knots >= minKnots) {
    stop("The number of knots cannot be larger than the minimum number of ",
         "required time points, which is ", minKnots)
  }
  ## Construct formula for fitting model.
  modForm <- as.formula(paste(trait, "~s(timeNumber, bs = 'ps', k = ",
                              knots, ")"))
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
  ## Fit splines.
  fitSp <- lapply(X = levels(plantGeno[["plotId"]]), FUN = function(plant) {
    ## Restrict data to current plant.
    dat <- inDat[inDat[["plotId"]] == plant, c("timeNumber", trait)]
    ## Manually select minimum number of time points.
    if (length(unique(dat[!is.na(dat[[trait]]), "timeNumber"])) >= minTP) {
      ## Fit the P-spline using gam() function in mgcv.
      ## Manually set the number of knots.
      ## Depends on the number of time points and shape of the curve.
      obj <- mgcv::gam(modForm, data = dat, method = "REML")
      ## Extract the spline coefficients.
      coeff <- data.frame(obj$coefficients, plotId = plant)
      coeff$type <- row.names(coeff)
      ## Restrict dense grid to points within observation range.
      timeRangePl <- timeRange[timeRange[["timeNumber"]] >=
                                 min(dat[["timeNumber"]]) &
                                 timeRange[["timeNumber"]] <=
                                 max(dat[["timeNumber"]]),
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
#' @family Fit splines
#'
#' @export
plot.HTPSpline <- function(x,
                           ...,
                           plotType = c("predictions", "derivatives"),
                           genotypes = NULL,
                           plotIds = NULL,
                           title = NULL,
                           output = TRUE) {
  plotType <- match.arg(plotType)
  plotVar <- if (plotType == "predictions") "pred.value" else "deriv"
  modDat <- attr(x, which = "modDat")
  trait <- attr(x, which = "trait")
  fitLevel <- attr(x, which = "fitLevel")
  useGenoDecomp <- attr(x, which = "useGenoDecomp")
  useTimeNumber <- attr(x, which = "useTimeNumber")
  predDat <- x$predDat
  if (!is.null(genotypes) &&
      (!is.character(genotypes) &&
       !all(genotypes %in% predDat[["genotype"]]))) {
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
    } else {
      title <- "P-spline first derivatives"
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
    rowPag <- min(plotsPag %/% 5 + 1, 5)
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
#' @examples # Format the timepoint
#' spatCorrVator$timePoint <- lubridate::as_datetime(spatCorrVator$timePoint)
#' # Run the function to fit P-spline using the mgcv package on a subset of
#' # genotypes
#' subGeno <- c("G70","G160","G151","G179","G175","G4","G55")
#' fit.spline <- fitSpline(inDat = spatCorrVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGeno,
#'                         knots = 50,
#'                         perMinTP = 0.8)
#'
#' # Estimate the maximum value of the trait at the beginning of the time course
#' paramVator1 <-
#'      estimateSplineParameters(HTPSpline = fit.spline,
#'                               estimate = "predictions",
#'                               what = "max",
#'                               timeMin = 1527784620,
#'                               timeMax = 1528500000,
#'                               genotypes = subGeno)
#'
#' @family Fit splines
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
  useGenoDecomp <- attr(HTPSpline, which = "useGenoDecomp")
  fitLevel <- attr(HTPSpline, which = "fitLevel")
  predDat <- HTPSpline$predDat
  if (!is.null(genotypes) &&
      (!is.character(genotypes) &&
       !all(genotypes %in% predDat[["genotype"]]))) {
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
  ## Get estimates.
  res <- aggregate(x = predDat[[estVar]],
                   by = predDat[c("genotype", if (useGenoDecomp) "geno.decomp",
                                  if (fitLevel == "plotId") "plotId")],
                   FUN = what)
  return(res)
}
