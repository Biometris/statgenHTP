#' Extract estimates from fitted splines.
#'
#' Function for extracting parameter estimates from fitted splines on a
#' specified interval.
#'
#' @param x An object of class HTPSpline, the output of the
#' \code{\link{fitSpline}} function, or class splineHDm, the output of the
#' \code{\link{fitSplineHDM}} function
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
#' which the splines were fitted is used. \code{timeMin} should be given as a
#' numerical value that corresponds to the time scale used for fitting the
#' splines. See the examples.
#' @param timeMax The upper bound of the time interval from which the
#' estimates should be extracted. If \code{NULL} the largest time value for
#' which the splines were fitted is used. \code{timeMin} should be given as a
#' numerical value that corresponds to the time scale used for fitting the
#' splines. See the examples.
#' @param genotypes A character vector indicating the genotypes for which
#' estimates should be extracted. If \code{NULL}, estimates will be extracted
#' for all genotypes for which splines where fitted.
#' @param plotIds A character vector indicating the plotIds for which
#' estimates should be extracted. If \code{NULL}, estimates will be extracted
#' for all plotIds for which splines where fitted.
#' @param fitLevel A character string indicating at which level of the data
#' the parameter estimates should be made. Only used for splines fitted using
#' \code{\link{fitSplineHDM}}.
#'
#' @returns An object of class splineEst, a data.frame containing the
#' estimated parameters.
#'
#' @examples
#' ### Estimate parameters for fitted P-splines.
#'
#' ## Run the function to fit P-splines on a subset of genotypes.
#' subGeno <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGeno,
#'                         knots = 50)
#'
#' ## Estimate the maximum value of the predictions at the beginning of the time course.
#' ## The spline was fitted at a timePoints scale, i.e. date-time so
#' ## timeMin and timeMax should be given at this scale as well.
#' paramVator <- estimateSplineParameters(x = fit.spline,
#'                                        estimate = "predictions",
#'                                        what = "max",
#'                                        timeMin = 1527784620,
#'                                        timeMax = 1528500000,
#'                                        genotypes = subGeno)
#' head(paramVator)
#'
#' ## Create a boxplot of the estimates.
#' plot(paramVator, plotType = "box")
#'
#' ## Estimate the minimum and maximum value of the predictions.
#' paramVator2 <- estimateSplineParameters(x = fit.spline,
#'                                         estimate = "predictions",
#'                                         what = c("min", "max"),
#'                                         genotypes = subGeno)
#' head(paramVator2)
#'
#'
#' ### Estimate parameters for fitted HDM-splines.
#'
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## We need to specify the genotype-by-treatment interaction.
#' ## Treatment: water regime (WW, WD).
#' spatCorrectedArch[["treat"]] <- substr(spatCorrectedArch[["geno.decomp"]],
#'                                       start = 1, stop = 2)
#' spatCorrectedArch[["genoTreat"]] <-
#'   interaction(spatCorrectedArch[["genotype"]],
#'              spatCorrectedArch[["treat"]], sep = "_")
#'
#' ## Fit P-Splines Hierarchical Curve Data Model for selection of genotypes.
#' fit.psHDM  <- fitSplineHDM(inDat = spatCorrectedArch,
#'                           trait = "LeafArea_corr",
#'                           genotypes = c("GenoA14_WD", "GenoA51_WD",
#'                                        "GenoB11_WW", "GenoB02_WD",
#'                                        "GenoB02_WW"),
#'                           time = "timeNumber",
#'                           pop = "geno.decomp",
#'                           genotype = "genoTreat",
#'                           plotId = "plotId",
#'                           difVar = list(geno = FALSE, plot = FALSE),
#'                           smoothPop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothGeno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothPlot = list(nseg = 4, bdeg = 3, pord = 2),
#'                           weights = "wt",
#'                           trace = FALSE)
#'
#' ## Estimate minimum, maximum, and mean for predictions at the genotype level.
#' ## The spline was fitted at the timeNumber scale, so timeMax
#' ## should be given at that scale as well.
#' paramArch <- estimateSplineParameters(x = fit.psHDM,
#'                                      what = c("min", "max", "mean"),
#'                                      fitLevel = "geno",
#'                                      estimate = "predictions",
#'                                      timeMax = 28)
#' head(paramArch)
#'
#' ## Create a boxplot of the estimates.
#' plot(paramArch, plotType = "box")
#'
#' ## Estimate area under the curve for predictions at the plot level.
#' paramArch2 <- estimateSplineParameters(x = fit.psHDM,
#'                                       what = "AUC",
#'                                       fitLevel = "plot",
#'                                       estimate = "predictions")
#' head(paramArch2)
#'
#' @family functions for spline parameter estimation
#'
#' @export
estimateSplineParameters <- function(x,
                                     estimate = c("predictions", "derivatives",
                                                  "derivatives2"),
                                     what = c("min", "max", "mean", "AUC", "p"),
                                     AUCScale = c("min", "hour", "day"),
                                     timeMin = NULL,
                                     timeMax = NULL,
                                     genotypes = NULL,
                                     plotIds = NULL,
                                     fitLevel = c("geno", "plot", "genoDev",
                                                  "plotDev")) {
  ## General settings.
  estimate <- match.arg(estimate)
  AUCScale <- match.arg(AUCScale)
  trait <- attr(x, which = "trait")
  if (inherits(x, "HTPSpline")) {
    ## HTP spline specific settings.
    estVar <- if (estimate == "predictions") "pred.value" else if
    (estimate == "derivatives") "deriv" else "deriv2"
    useTimeNumber <- attr(x, which = "useTimeNumber")
    useGenoDecomp <- attr(x, which = "useGenoDecomp")
    fitLevel <- attr(x, which = "fitLevel")
    predDat <- x$predDat
    ## Construct levels for aggregating.
    aggLevs <- unique(c("genotype", if (useGenoDecomp) "geno.decomp",
                        fitLevel))
  } else if (inherits(x, "psHDM")) {
    ## HDM spline specific settings.
    fitLevel <- match.arg(fitLevel)
    ## To title case doesn't work for genoDev and plotDev.
    ## It leaves them as is.
    estVarBase <- paste0("f", toupper(substring(fitLevel, first = 1, last = 1)),
                         substring(fitLevel, first = 2))
    estVar <- if (estimate == "predictions") estVarBase else if
    (estimate == "derivatives") paste0(estVarBase, "Deriv1") else
      paste0(estVarBase, "Deriv2")
    useTimeNumber <- TRUE
    useGenoDecomp <- FALSE
    ## Check that predictions were made at the level estimates are computed.
    fitLevelName <- paste0(substr(fitLevel, 1, 4), "Level")
    if (is.null(x[[fitLevelName]])) {
      stop("No predictions were made at ", substr(fitLevel, 1, 4), "level.\n")
    }
    predDat <- x[[fitLevelName]]
    ## Construct levels for aggregating.
    aggLevs <- unique(c("pop", "genotype",
                        if (fitLevel %in% c("plot", "plotDev")) "plotId"))
  }
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
  ## Construct data for aggregating.
  aggDat <- predDat[, aggLevs, drop = FALSE]
  ## Loop over what
  whatRes <- lapply(X = what, FUN = function(w) {
    if (substr(w, 1, 1) == "p") {
      percentile <- suppressWarnings(as.numeric(substring(w, 2))) / 100
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
        ## For converting to hours/days divide by appropriate factor.
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
                      by = aggDat,
                      FUN = estFun,
                      probs = if (is.character(estFun) &&
                                  estFun == "quantile") percentile)
    colnames(resW)[colnames(resW) == "x"] <- paste0(w, "_", estimate)
    ## For min and max get corresponding time point.
    if (w %in% c("min", "max")) {
      resW <- merge(resW, predDat, by.x = colnames(resW),
                    by.y = c(colnames(resW)[-ncol(resW)], estVar))
      nTimeCols <- sum(c("timeNumber", "timePoint") %in% colnames(predDat))
      resW <- resW[, 1:(1 + nTimeCols + length(aggLevs))]
      colnames(resW)[colnames(resW) %in% c("timeNumber", "timePoint")] <-
        paste0(w, "_", colnames(resW)[colnames(resW) %in% c("timeNumber", "timePoint")])
    }
    return(resW)
  })
  if (length(whatRes) > 1) {
    res <- Reduce(f = merge, x = whatRes)
  } else {
    res <- whatRes[[1]]
  }
  res <- structure(res,
                   what = what,
                   estimate = estimate,
                   trait = trait,
                   useGenoDecomp = useGenoDecomp,
                   aggLevs = aggLevs,
                   class = c("splineEst", "data.frame"))
  return(res)
}

#' Plot the results of estimated spline parameters.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class \code{splineEst}
#' @param ... Ignored.
#' @param plotType A character string indicating the type of plot to be made.
#' @param what The types of estimate that should be plotted.
#'
#' @returns A list of objects of class ggplot is invisibly returned.
#'
#' @family functions for spline parameter estimation
#'
#' @export
plot.splineEst <- function(x,
                           ...,
                           plotType = c("box", "hist"),
                           what = attr(x, "what"),
                           title = NULL,
                           output = TRUE,
                           outFile = NULL,
                           outFileOpts = NULL) {
  plotType <- match.arg(plotType)
  what <- match.arg(what, several.ok = TRUE)
  trait <- attr(x, which = "trait")
  estimate <- attr(x, which = "estimate")
  useGenoDecomp <- attr(x, which = "useGenoDecomp")
  aggLevs <- attr(x, which = "aggLevs")
  plotVar <- if (length(aggLevs) == 1) aggLevs else aggLevs[length(aggLevs) - 1]
  for (aggLev in aggLevs) {
    if (!is.factor(x[[aggLev]])) {
      x[[aggLev]] <- as.factor(x[[aggLev]])
    }
  }
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    on.exit(dev.off(), add = TRUE)
    do.call(pdf, args = outFileOpts)
  }
  if (plotType == "box") {
    pTot <- lapply(X = what, FUN = function(w) {
      estVar <- paste0(w, "_", estimate)
      pWhat <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[plotVar]],
                                               y = .data[[estVar]])) +
        ggplot2::geom_boxplot(na.rm = TRUE) +
        ggplot2::labs(title = title, y = (paste(w, "of", trait))) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5))
      if (useGenoDecomp) {
        pWhat <- pWhat +
          ggplot2::facet_grid(. ~ geno.decomp, scales = "free_x",
                              space = "free_x")
      }
      return(pWhat)
    })
  } else if (plotType == "hist") {
    pTot <- lapply(X = what, FUN = function(w) {
      estVar <- paste0(w, "_", estimate)
      pWhat <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[estVar]])) +
        ggplot2::geom_histogram(na.rm = TRUE, bins = 10) +
        ggplot2::labs(title = title, y = (paste(w, "of", trait)))
      if (useGenoDecomp) {
        pWhat <- pWhat +
          ggplot2::facet_grid(. ~ geno.decomp, scales = "free_x",
                              space = "free_x")
      }
      return(pWhat)
    })
  }
  if (output) {
    for (p in pTot) {
      plot(p)
    }
  }
  invisible(pTot)
}

