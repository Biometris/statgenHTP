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
#' @return An object of class HTPSplineEst, a data.frame containing the
#' estimated parameters.
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
#' @family functions for spline parameter estimation
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
  AUCScale <- match.arg(AUCScale)
  estVar <- if (estimate == "predictions") "pred.value" else if
  (estimate == "derivatives") "deriv" else "deriv2"
  useTimeNumber <- attr(HTPSpline, which = "useTimeNumber")
  useGenoDecomp <- attr(HTPSpline, which = "useGenoDecomp")
  fitLevel <- attr(HTPSpline, which = "fitLevel")
  trait <- attr(HTPSpline, which = "trait")
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
                   estimate = estimate,
                   trait = trait,
                   useGenoDecomp = useGenoDecomp,
                   class = c("HTPSplineEst", "data.frame"))
  return(res)
}

#' Plot the results of estimated spline parameters.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class \code{HTPSpline}.
#' @param ... Ignored.
#' @param plotType A character string indicating the type of plot to be made.
#' @param what The types of estimate that should be extracted.
#'
#' @return A list of object of class ggplot is invisibly returned.
#'
#' @family functions for spline parameter estimation
#'
#' @export
plot.HTPSplineEst <- function(x,
                              ...,
                              plotType = c("box", "hist"),
                              what = attr(x, "what"),
                              title = NULL,
                              output = TRUE,
                              outFile = NULL,
                              outFileOpts = NULL) {
  plotType <- match.arg(plotType)
  what <- match.arg(what)
  trait <- attr(x, which = "trait")
  estimate <- attr(x, which = "estimate")
  useGenoDecomp <- attr(x, which = "useGenoDecomp")
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
      pWhat <- ggplot2::ggplot(x,
                               ggplot2::aes_string(x = "genotype", y = estVar)) +
        ggplot2::geom_boxplot(na.rm = TRUE) +
        ggplot2::labs(title = title, y = (paste(w, "of", trait))) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5))
      if (useGenoDecomp) {
        pWhat <- pWhat + ggplot2::facet_grid(. ~ geno.decomp, scales = "free_x",
                                             space = "free_x")
      }
      if (output) {
        plot(pWhat)
      }
    })
  } else if (plotType == "hist") {
    pTot <- lapply(X = what, FUN = function(w) {
      estVar <- paste0(w, "_", estimate)
      pWhat <- ggplot2::ggplot(x, ggplot2::aes_string(x = estVar)) +
        ggplot2::geom_histogram(na.rm = TRUE, bins = 10) +
        ggplot2::labs(title = title, y = (paste(w, "of", trait)))
      if (useGenoDecomp) {
        pWhat <- pWhat + ggplot2::facet_grid(. ~ geno.decomp, scales = "free_x",
                                             space = "free_x")
      }
      if (output) {
        plot(pWhat)
      }
    })
  }
  invisible(pTot)
}

