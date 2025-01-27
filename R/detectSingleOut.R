#' Detect outliers for single observations
#'
#' Detect outlying observations in a time series by modeling each plotId using
#' a local regression.
#'
#' See locfit() help function from the locfit R library. The user can act on:
#' \describe{
#'   \item{nnLocfit}{the constant of the smoothing parameter. Increase nnLocfit
#'   to have a very smooth curve}
#'   \item{confIntSize}{the level to calculate the confidence interval. Increase
#'   confIntSize to exclude less outliers}
#' }
#'
#' @param TP An object of class \code{TP}.
#' @param trait A character vector indicating the trait to model in \code{TP}.
#' @param plotIds A character vector of plotIds for which the outliers should be
#' detected. If \code{NULL}, all plotIds in \code{TP} are used.
#' @param checkEdges Before fitting the local regression should a check be done
#' if the first and last time point for a plot are outlying observations?
#' @param confIntSize A numeric value defining the confidence interval (see
#' Details).
#' @param nnLocfit A numeric value defining the constant component of the
#' smoothing parameter nn (see Details).
#'
#' @returns An object of class singleOut, a \code{data.frame} with the following
#' columns.
#' \describe{
#'   \item{plotId}{plotId}
#'   \item{timePoint}{time point}
#'   \item{trait}{modeled trait}
#'   \item{yPred}{prediction from the local regression}
#'   \item{sd_yPred}{standard deviation of the prediction}
#'   \item{lwr}{lower bound of the confidence interval}
#'   \item{upr}{upper bound of the confidence interval}
#'   \item{outlier}{flag for detected outlier (a value of 1 indicates the
#'   observation is an outlier)}
#' }
#'
#' @examples
#' ## Create a TP object containing the data from the Phenovator.
#' PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in%
#'                                  c("c24r41", "c7r18", "c7r49"), ]
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## First select a subset of plants, for example here 9 plants
#' plantSel <- phenoTP[[1]]$plotId[1:9]
#' # Then run on the subset
#' resuVatorHTP <- detectSingleOut(TP = phenoTP,
#'                                 trait = "EffpsII",
#'                                 plotIds = plantSel,
#'                                 confIntSize = 3,
#'                                 nnLocfit = 0.1)
#'
#' @family functions for detecting outliers for single observations
#'
#' @export
detectSingleOut <- function(TP,
                            trait,
                            plotIds = NULL,
                            checkEdges = TRUE,
                            confIntSize = 5,
                            nnLocfit = 0.5) {
  ## Checks.
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  TPTot <- do.call(rbind, TP)
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!hasName(x = TPTot, name = trait)) {
    stop("TP should contain a column ", trait, ".\n")
  }
  if (!is.null(plotIds)) {
    if (!all(plotIds %in% TPTot[["plotId"]])) {
      stop("All plotIds should be in TP.\n")
    }
    ## Restrict TPTot to selected plotIds
    TPTot <- TPTot[TPTot[["plotId"]] %in% plotIds, ]
  }
  ## Remove missing values for trait.
  TPTotRest <- na.omit(TPTot[c("plotId", "timePoint", trait)])
  ## Split data into a list of data.frames per plot.
  TPPlot <- split(x = TPTotRest, f = TPTotRest[["plotId"]], drop = TRUE)
  ## Compute outliers per plot.
  plotPreds <- lapply(X = TPPlot, FUN = function(plotDat) {
    ## Only makes sense for at least 6 time points.
    if (nrow(plotDat) <= 6) {
      warning("Not enough data points (at least 7) to fit a model for: ",
              plotDat[1, "plotId"], ".\n", call. = FALSE)
      return(NULL)
    }
    ## Fit a model using locfit.
    y <- plotDat[[trait]]
    x <- plotDat[["timePoint"]]
    ## Fit model with all time points included.
    fitMod <- locfit::locfit(y ~ locfit::lp(x, nn = nnLocfit, deg = 2))
    yPred <- predict(fitMod, newdata = x, se.fit = TRUE)
    if (checkEdges) {
      ## Check if first timepoint is an outlier.
      ## Fit model excluding first timepoint.
      fitMod0 <- locfit::locfit(y[-1] ~ locfit::lp(x[-1], nn = nnLocfit, deg = 2))
      ## Get predictions for new models.
      yPred0 <- predict(fitMod0, newdata = x, se.fit = TRUE)
      ## Compute mean and standard deviation for first 5 time points.
      m0 <- mean(y[1:5])
      s0 <- sqrt(sum((y[1:5] - m0) ^ 2) / 5)
      ## Remove first time point if it is outside mean +/- 1.5 sd.
      if (y[1] < m0 -  1.5 * s0 || y[1] > m0 + 1.5 * s0) {
        y <- y[-1]
        x <- x[-1]
        yPred <- yPred0
      }
      ## Check if last timepoint is an outlier.
      ## Get position of last timepoint.
      posL <- length(x)
      ## Fit model excluding last timepoint.
      fitMod1 <- locfit::locfit(y[-posL] ~ locfit::lp(x[-posL],
                                                      nn = nnLocfit, deg = 2))
      ## Get predictions for the model.
      yPred1 <- predict(fitMod1, newdata = plotDat[["timePoint"]], se.fit = TRUE)
      ## Compute mean and standard deviation for last 5 time points.
      m1 <- mean(y[(posL - 4):posL])
      s1 <- sqrt(sum((y[(posL - 4):posL] - m1) ^ 2) / 5)
      ## Remove last time point if it is outside mean +/- 1.5 sd.
      if (y[posL] < m1 - 1.5 * s1 || y[posL] > m1 + 1.5 * s1) {
        yPred <- yPred1
      }
    }
    ## Compute upper and lower boundaries.
    lwr <- yPred$fit - confIntSize * yPred$se.fit
    upr <- yPred$fit + confIntSize * yPred$se.fit
    ## Add results to plotDat.
    plotDat[["yPred"]] <- yPred$fit
    plotDat[["sd_yPred"]] <- yPred$se.fit
    plotDat[["lwr"]] <- lwr
    plotDat[["upr"]] <- upr
    ## All values outside the [lwr, upr] range are marked as outlier.
    plotDat[["outlier"]] <- ifelse(plotDat[[trait]] < plotDat[["upr"]] &
                                     plotDat[[trait]] > plotDat[["lwr"]], 0, 1)
    return(plotDat)
  })
  ## Bind everything together in a single data.frame.
  plotPred <- do.call(rbind, plotPreds)
  if (is.null(plotPred)) {
    stop("Not enough data points (<= 6) for any of the plots.\n")
  }
  ## Rownames are confusing and redundant.
  rownames(plotPred) <- NULL
  ## Add class and trait as attribute.
  class(plotPred) <- c("singleOut", class(plotPred))
  attr(plotPred, which = "trait") <- trait
  return(plotPred)
}

#' Plot outliers for single observations
#'
#' Plot the fitted local regression, confidence intervals and detected outliers
#' for each plotId.
#'
#' @inheritParams detectSingleOut
#' @inheritParams plot.TP
#'
#' @param x An object of class singleOut.
#' @param ... Ignored.
#' @param outOnly Should only plots containing outliers be plotted?
#'
#' @returns A list of ggplot objects is invisibly returned.
#'
#' @examples
#' ## Create a TP object containing the data from the Phenovator.
#' PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in%
#'                                  c("c24r41", "c7r18", "c7r49"), ]
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## Select a subset of plants, for example here 9 plants.
#' plantSel <- phenoTP[[1]]$plotId[1:9]
#' # Then run on the subset.
#' resuVatorHTP <- detectSingleOut(TP = phenoTP,
#'                                trait = "EffpsII",
#'                                plotIds = plantSel,
#'                                confIntSize = 3,
#'                                nnLocfit = 0.1)
#'
#' ## Visualize the prediction by choosing a single plant...
#' plot(resuVatorHTP, plotIds = "c21r24", outOnly = FALSE)
#' ## ...or a subset of plants.
#' plot(resuVatorHTP, plotIds = plantSel, outOnly = FALSE)
#'
#' @family functions for detecting outliers for single observations
#'
#' @export
plot.singleOut <- function(x,
                           ...,
                           plotIds = NULL,
                           outOnly = TRUE,
                           output = TRUE) {
  plotDat <- x
  if (!is.null(plotIds)) {
    if (!all(plotIds %in% plotDat[["plotId"]])) {
      stop("All plotIds should be in x\n")
    }
    plotDat <- plotDat[plotDat[["plotId"]] %in% plotIds, ]
  }
  ## Select outliers.
  outliers <- plotDat[plotDat[["outlier"]] == 1, ]
  if (outOnly) {
    plotDat <- plotDat[plotDat[["plotId"]] %in% outliers[["plotId"]], ]
    if (nrow(plotDat) == 0) {
      stop("No outliers present for selected plotIds.\n")
    }
  }
  plotDat <- droplevels(plotDat)
  trait <- attr(x = x, which = "trait")
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes(x = .data[["timePoint"]],
                                    y = .data[[trait]])) +
    ggplot2::geom_point(na.rm = TRUE)  +
    ggplot2::geom_line(mapping = ggplot2::aes(y = .data[["lwr"]]),
                       col = "green", linewidth = .8) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = .data[["upr"]]),
                       col = "green", linewidth = .8) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = .data[["yPred"]]),
                       col = "red", linewidth = .8) +
    ggplot2::geom_point(data = outliers, col = "blue", size = 2) +
    ggplot2::theme(legend.position = "none")
  nTp <- length(unique(plotDat[["timePoint"]]))
  if (nTp < 5) {
    p <- p + ggplot2::scale_x_datetime(breaks = unique(plotDat[["timePoint"]]),
                                       labels = scales::date_format("%B %d"))
  } else {
    ## Format the time scale to Month + day.
    p <- p + ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                       labels = scales::date_format("%B %d"))
  }
  ## Calculate the total number of plots.
  nPlots <- length(unique(plotDat[["plotId"]]))
  ## 25 Plots per page.
  nPag <- ceiling(nPlots / 25)
  if (nPlots >= 25) {
    ## More than 25 plots.
    ## For identical layout on all pages use 5 x 5 plots throughout.
    rowPag <- colPag <- rep(x = 5, times = nPag)
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
      ggforce::facet_wrap_paginate(facets = "plotId", nrow = rowPag[i],
                                   ncol = colPag[i], page = i)
    if (output) {
      suppressMessages(plot(pPag[[i]]))
    }
  }
  invisible(pPag)
}

#' Replace outliers for single observations by NA
#'
#' Function for replacing outliers for single observations by NA.
#'
#' @param TP An object of class TP.
#' @param singleOut A data.frame with at least the columns plotId and
#' timePoint with values corresponding to those in TP. If a column outlier is
#' present, as in the output of \code{detectSingleOut}, only plotId x
#' timePoint combinations for which outlier = 1 will be set to NA. If no
#' column outlier is present, all observations in singleOut will be set to NA.
#' @param trait The trait that should be set to NA. Can be ignored when using
#' the output of \code{detectSingleOut} as input.
#'
#' @returns An object of class TP, the input with the outlier replaced by NA.
#'
#' @examples
#' ## Create a TP object containing the data from the Phenovator.
#' PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in%
#'                                  c("c24r41", "c7r18", "c7r49"), ]
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## First select a subset of plants, for example here 9 plants.
#' plantSel <- phenoTP[[1]]$plotId[1:9]
#' # Then run on the subset
#' resuVatorHTP <- detectSingleOut(TP = phenoTP,
#'                                 trait = "EffpsII",
#'                                 plotIds = plantSel,
#'                                 confIntSize = 3,
#'                                 nnLocfit = 0.1)
#'
#' ## Replace the studied trait by NA for the plants marked as outliers.
#' phenoTPOut <- removeSingleOut(phenoTP, resuVatorHTP)
#'
#' @family functions for detecting outliers for single observations
#'
#' @export
removeSingleOut <- function(TP,
                            singleOut,
                            trait = attr(x = singleOut,
                                         which = "trait")) {
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  if (!inherits(singleOut, "data.frame")) {
    stop("singleOut should be a data.frame.\n")
  }
  if (!all(hasName(singleOut, c("plotId", "timePoint")))) {
    stop("singleOut should at least contain the columns plotId ",
         "and timePoint.\n")
  }
  if (!all(as.character(singleOut[["timePoint"]]) %in% names(TP))) {
    stop("All time points in singleOut should be in TP.\n")
  }
  if (!any(singleOut[["outlier"]] == 1)) {
    stop("There are no outlying points in TP.\n")
  }
  if (hasName(x = singleOut, "outlier")) {
    ## Remove observations that are not actually outliers.
    singleOut <- singleOut[singleOut[["outlier"]] == 1, ]
  }
  for (i in seq_len(nrow(singleOut))) {
    plotI <- singleOut[i, "plotId"]
    timeI <- as.character(singleOut[i, "timePoint"])
    if (!plotI %in% TP[[timeI]][["plotId"]]) {
      warning(plotI, " not present in timePoint ", timeI, ".\n", call. = FALSE)
    }
    TP[[timeI]][TP[[timeI]][["plotId"]] == plotI, trait] <- NA
  }
  return(TP)
}

