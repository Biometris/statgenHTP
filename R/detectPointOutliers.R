#' detectPointOutliers
#'
#' Function to model each curve of a dataset using a local regression. This is the
#' first step of the detection of outlying points in a curve.
#'
#' see locfit() help function from the locfit R library. The user can act on:
#' \describe{
#'   \item{mylocfit}{the constant of the smoothing parameter. Increase mylocfit
#'   to have a very smooth curve}
#'   \item{confIntSize}{the level to calculate the confidence interval. Increase
#'   confIntSize to exclude less outliers}
#' }
#' to produce the grahics of the prediction and detected outliers, please use
#' plotDetectPointOutlierLocFit() function.
#'
#' @param TP An object of class TP.
#' @param trait A character vector indicating the trait to model in TP.
#' @param plotIds A character vector of plotIds for which the outliers should be
#' detected. If \code{NULL}, all plotId in TP are used.
#' @param confIntSize A numeric value defining the confidence interval.
#' @param mylocfit A numeric value defining the constant component of the smoothing parameter nn.
#' (see the locfit())
#'
#' @return An object of class pointOutliers, a data.frame with the following
#' columns.
#' \describe{
#'   \item{plotId}{the id variable}
#'   \item{timePoint}{time variable in datain}
#'   \item{trait}{name of the modeled variable in datain}
#'   \item{yPred}{the locfit prediction}
#'   \item{sd_yPred}{standard deviation of the prediction}
#'   \item{lwr}{lower bound of the confidence interval}
#'   \item{upr}{upper bound of the confidence interval}
#'   \item{outlier}{flag of detected outlier (0 is outlier, 1 is not)}
#' }
#'
#' @seealso plot.pointOutliers
#'
#' @examples ## Create a TP object containing the data from the Phenovator.
#' PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in% c("c24r41", "c7r18","c7r49"),]
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1","check2","check3","check4"))
#'
#' ## First select a subset of plants, for example here 10 plants
#' plantSel <- phenoTP[[1]]$plotId[1:9]
#' # Then run on the subset
#' resuVatorHTP <- detectPointOutliers(TP = phenoTP,
#'                                     trait = "EffpsII",
#'                                     plotIds = plantSel,
#'                                     confIntSize = 3,
#'                                     mylocfit = 0.1)
#'
#' @export
detectPointOutliers <- function(TP,
                                trait,
                                plotIds = NULL,
                                confIntSize = 5,
                                mylocfit = 0.5) {
  TPTot <- do.call(rbind, TP)
  if (!is.null(plotIds)) {
    if (!all(plotIds %in% TPTot[["plotId"]])) {
      stop("All plotIds should be in TP.\n")
    }
    TPTot <- TPTot[TPTot[["plotId"]] %in% plotIds, ]
  }
  TPTotRest <- na.omit(TPTot[c("plotId", "timePoint", trait)])
  TPPlot <- split(x = TPTotRest, f = TPTotRest[["plotId"]], drop = TRUE)
  plotPreds <- lapply(X = TPPlot, FUN = function(plotDat) {
    if (nrow(plotDat) <= 6) {
      warning("Not enough data points (<= 6) to fit a model for: ",
              plotDat[1, "plotId"], ".\n", call. = FALSE)
      return(NULL)
    }
    y <- plotDat[[trait]]
    x <- plotDat[["timePoint"]]
    fitMod <- locfit::locfit(y ~ locfit::lp(x, nn = mylocfit, deg = 2))
    ## Retrieving predictions for the x input interval.
    yPred <- predict(fitMod, newdata = x, se.fit = TRUE)
    lwr <- yPred$fit - confIntSize * yPred$se.fit
    upr <- yPred$fit + confIntSize * yPred$se.fit
    ## Add results to plotDat.
    plotDat[["yPred"]] <- yPred$fit
    plotDat[["sd_yPred"]] <- yPred$se.fit
    plotDat[["lwr"]] <- lwr
    plotDat[["upr"]] <- upr
    plotDat[["outlier"]] <- ifelse(y < upr & y > lwr, 0, 1)
    return(plotDat)
  })
  plotPred <- do.call(rbind, plotPreds)
  rownames(plotPred) <- NULL
  class(plotPred) <- c("pointOutliers", class(plotPred))
  attr(plotPred, which = "trait") <- trait
  return(plotPred)
}

#' plotDetectPointOutlierLocFit
#'
#' graphical function to produce the modeled smoothing and detected outliers
#' for each curve of a dataset using a local regression
#'
#' @inheritParams detectPointOutliers
#' @inheritParams plot.TP
#'
#' @param x An object of class pointOutliers.
#' @param outOnly Should only plots containing outliers be plotted?
#'
#' @examples ## Create a TP object containing the data from the Phenovator.
#' PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in% c("c24r41", "c7r18","c7r49"),]
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1","check2","check3","check4"))
#'
#' ## First select a subset of plants, for example here 10 plants
#' plantSel <- phenoTP[[1]]$plotId[1:9]
#' # Then run on the subset
#' resuVatorHTP <- detectPointOutliers(TP = phenoTP,
#'                                     trait = "EffpsII",
#'                                     plotIds = plantSel,
#'                                     confIntSize = 3,
#'                                     mylocfit = 0.1)
#'
#' ## We can then visualize the prediction by choosing a single plant...
#' plot(resuVatorHTP, plotIds = "c21r24", outOnly = FALSE)
#' ## ...or a subset of plants.
#' plot(resuVatorHTP, plotIds = plantSel, outOnly = FALSE)
#'
#' @export
plot.pointOutliers <- function(x,
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
  ## Compute the number of breaks for the time scale.
  ## If there are less than 3 time points use the number of time points.
  ## Otherwise use 3.
  nBr <- min(length(unique(plotDat[["timePoint"]])), 3)
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes_string(x = "timePoint", y = trait)) +
    ggplot2::geom_point(na.rm = TRUE)  +
    ggplot2::geom_line(mapping = ggplot2::aes_string(y = "lwr"),
                       col = "green", size = .8) +
    ggplot2::geom_line(mapping = ggplot2::aes_string(y = "upr"),
                       col = "green", size = .8) +
    ggplot2::geom_line(mapping = ggplot2::aes_string(y = "yPred"),
                       col = "red", size = .8) +
    ggplot2::geom_point(data = outliers, col = "blue", size = 2) +
    ggplot2::theme(legend.position = "none") +
    ## Format the time scale to Month + day.
    ggplot2::scale_x_datetime(breaks = prettier(n = nBr),
                              labels = scales::date_format("%B %d"))
  ## Calculate the total number of plots.
  nPlots <- length(unique(plotDat[["plotId"]]))
  ## 25 Plots per page.
  nPag <- ceiling(nPlots / 25)
  if (nPlots >= 25) {
    ## More than 25 plots.
    ## For identical layout on all pages use 5 x 5 plots throughout.
    #rowPag <- colPag <- rep(x = 5, times = nPag)

    # 28-7-2020. ggforce has a bug that prevents this identical layout
    # https://github.com/thomasp85/ggforce/issues/201
    # When fixed the code above can be reactivated and the three lines below
    # removed.
    plotsLastPag <- nPlots %% 25
    rowPag <- c(rep(x = 5, times = nPag - 1), min(plotsLastPag %/% 5 + 1, 5))
    colPag <- c(rep(x = 5, times = nPag - 1),
                ifelse(plotsLastPag >= 5, 5, plotsLastPag))
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
      ggforce::facet_wrap_paginate(facets = "plotId", nrow = rowPag[i],
                                   ncol = colPag[i], page = i)
    if (output) {
      suppressMessages(plot(pPag[[i]]))
    }
  }
  invisible(pPag)
}

#' Remove point outliers
#'
#' Function for setting point outliers to NA.
#'
#' @param TP An object of class TP.
#' @param pointOutliers A data.frame with at least the columns plotId and
#' timePoint with values corresponding to those in TP. If a column outlier is
#' present, as in the output of \code{detectPointOutliers}, only plot x time
#' combinations for which outlier = 1 will be set to NA. If no column outlier is
#' present, all observations in pointOutliers will be set to NA.
#' @param trait The trait that should be set to NA. Can be ignored when using
#' the output of \code{detectPointOutliers} as input.

#' @examples ## Create a TP object containing the data from the Phenovator.
#' PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in% c("c24r41", "c7r18","c7r49"),]
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1","check2","check3","check4"))
#'
#' ## First select a subset of plants, for example here 10 plants
#' plantSel <- phenoTP[[1]]$plotId[1:9]
#' # Then run on the subset
#' resuVatorHTP <- detectPointOutliers(TP = phenoTP,
#'                                     trait = "EffpsII",
#'                                     plotIds = plantSel,
#'                                     confIntSize = 3,
#'                                     mylocfit = 0.1)
#'
#' ## The annotated points can be replaced by NA for the studied trait
#' phenoTPOut <- removePointOutliers(phenoTP, resuVatorHTP)
#'
#' @export
removePointOutliers <- function(TP,
                                pointOutliers,
                                trait = attr(x = pointOutliers,
                                             which = "trait")) {
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  if (!inherits(pointOutliers, "data.frame")) {
    stop("pointOutliers should be a data.frame.\n")
  }
  if (!all(hasName(pointOutliers, c("plotId", "timePoint")))) {
    stop("pointOutliers should at least contain the columns plotId ",
         "and timePoint.\n")
  }
  if (!all(as.character(pointOutliers[["timePoint"]]) %in% names(TP))) {
    stop("All time points in pointOutliers should be in TP.\n")
  }
  if (!any(pointOutliers[["outlier"]] == 1)) {
    stop("There are no outlying points in TP.\n")
  }
  if (hasName(x = pointOutliers, "outlier")) {
    ## Remove observations that are not actually outliers.
    pointOutliers <- pointOutliers[pointOutliers[["outlier"]] == 1, ]
  }
  for (i in 1:nrow(pointOutliers)) {
    plotI <- pointOutliers[i, "plotId"]
    timeI <- as.character(pointOutliers[i, "timePoint"])
    if (!plotI %in% TP[[timeI]][["plotId"]]) {
      warning(plotI, " not present in timePoint ", timeI, ".\n", call. = FALSE)
    }
    TP[[timeI]][TP[[timeI]][["plotId"]] == plotI, trait] <- NA
  }
  return(TP)
}

