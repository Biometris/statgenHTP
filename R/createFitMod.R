#' @keywords internal
createFitMod <- function(models,
                         timePoints) {
  fitMod <- structure(models,
                      timePoints = timePoints,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

#' Plot function for class fitMod
#'
#' Plotting function for objects of class fitMod.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class fitMod.
#' @param genotypes A character vector indicating the genotypes to be plotted.
#' Only used if \code{plotType} = "rawPred" or "corrPred".
#' @param title A character string used as title for the plot.
#' @param outFile A character string indicting the .pdf file to which the
#' plots should be written. If \code{NULL} no file is written.
#' @param outFileOpts A named list of extra options for the pdf outfile, e.g.
#' width and height. See \code{\link[grDevices]{pdf}} for all possible options.
#'
#' @export
plot.fitMod <- function(x,
                        ...,
                        plotType = c("rawPred", "corrPred", "herit", "effDim",
                                     "variance", "rowPred", "colPred",
                                     "timeLapse"),
                        timePoints = names(x),
                        genotypes = NULL,
                        title = NULL,
                        output = TRUE,
                        outFile = NULL,
                        outFileOpts = NULL) {
  ## Checks.
  timePoints <- chkTimePoints(x, timePoints)
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  ## Restrict x to selected time points.
  fitMods <- x[timePoints]
  if (!is.null(outFile) && plotType != "timeLapse") {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    do.call(pdf, args = outFileOpts)
  }
  if (plotType == "rawPred") {
    if (is.null(title)) title <- "Genomic predictions + raw data"
    ## Get trait from model.
    trait <- fitMods[[1]]$model$response
    ## Get genomic predictions.
    preds <- getGenoPred(fitMods)
    ## Construct full raw data from models.
    raw <- Reduce(f = rbind, x = lapply(fitMods, `[[`, "data"))
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      raw <- raw[raw[["genotype"]] %in% genotypes, ]
      raw <- droplevels(raw)
    }
    ## Add combinations missing in data to raw.
    raw <- addMissVals(dat = raw, trait = trait)
    p <- xyFacetPlot(baseDat = raw, overlayDat = preds, yVal = trait,
                     yValOverlay = "predicted.values", title = title,
                     yLab = trait, output = output)
  } else if (plotType == "corrPred") {
    if (is.null(title)) title <- "Genomic predictions + spatial corrected data"
    ## Get trait from model.
    trait <- fitMods[[1]]$model$response
    ## Get genomic predictions.
    preds <- getGenoPred(fitMods)
    ## Get spatial corrected values.
    corrected <- getCorrected(fitMods)
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      corrected <- corrected[corrected[["genotype"]] %in% genotypes, ]
      corrected <- droplevels(corrected)
    }
    ## Add combinations missing in data to corrected.
    corrected <- addMissVals(dat = corrected, trait = "newTrait")
    p <- xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = "newTrait",
                     yValOverlay = "predicted.values", title = title,
                     yLab = trait, output = output)
  } else if (plotType == "herit") {
    if (is.null(title)) title <- "Heritabilities"
    ## Get heritabilities.
    herit <- getHerit(fitMods)
    ## Convert to long format needed by ggplot.
    herit <- reshape2::melt(herit, measure.vars = setdiff(colnames(herit),
                                                          "timePoint"),
                            variable.name = "herit", value.name = "h2")
    p <- ggplot2::ggplot(herit,
                         ggplot2::aes_string(x = "timePoint", y = "h2",
                                             group = "herit", color = "herit")) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title)
    if (output) {
      plot(p)
    }
  } else if (plotType == "effDim") {
    if (is.null(title)) title <- "Effective dimensions"
    ## Get effective dimensions.
    effDim <- getEffDims(fitMods)
    ## Convert to long format needed by ggplot.
    effDim <- reshape2::melt(effDim, measure.vars = c("effDimSurface",
                                                      "effDimCol", "effDimRow"),
                             variable.name = "effDim", value.name = "ED")
    p <- ggplot2::ggplot(effDim,
                         ggplot2::aes_string(x = "timePoint", y = "ED",
                                             group = "effDim", color = "effDim")) +
      ggplot2::geom_line() +
      ggplot2::scale_color_discrete(labels = c("Surface", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title, color = "Effective dimension")
    if (output) {
      plot(p)
    }
  } else if (plotType == "variance") {
    if (is.null(title)) title <- "Variances"
    ## Get variances.
    variance <- getVar(fitMods)
    ## Convert to long format needed by ggplot.
    variance <- reshape2::melt(variance, measure.vars = c("varRes", "varCol",
                                                          "varRow"),
                               variable.name = "var")
    p <- ggplot2::ggplot(variance,
                         ggplot2::aes_string(x = "timePoint", y = "value",
                                             group = "var", color = "var")) +
      ggplot2::geom_line() +
      ggplot2::scale_color_discrete(labels = c("Residual", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title, color = "variance",
                    y = expression(sigma ^ 2))
    if (output) {
      plot(p)
    }
  } else if (plotType == "timeLapse") {
    chkFile(outFile, fileType = "gif")
    timeLapsePlot(fitMods, outFile = outFile)
  }
  if (!is.null(outFile) && plotType != "timeLapse") {
    dev.off()
  }
  invisible(p)
}

#' Helper function for creating field plots.
#'
#' @noRd
#' @keywords internal
fieldPlot <- function(plotDat,
                      fillVar,
                      title,
                      colors,
                      zlim = range(plotDat[fillVar]),
                      xTicks = ggplot2::waiver(),
                      ...) {
  p <- ggplot2::ggplot(data = plotDat,
                       ggplot2::aes_string(x = "colNum", y = "rowNum",
                                           fill = fillVar)) +
    ggplot2::geom_raster(na.rm = TRUE) +
    ## Remove empty space between ticks and actual plot.
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = xTicks) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ## Adjust plot colors.
    ggplot2::scale_fill_gradientn(limits = zlim, colors = colors) +
    ## No background. Center and resize title. Resize axis labels.
    ## Remove legend title and resize legend entries.
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.title = ggplot2::element_text(size = 9),
                   legend.title = ggplot2::element_blank(),
                   legend.text =
                     ggplot2::element_text(size = 8,
                                           margin = ggplot2::margin(l = 5))) +
    ggplot2::ggtitle(title)
  return(p)
}

#' Helper function for creating time lapse plots.
#'
#' @noRd
#' @keywords internal
timeLapsePlot <- function(fitMods,
                          outFile = "spatialTrends.gif") {
  animation::saveGIF({
    ## Get trait from fitted models.
    trait <- fitMods[[1]]$model$response
    ## First get plot data for all fields so a single zLim can be extracted
    ## for all plots.
    plotSpatDats <- lapply(X = fitMods, FUN = function(fitMod) {
      ## Get data from fitted model
      plotDat <- fitMod$data
      ## Get min and max values for rows and columns.
      yMin <- min(plotDat$rowNum)
      yMax <- max(plotDat$rowNum)
      xMin <- min(plotDat$colNum)
      xMax <- max(plotDat$colNum)
      ## Execute this part first since it needs plotData without missings
      ## removed.
      ## Code mimickes code from SpATS package but is adapted to create a
      ## data.frame useable by ggplot.
      plotDat <- plotDat[order(plotDat$colNum, plotDat$rowNum), ]
      nCol <- xMax - xMin + 1
      nRow <- yMax - yMin + 1
      p1 <- 100 %/% nCol + 1
      p2 <- 100 %/% nRow + 1
      ## Get spatial trend from SpATS object.
      spatTr <- SpATS::obtain.spatialtrend(fitMod,
                                           grid = c(nCol * p1, nRow * p2))
      ## spatial trend contains values for all data points, so NA in original
      ## data need to be removed. The kronecker multiplication is needed to
      ## convert the normal row col pattern to the smaller grid extending the
      ## missing values.
      ## First a matrix M is created containing information for all
      ## columns/rows in the field even if they are completely empty.
      M <- matrix(nrow = nRow, ncol = nCol,
                  dimnames = list(yMin:yMax, xMin:xMax))
      for (i in 1:nrow(plotDat)) {
        M[as.character(plotDat[i, "rowNum"]),
          as.character(plotDat[i, "colNum"])] <-
          ifelse(is.na(plotDat[i, trait]), NA, 1)
      }
      spatTrDat <- kronecker(M, matrix(data = 1, ncol = p1, nrow = p2)) *
        spatTr$fit
      ## Melt to get the data in ggplot shape. Rows and columns in the
      ## spatial trend coming from SpATS are swapped so therefore use t()
      plotDatSpat <- reshape2::melt(t(spatTrDat),
                                    varnames = c("colNum", "rowNum"))
      ## Add true values for columns and rows for plotting.
      plotDatSpat$colNum <- spatTr$col.p
      plotDatSpat$rowNum <- rep(x = spatTr$row.p, each = p1 * nCol)
      ## Remove missings from data.
      plotDatSpat <- ggplot2::remove_missing(plotDatSpat, na.rm = TRUE)
      return(plotDatSpat)
    })
    ## Extract all zVals to use identical limits for spatial pattern
    ## This enables proper comparison of plots over timePoints.
    zVals <- unlist(sapply(X = plotSpatDats, `[[`, "value"))
    zLim <- c(min(zVals), max(zVals))
    ## Create a plot of the spatial trend per time point.
    for (i in seq_along(plotSpatDats)) {
      p <- fieldPlot(plotDat = plotSpatDats[[i]], fillVar = "value",
                     title = paste(trait, names(plotSpatDats)[i]),
                     colors = topo.colors(100), zlim = zLim)
      plot(p)
    }
  }, movie.name = outFile, autobrowse = FALSE)
}

#' Function for extracting for objects of class fitMod that keeps class.
#'
#' @param x An object of class fitMod.
#' @param i An index specifying the element to extract of replace.
#' @param ... Ignored.
#'
#' @export
`[.fitMod` <- function(x, i, ...) {
  timePoints <- chkTimePoints(x, i)
  timePointsX <- attr(x, which = "timePoints")
  timePointsR <- timePointsX[timePointsX[["timePoint"]] %in% timePoints, ]
  if (nrow(timePointsR) > 0) {
    class(x) <- "list"
    r <- x[timePointsR[["timePoint"]]]
    attr(r, "timePoints") <- timePointsR
    attr(r, "class") <- c("fitMod", "list")
    attr(r, "timestamp") <- attr(x, "timestamp")
  } else {
    r <- NULL
  }
  return(r)
}

