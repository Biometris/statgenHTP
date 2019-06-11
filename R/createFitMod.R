#' @keywords internal
createFitMod <- function(models,
                         timePoints) {
  fitMod <- structure(models,
                      timePoints = timePoints,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

#' @export
plot.fitMod <- function(x,
                        ...,
                        plotType = c("rawPred", "corrPred", "herit", "effDim",
                                     "variance", "rowPred", "colPred",
                                     "timeLapse"),
                        timePoints = names(x),
                        genotypes = NULL,
                        title = NULL) {
  ## Checks.
  if (!is.character(timePoints) || !all(hasName(x = x, name = timePoints))) {
    stop(paste0("All timePoints should be in ", deparse(substitute(x)), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (plotType == "rawPred") {
    if (is.null(title)) title <- "Genomic predictions + raw data"
    ## Get trait from model.
    trait <- x[[1]]$model$response
    ## Get genomic predictions.
    preds <- getGenoPred(x)
    ## Construct full raw data from models.
    raw <- Reduce(f = rbind, x = lapply(x, `[[`, "data"))
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      raw <- raw[raw[["genotype"]] %in% genotypes, ]
      raw <- droplevels(raw)
    }
    ## Add combinations missing in data to raw.
    raw <- addMissVals(dat = raw, trait = trait)
    xyFacetPlot(baseDat = raw, overlayDat = preds, yVal = trait,
                yValOverlay = "predicted.values", title = title, yLab = trait)
  } else if (plotType == "corrPred") {
    if (is.null(title)) title <- "Genomic predictions + spatial corrected data"
    ## Get trait from model.
    trait <- x[[1]]$model$response
    ## Get genomic predictions.
    preds <- getGenoPred(x)
    ## Get spatial corrected values.
    corrected <- getCorrected(x)
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      corrected <- corrected[corrected[["genotype"]] %in% genotypes, ]
      corrected <- droplevels(corrected)
    }
    ## Add combinations missing in data to corrected.
    corrected <- addMissVals(dat = corrected, trait = "newTrait")
    xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = "newTrait",
                yValOverlay = "predicted.values", title = title, yLab = trait)
  } else if (plotType == "herit") {
    if (is.null(title)) title <- "Heritabilities"
    ## Get heritabilities.
    herit <- getHerit(x)
    ## Convert to long format needed by ggplot.
    herit <- reshape2::melt(herit, measure.vars = setdiff(colnames(herit),
                                                          "timePoint"),
                            variable.name = "herit", value.name = "h2")
    ggplot2::ggplot(herit,
                    ggplot2::aes_string(x = "timePoint", y = "h2",
                                        group = "herit", color = "herit")) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title)
  } else if (plotType == "effDim") {
    if (is.null(title)) title <- "Effective dimensions"
    ## Get effective dimensions.
    effDim <- getEffDims(x)
    ## Convert to long format needed by ggplot.
    effDim <- reshape2::melt(effDim, measure.vars = c("effDimSurface",
                                                      "effDimCol", "effDimRow"),
                             variable.name = "effDim", value.name = "ED")
    ggplot2::ggplot(effDim,
                    ggplot2::aes_string(x = "timePoint", y = "ED",
                                        group = "effDim", color = "effDim")) +
      ggplot2::geom_line() +
      ggplot2::scale_color_discrete(labels = c("Surface", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title, color = "Effective dimension")
  } else if (plotType == "variance") {
    if (is.null(title)) title <- "Variances"
    ## Get variances.
    variance <- getVar(x)
    ## Convert to long format needed by ggplot.
    variance <- reshape2::melt(variance, measure.vars = c("varRes", "varCol",
                                                          "varRow"),
                               variable.name = "var")
    ggplot2::ggplot(variance,
                    ggplot2::aes_string(x = "timePoint", y = "value",
                                        group = "var", color = "var")) +
      ggplot2::geom_line() +
      ggplot2::scale_color_discrete(labels = c("Residual", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title, color = "variance",
                    y = expression(sigma ^ 2))
  } else if (plotType == "timeLapse") {
    outFile <- dotArgs$outFile
    timeLapsePlot(x, outFile = outFile)
  }
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



