#' @keywords internal
createFitMod <- function(models,
                         what,
                         useRepId,
                         timePoints) {
  fitMod <- structure(models,
                      timePoints = timePoints,
                      what = what,
                      useRepId = useRepId,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

#' Plot function for class fitMod
#'
#' Plotting function for objects of class fitMod. A detailed description and
#' optional extra parameters of the different plots is given in the sections
#' below.
#'
#' @section rawPred plot:
#' Plots the raw data overlayed with the predicted values from the fitted model.
#' For each genotype a plot is made per plot/plant over time. These plots are
#' put together in a 5x5 grid. By using the parameter \code{genotypes} a
#' selection of genotypes can be plotted. Extra parameter options:
#' \describe{
#' \item{genotypes}{A character vector indicating the genotypes to be plotted.}
#' }
#'
#' @section corrPred plot:
#' Plots the spatially corrected data overlayed with the predicted values from
#' the fitted model. For each genotype a plot is made per plot/plant over time.
#' These plots are put together in a 5x5 grid. By using the parameter
#' \code{genotypes} a selection of genotypes can be plotted. Extra parameter
#' options:
#' \describe{
#' \item{genotypes}{A character vector indicating the genotypes to be plotted.}
#' }
#'
#' @section herit plot:
#' Plots the heritability over time. If \code{geno.decomp} is used when fitting
#' the model, heritabilities are plotted for each level of geno.decomp in a
#' single plot. Extra parameter options:
#' \describe{
#' \item{yLim}{A numerical vector of length two used for setting the limits of
#' the y-axis of the plot. If values outside of the plotting range are given,
#' then these are ignored.}
#' }
#'
#' @section effDim plot:
#' Plots the effective dimension from models fitted using SpATS over time.
#' Extra parameter options:
#' \describe{
#' \item{whichED}{A character vector indicating which effective dimensions
#' shoul be plotted. This should be a subset of "colId", "rowId", "fCol",
#' "fRow", "fColRow", "colfRow", "fColfRow" and "surface". Default all
#' effective dimensions are plotted.}
#' \item{yLim}{A numerical value used for setting the upper limit of the y-axis
#' of the plot. If a value lower than the highest value to be plotted is
#' given, then it is ignored.}
#' }
#'
#' @section variance plot:
#' Plots the residual, column and row variance for the fitted model over time.
#' Extra parameter options:
#' \describe{
#' \item{yLim}{A numerical value used for setting the upper limit of the y-axis
#' of the plot. If a value lower than the highest value to be plotted is
#' given, then it is ignored.}
#' }
#'
#' @section timeLapse plot:
#' Creates a time lapse of the spatial trends of models fitted using SpATS over
#' time.
#'
#' @section spatial plot:
#' Create a series of six spatial plots. Extra parameter options:
#' \describe{
#' \item{spaTrend}{A character string indicating how the spatial trend should
#' be displayed. Either "raw" for raw values, or "percentage" for displaying
#' as a percentage of the original phenotypic values.}
#' }
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class fitMod.
#' @param title A character string used as title for the plot. If \code{NULL} a
#' default title is added to the plot depending on \code{plotType}.
#' @param outFile A character string indicating the .pdf file or .gif file
#' (For \code{plotType} = "timeLapse") to which the plots should be written.
#' @param outFileOpts A named list of extra options for the pdf outfile, e.g.
#' width and height. See \code{\link[grDevices]{pdf}} for all possible options.
#'
#' @return Depending on the plot type either a ggplot object or a list of
#' ggplot objects is invisibly returned.
#'
#' @import ggplot2
#' @export
plot.fitMod <- function(x,
                        ...,
                        plotType = c("rawPred", "corrPred", "herit", "effDim",
                                     "variance", "timeLapse", "spatial"),
                        timePoints = names(x),
                        title = NULL,
                        output = TRUE,
                        outFile = NULL,
                        outFileOpts = NULL) {
  ## Checks.
  timePoints <- chkTimePoints(x, timePoints)
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (!is.null(title) && (!is.character(title) || length(title) > 1)) {
    stop("title should be NULL or a character string.\n")
  }
  ## Restrict x to selected time points.
  fitMods <- x[timePoints]
  ## Get engine from fitted models.
  engine <- class(fitMods[[1]])
  if (engine == "asreml" && plotType == "effDim") {
    stop("Effective dimensions can only be plotted for models fitted with SpATS")
  }
  ## Get geno.decomp from fitted models.
  if (engine == "SpATS") {
    geno.decomp <- fitMods[[1]]$model$geno$geno.decomp
  } else if (engine == "asreml") {
    if ("geno.decomp" %in% all.vars(fitMods[[1]]$formulae$random)) {
      geno.decomp <- "geno.decomp"
    } else {
      geno.decomp <- NULL
    }
  }
  ## Get trait from fitted models.
  if (engine == "SpATS") {
    trait <- fitMods[[1]]$model$response
  } else if (engine == "asreml") {
    ## Trait is always the only lhs variable in the fixed part.
    trait <- all.vars(update(fitMods[[1]]$formulae$fixed, .~0))
  }
  ## Get check from fitted models.
  if (engine == "SpATS") {
    useCheck <- grepl(pattern = "check", x = deparse(fitMods[[1]]$model$fixed))
  } else if (engine == "asreml") {
    useCheck <- "check" %in% all.vars(update(fitMods[[1]]$formulae$fixed, 0~.))
  }
  if (!is.null(outFile) && plotType != "timeLapse") {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    on.exit(dev.off(), add = TRUE)
    do.call(pdf, args = outFileOpts)
  }
  if (plotType == "rawPred") {
    genotypes <- dotArgs$genotypes
    if (!is.null(genotypes) && !is.character(genotypes)) {
      stop("genotypes should be NULL or a character vector.\n")
    }
    if (is.null(title)) title <- "Genotypic predictions + raw data"
    ## Get genotypic predictions.
    preds <- getGenoPred(fitMods)
    ## Construct full raw data from models.
    if (engine == "SpATS") {
      raw <- Reduce(f = rbind, x = lapply(X = fitMods, FUN = `[[`, "data"))
    } else if (engine == "asreml") {
      raw <- Reduce(f = rbind, x = lapply(X = fitMods, FUN = function(fitMod) {
        fitMod$call$data
      }))
    }
    if (!is.null(geno.decomp) && engine == "SpATS") {
      ## Genotype was converted to an interaction term of genotype and
      ## geno.decomp in the proces of fitting the model. That needs to be
      ## undone to get the genotype back in the output again.
      genoStart <- nchar(as.character(raw[["geno.decomp"]])) + 2
      raw[["genotype"]] <- as.factor(substring(raw[["genotype"]],
                                               first = genoStart))
    }
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      if (!all(genotypes %in% preds[["genotype"]])) {
        stop("All genotypes should be in ", deparse(substitute(x)))
      }
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      raw <- raw[raw[["genotype"]] %in% genotypes, ]
      raw <- droplevels(raw)
    }
    ## Restrict raw to neccessary columns so adding missing combinations works
    ## properly. With extra columns unneeded combinations might be added.
    raw <- raw[colnames(raw) %in% c("timePoint", "genotype", "plotId", trait,
                                    if (!is.null(geno.decomp)) "geno.decomp",
                                    if (useCheck) c("check", "genoCheck"))]
    ## Add combinations missing in data to raw.
    raw <- addMissVals(dat = raw, trait = trait)
    p <- xyFacetPlot(baseDat = raw, overlayDat = preds, yVal = trait,
                     yValOverlay = "predicted.values",
                     facetVal = c("genotype",
                                  if (!is.null(geno.decomp)) "geno.decomp"),
                     title = title, yLab = trait, output = output)
  } else if (plotType == "corrPred") {
    genotypes <- dotArgs$genotypes
    if (!is.null(genotypes) && !is.character(genotypes)) {
      stop("genotypes should be NULL or a character vector.\n")
    }
    if (is.null(title)) title <- "Genotypic predictions + spatial corrected data"
    ## Get genotypic predictions.
    preds <- getGenoPred(fitMods)
    ## Get spatial corrected values.
    corrected <- suppressWarnings(getCorrected(fitMods))
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      if (!all(genotypes %in% preds[["genotype"]])) {
        stop("All genotypes should be in ", deparse(substitute(x)))
      }
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      corrected <- corrected[corrected[["genotype"]] %in% genotypes, ]
      corrected <- droplevels(corrected)
    }
    ## Add combinations missing in data to corrected.
    corrected <- addMissVals(dat = corrected, trait = "newTrait")
    p <- xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = "newTrait",
                     yValOverlay = "predicted.values",
                     facetVal = c("genotype",
                                  if (!is.null(geno.decomp)) "geno.decomp"),
                     title = title,
                     yLab = trait, output = output)
  } else if (plotType == "herit") {
    if (is.null(title)) title <- "Heritabilities"
    ## Get heritabilities.
    herit <- getHerit(fitMods)
    ## Convert to long format needed by ggplot.
    herit <- reshape2::melt(herit, measure.vars = setdiff(colnames(herit),
                                                          c("timeNumber",
                                                            "timePoint")),
                            variable.name = "herit", value.name = "h2")
    ## Manually modify limit of y-axis.
    yLim <- c(min(dotArgs$yLim[1], herit[["h2"]]),
              max(dotArgs$yLim[2], herit[["h2"]]))
    p <- ggplot(herit, aes_string(x = "timePoint", y = "h2",
                                  group = "herit", color = "herit")) +
      geom_line(size = 0.5, na.rm = TRUE) +
      geom_point(size = 3, na.rm = TRUE) +
      plotTheme() +
      ylim(yLim) +
      labs(title = title)
    if (output) {
      plot(p)
    }
  } else if (plotType == "effDim") {
    whichEDopts <- c("colId", "rowId", "fCol", "fRow", "fColRow", "colfRow",
                     "fColfRow", "surface")
    if (is.null(dotArgs$which)) {
      whichED <- whichEDopts
    } else {
      whichED <- match.arg(dotArgs$whichED, choices = whichEDopts,
                           several.ok = TRUE)
    }
    if (is.null(title)) title <- "Effective dimensions"
    ## Get effective dimensions.
    effDim <- getEffDims(fitMods)
    ## Convert to long format needed by ggplot.
    effDim <- reshape2::melt(effDim, measure.vars = whichED,
                             variable.name = "effDim", value.name = "ED")
    ## Manually modify limit of y-axis.
    yLim <- c(min(dotArgs$yLim[1], effDim[["ED"]]),
              max(dotArgs$yLim[2], effDim[["ED"]]))
    p <- ggplot(effDim, aes_string(x = "timePoint", y = "ED",
                                   group = "effDim", color = "effDim")) +
      geom_line(size = 0.5, na.rm = TRUE) +
      geom_point(size = 3, na.rm = TRUE) +
      plotTheme() +
      ylim(yLim) +
      labs(title = title, color = "Effective dimension")
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
    ## Manually modify limit of y-axis.
    yLim <- c(min(dotArgs$yLim[1], variance[["value"]]),
              max(dotArgs$yLim[2], variance[["value"]]))
    p <- ggplot(variance, aes_string(x = "timePoint", y = "value",
                                     group = "var", color = "var")) +
      geom_line(size = 0.5, na.rm = TRUE) +
      geom_point(size = 3, na.rm = TRUE) +
      scale_color_discrete(labels = c("Residual", "Columns", "Rows")) +
      plotTheme() +
      ylim(yLim) +
      labs(title = title, color = "variance",
           y = expression(sigma ^ 2))
    if (output) {
      plot(p)
    }
  } else if (plotType == "spatial") {
    p <- lapply(X = fitMods, FUN = spatPlot, output = output, ... = ...)
  } else if (plotType == "timeLapse") {
    chkFile(outFile, fileType = "gif")
    timeLapsePlot(fitMods, outFile = outFile, ...)
  }
  if (!plotType == "timeLapse") {
    invisible(p)
  }
}

#' Helper function for minimal plot theme.
#'
#' @noRd
#' @keywords internal
plotTheme <- function() {
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 0.5))
}

#' Helper function for creating spatial plots.
#'
#' @noRd
#' @keywords internal
spatPlot <- function(fitMod,
                     output = TRUE,
                     ...) {
  dotArgs <- list(...)
  ## Get plot type for spatial trend from args.
  spaTrend <- match.arg(dotArgs$spaTrend, choices = c("raw", "percentage"))
  ## Get engine from fitted model.
  engine <- class(fitMod)
  ## Extract data from model.
  modDat <- fitMod$data
  ## Extract time point from model data.
  timePoint <- modDat[["timePoint"]][1]
  ## Extract trait from model.
  trait <- fitMod$model$response
  ## Get geno.decomp from fitted models.
  geno.decomp <- fitMod$model$geno$geno.decomp
  ## Extract what from model.
  what <- ifelse(fitMod$model$geno$as.random, "random", "fixed")
  ## Extract raw data.
  raw <- modDat[c("genotype", trait, "rowNum", "colNum", geno.decomp)]
  ## Extract fitted values from model.
  fitted <- fitted(fitMod)
  ## Extract predictions (BLUEs or BLUPs) from model.
  pred <- predictGeno(fitMod)[c("genotype", "predicted.values", geno.decomp)]
  if (!is.null(geno.decomp) && engine == "SpATS") {
    ## Genotype was converted to an interaction term of genotype and
    ## geno.decomp in the proces of fitting the model. That needs to be
    ## undone to get the genotype back in the output again.
    genoStart <- nchar(as.character(raw[["geno.decomp"]])) + 2
    raw[["genotype"]] <- as.factor(substring(raw[["genotype"]],
                                             first = genoStart))
  }
  ## Create plot data by merging extracted data together and renaming some
  ## columns.
  plotDat <- cbind(raw, fitted)
  plotDat <- merge(plotDat, pred, by = c("genotype", geno.decomp))
  plotDat[["predicted.values"]][is.na(plotDat[["fitted"]])] <- NA
  plotDat[["resid"]] <- plotDat[[trait]] - plotDat[["fitted"]]
  ## Get limits for row and columns.
  yMin <- min(plotDat[["rowNum"]])
  yMax <- max(plotDat[["rowNum"]])
  xMin <- min(plotDat[["rowNum"]])
  xMax <- max(plotDat[["colNum"]])
  ## Execute this part first since it needs plotData without missings
  ## removed.
  ## Code mimickes code from SpATS package but is adapted to create a
  ## data.frame useable by ggplot.
  plotDat <- plotDat[order(plotDat[["colNum"]], plotDat[["rowNum"]]), ]
  nCol <- xMax - xMin + 1
  nRow <- yMax - yMin + 1
  p1 <- 100 %/% nCol + 1
  p2 <- 100 %/% nRow + 1
  ## Get spatial trend from SpATS object.
  spatTr <- SpATS::obtain.spatialtrend(fitMod, grid = c(nCol * p1, nRow * p2))
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
  plotDatSpat <- reshape2::melt(t(spatTrDat), varnames = c("colNum", "rowNum"))
  ## Add true values for columns and rows for plotting.
  plotDatSpat[["colNum"]] <- spatTr$col.p
  plotDatSpat[["rowNum"]] <- rep(x = spatTr$row.p, each = p1 * nCol)
  ## Remove missings from data.
  plotDatSpat <- remove_missing(plotDatSpat, na.rm = TRUE)
  ## Now missing values can be removed from plotDat.
  plotDat <- remove_missing(plotDat, na.rm = TRUE)
  ## Code taken from plot.SpATS and simplified.
  ## Set colors and legends.
  colors <- topo.colors(100)
  legends <- c("Raw data", "Fitted data", "Residuals",
               "Fitted Spatial Trend",
               ifelse(what == "fixed", "Genotypic BLUEs",
                      "Genotypic BLUPs"), "Histogram")
  ## Compute range of values in response + fitted data so same scale
  ## can be used over plots.
  zlim <- range(c(plotDat[[trait]], plotDat[["fitted"]]), na.rm = TRUE)
  ## Create empty list for storing plots
  plots <- vector(mode = "list")
  ## Create main plot title.
  plotTitle <- ifelse(!is.null(dotArgs$title), dotArgs$title,
                      paste("Timepoint:", timePoint, "Trait:", trait))
  ## Create separate plots.
  plots$p1 <- fieldPlot(plotDat = plotDat, fillVar = trait,
                        title = legends[1], colors = colors, zlim = zlim)
  plots$p2 <- fieldPlot(plotDat = plotDat, fillVar = "fitted",
                        title = legends[2], colors = colors, zlim = zlim)
  plots$p3 <- fieldPlot(plotDat = plotDat, fillVar = "resid",
                        title = legends[3], colors = colors)
  if (spaTrend == "raw") {
    ## Get tickmarks from first plot to be used as ticks.
    ## Spatial plot tends to use different tickmarks by default.
    xTicks <- ggplot_build(plots[[1]])$layout$panel_params[[1]]$x.major_source
    plots$p4 <- fieldPlot(plotDat = plotDatSpat, fillVar = "value",
                          title = legends[4], colors = colors,
                          xTicks = xTicks)
  } else {
    phenoMean <- mean(modDat[[trait]], na.rm = TRUE)
    plotDatSpat[["value"]] <- plotDatSpat[["value"]] / phenoMean
    zlim <- c(-1, 1) * max(c(abs(plotDatSpat[["value"]]), 0.1))
    plots$p4 <- fieldPlotPcts(plotDat = plotDatSpat, fillVar = "value",
                              title = legends[4], zlim = zlim,
                              colors = colorRampPalette(c("red", "yellow", "blue"),
                                                        space = "Lab")(100))
  }
  plots$p5 <- fieldPlot(plotDat = plotDat, fillVar = "predicted.values",
                        title = legends[5], colors = colors)
  plots$p6 <- ggplot(data = plotDat) +
    geom_histogram(aes_string(x = "predicted.values"),
                   fill = "white", col = "black", bins = 10,
                   boundary = 0) +
    ## Remove empty space between ticks and actual plot.
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ## No background. Center and resize title. Resize axis labels.
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 9)) +
    labs(y = "Frequency", x = legends[5], title = legends[6])
  if (output) {
    ## do.call is needed since grid.arrange doesn't accept lists as input.
    do.call(gridExtra::grid.arrange,
            args = c(Filter(f = Negate(f = is.null), x = plots),
                     list(ncol = 3, top = plotTitle)))
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
                      xTicks = waiver(),
                      ...) {
  p <- ggplot(data = plotDat, aes_string(x = "colNum", y = "rowNum",
                                         fill = fillVar)) +
    geom_raster() +
    ## Remove empty space between ticks and actual plot.
    scale_x_continuous(expand = c(0, 0), breaks = xTicks) +
    scale_y_continuous(expand = c(0, 0)) +
    ## Adjust plot colors.
    scale_fill_gradientn(limits = zlim, colors = colors) +
    ## No background. Center and resize title. Resize axis labels.
    ## Remove legend title and resize legend entries.
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 9),
          legend.title = element_blank(),
          legend.text = element_text(size = 8, margin = margin(l = 5))) +
    ggtitle(title)
  return(p)
}


#' Helper function for creating field plots.
#'
#' @noRd
#' @keywords internal
fieldPlotPcts <- function(plotDat,
                          fillVar,
                          title,
                          colors,
                          zlim = range(plotDat[fillVar]),
                          xTicks = waiver(),
                          scaleLim = Inf,
                          ...) {
  p <- ggplot(
    data = plotDat,
    aes_string(x = "colNum", y = "rowNum",
               fill = fillVar,
               color = if (is.infinite(scaleLim)) NULL else "''")) +
    geom_tile(na.rm = TRUE) +
    ## Remove empty space between ticks and actual plot.
    scale_x_continuous(expand = c(0, 0), breaks = xTicks) +
    scale_y_continuous(expand = c(0, 0)) +
    ## Adjust plot colors.
    scale_fill_gradientn(limits = zlim, colors = colors, name = NULL,
                         labels = scales::percent,
                         breaks = seq(zlim[1], zlim[2],
                                      length.out = 5)) +
    scale_color_manual(values = NA) +
    guides(
      fill = guide_colorbar(order = 1),
      color = guide_legend("Larger than scale limit",
                           override.aes = list(fill = "grey50",
                                               color = "grey50"))) +
    ## No background. Center and resize title. Resize axis labels.
    ## Remove legend title and resize legend entries.
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title = element_text(size = 9)) +
    ggtitle(title)
  return(p)
}

#' Helper function for creating time lapse plots.
#'
#' @importFrom grDevices colorRampPalette
#' @noRd
#' @keywords internal
timeLapsePlot <- function(fitMods,
                          outFile = "spatialTrends.gif",
                          ...) {
  dotArgs <- list(...)
  if (!is.null(dotArgs$scaleLim)) {
    scaleLim <- dotArgs$scaleLim / 100
  } else {
    scaleLim <- Inf
  }
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
        spatTr[["fit"]]
      ## Melt to get the data in ggplot shape. Rows and columns in the
      ## spatial trend coming from SpATS are swapped so therefore use t()
      plotDatSpat <- reshape2::melt(t(spatTrDat),
                                    varnames = c("colNum", "rowNum"))
      ## Add true values for columns and rows for plotting.
      plotDatSpat[["colNum"]] <- spatTr[["col.p"]]
      plotDatSpat[["rowNum"]] <- rep(x = spatTr[["row.p"]], each = p1 * nCol)
      ## Remove missings from data.
      plotDatSpat <- remove_missing(plotDatSpat, na.rm = TRUE)
      ## Divide by mean value of trait to get trend as percentage.
      plotDatSpat[["mean"]] <- mean(plotDat[[trait]], na.rm = TRUE)
      plotDatSpat[["value"]] <- plotDatSpat[["value"]] / plotDatSpat[["mean"]]
      ## Set values outside scale limits to NA.
      plotDatSpat[["value"]][plotDatSpat[["value"]] > scaleLim |
                               plotDatSpat[["value"]] < -scaleLim] <- NA
      return(plotDatSpat)
    })
    ## Extract all zVals to use identical limits for spatial pattern
    ## This enables proper comparison of plots over timePoints.
    zVals <- unlist(sapply(X = plotSpatDats, `[[`, "value"))
    if (is.infinite(scaleLim)) {
      zLim <- c(min(c(zVals, -0.1), na.rm = TRUE),
                max(c(zVals, 0.1), na.rm = TRUE))
    } else {
      zLim <- c(-scaleLim, scaleLim)
    }
    ## Create a plot of the spatial trend per time point.
    for (i in seq_along(plotSpatDats)) {
      p <- fieldPlotPcts(plotDat = plotSpatDats[[i]], fillVar = "value",
                         title = paste(trait, names(plotSpatDats)[i]),
                         colors = colorRampPalette(c("red", "yellow", "blue"),
                                                   space = "Lab")(100),
                         zlim = zLim, scaleLim = scaleLim)
      plot(p)
    }
  }, movie.name = outFile, autobrowse = FALSE, loop = 1)
}

#' Function for extracting objects of class fitMod that keeps class.
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
    attr(r, "what") <- attr(x, "what")
    attr(r, "useRepId") <- attr(x, "useRepId")
    attr(r, "class") <- c("fitMod", "list")
    attr(r, "timestamp") <- attr(x, "timestamp")
  } else {
    r <- NULL
  }
  return(r)
}

