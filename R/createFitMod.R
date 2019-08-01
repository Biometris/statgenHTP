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
#' Plotting function for objects of class fitMod.
#'
#' @section rawPred plot:
#' Plots the raw data overlayed with the predicted values from the fitted model.
#' For each genotype a plot is made per plot/plant over time. These plots are
#' put together in a 5x5 grid. By using the parameter \code{genotypes} a
#' selection of genotypes can be plotted.
#'
#' @section corrPred plot:
#' Plots the spatially corrected data overlayed with the predicted values from
#' the fitted model. For each genotype a plot is made per plot/plant over time.
#' These plots are put together in a 5x5 grid. By using the parameter
#' \code{genotypes} a selection of genotypes can be plotted.
#'
#' @section herit plot:
#' Plots the heritability over time. If \code{geno.decomp} is used when fitting
#' the model, heritabilities are plotted for each level of geno.decomp in a
#' single plot.
#'
#' @section effDim plot:
#' Plots the effective dimension from models fitted using SpATS over time.
#'
#' @section variance plot:
#' Plots the residual, column and row variance for the fitted model over time.
#'
#' @section timeLapse plot:
#' Creates a time lapse of the spatial trends of models fitted using SpATS over
#' time.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class fitMod.
#' @param whichED A character vector indicating which effective dimensions
#' shoul be plotted. Only used if \code{plotType} = "effDim".
#' @param genotypes A character vector indicating the genotypes to be plotted.
#' Only used if \code{plotType} = "rawPred" or "corrPred".
#' @param title A character string used as title for the plot. If \code{NULL} a
#' default title is added to the plot depending on \code{plotType}.
#' @param outFile A character string indicating the .pdf file or .gif file
#' (For \code{plotType} = "timeLapse") to which the plots should be written.
#' @param outFileOpts A named list of extra options for the pdf outfile, e.g.
#' width and height. See \code{\link[grDevices]{pdf}} for all possible options.
#'
#' @return Depending on the plot type either a ggplot object or a list of ggplot
#' objects is invisibly returned.
#'
#' @export
plot.fitMod <- function(x,
                        ...,
                        plotType = c("rawPred", "corrPred", "herit", "effDim",
                                     "variance", "timeLapse"),
                        whichED = c("colId", "rowId", "col", "row", "rowCol",
                                   "fCol", "fRow", "fColRow", "colfRow",
                                   "fColfRow", "surface"),
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
  if (!is.null(genotypes) && !is.character(genotypes)) {
    stop("genotypes should be NULL or a character vector.\n")
  }
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
    if (is.null(title)) title <- "Genotypic predictions + spatial corrected data"
    ## Get genotypic predictions.
    preds <- getGenoPred(fitMods)
    ## Get spatial corrected values.
    corrected <- suppressWarnings(getCorrected(fitMods))
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
    whichED <- match.arg(whichED, several.ok = TRUE)
    if (is.null(title)) title <- "Effective dimensions"
    ## Get effective dimensions.
    effDim <- getEffDims(fitMods)
    ## Convert to long format needed by ggplot.
    effDim <- reshape2::melt(effDim, measure.vars = whichED,
                             variable.name = "effDim", value.name = "ED")
    p <- ggplot2::ggplot(effDim,
                         ggplot2::aes_string(x = "timePoint", y = "ED",
                                             group = "effDim", color = "effDim")) +
      ggplot2::geom_line() +
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
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::scale_color_discrete(labels = c("Residual", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = title, color = "variance",
                    y = expression(sigma ^ 2))
    if (output) {
      plot(p)
    }
  } else if (plotType == "timeLapse") {
    chkFile(outFile, fileType = "gif")
    timeLapsePlot(fitMods, outFile = outFile, ...)
  }
  if (!plotType == "timeLapse") {
    invisible(p)
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
                      scaleLim = Inf,
                      ...) {
  p <- ggplot2::ggplot(
    data = plotDat,
    ggplot2::aes_string(x = "colNum", y = "rowNum",
                        fill = fillVar,
                        color = if (is.infinite(scaleLim)) NULL else "''")) +
    ggplot2::geom_tile(na.rm = TRUE) +
    ## Remove empty space between ticks and actual plot.
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = xTicks) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ## Adjust plot colors.
    ggplot2::scale_fill_gradientn(limits = zlim, colors = colors, name = NULL,
                                  labels = scales::percent,
                                  breaks = seq(zlim[1], zlim[2],
                                               length.out = 5)) +
    ggplot2::scale_color_manual(values = NA) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1),
      color = ggplot2::guide_legend("Larger than scale limit",
                                    override.aes = list(fill = "grey50",
                                                        color = "grey50"))) +
    ## No background. Center and resize title. Resize axis labels.
    ## Remove legend title and resize legend entries.
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.title = ggplot2::element_text(size = 9)) +
    ggplot2::ggtitle(title)
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
      plotDatSpat <- ggplot2::remove_missing(plotDatSpat, na.rm = TRUE)
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
      p <- fieldPlot(plotDat = plotSpatDats[[i]], fillVar = "value",
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

