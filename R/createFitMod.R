#' @keywords internal
createFitMod <- function(models,
                         experimentName,
                         what,
                         useRepId,
                         useCheck,
                         spatial,
                         timePoints) {
  fitMod <- structure(models,
                      experimentName = experimentName,
                      timePoints = timePoints,
                      what = what,
                      useRepId = useRepId,
                      useCheck = useCheck,
                      spatial = spatial,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

#' Summary function for fitMod objects
#'
#' Function for creating a short summary of the contents of a TP object. The
#' summary consists of the name of the experiment, the number of time points,
#' the engine used to fit the models and, in case spatial models where fitted
#' using asreml, the selected spatial model.
#'
#' @param object An object of class fitMod.
#' @param ... Ignored.
#'
#' @return No return value, a summary is printed.
#'
#' @examples
#' \donttest{
#' ## Using the first example dataset (PhenovatorDat1):
#' ## Create an object of class TP.
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
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 36))
#'
#' ## Create a summary.
#' summary(modPhenoSp)
#' }
#'
#' @family functions for spatial modeling
#'
#' @export
summary.fitMod <- function(object,
                           ...)  {
  experimentName <- attr(x = object, which = "experimentName")
  noTP <- length(object)
  engine <- class(object[[1]])
  if (engine == "asreml" && attr(x = object, which = "spatial")) {
    tpUsed <- min(max(noTP / 5, 10), noTP)
    bestSpat <- attr(x = object[[1]], which = "sumTab")[[1]][, "spatial"]
  }
  cat("Models in ", deparse(substitute(object)),
      " where fitted for experiment ", experimentName, ".\n\n", sep = "")
  cat("It contains", noTP, "time points.\n")
  cat("The models were fitted using ", engine, ".\n\n", sep = "")
  if (engine == "asreml" && attr(x = object, which = "spatial")) {
    cat("The selected spatial model is ", bestSpat, ".\n", sep = "")
    cat(tpUsed, "time points were used to select the best spatial model.\n")
  }
}

#' Plot function for class fitMod
#'
#' Plotting function for objects of class \code{fitMod}. Seven different types
#' of plots can be made for an object of class \code{fitMod}. A detailed
#' description and optional extra parameters for the different plots is given
#' in the sections below.
#'
#' @section rawPred plot:
#' Plots the raw data (colored dots) overlayed with the predicted values from
#' the fitted model (black dots). For each genotype a plot is made per
#' plot/plant over time. These plots are put together in a 5x5 grid. By using
#' the parameter \code{genotypes} a selection of genotypes can be plotted.
#' Extra parameter options:
#' \describe{
#' \item{genotypes}{A character vector indicating the genotypes to be plotted.}
#' \item{plotChecks}{Should the check genotypes be included in the plot?}
#' \item{plotLine}{Should the data be displayed as lines? Default is
#' \code{FALSE}.}
#' }
#'
#' @section corrPred plot:
#' Plots the spatially corrected data (colored dots) overlayed with the
#' predicted values from the fitted model (black dors). For each genotype a plot
#' is made per plot/plant over time. These plots are put together in a 5x5 grid.
#' By using the parameter \code{genotypes} a selection of genotypes can be
#' plotted. Extra parameter options:
#' \describe{
#' \item{genotypes}{A character vector indicating the genotypes to be plotted.}
#' \item{plotChecks}{Should the check genotypes be included in the plot?}
#' \item{plotLine}{Should the data be displayed as lines? Default is
#' \code{FALSE}.}
#' }
#'
#' @section herit plot:
#' Plots the heritability over time. This plot is only available when genotype
#' is fitted as random factor in the model. If \code{geno.decomp} is used when
#' fitting the model, heritabilities are plotted for each level of geno.decomp
#' in a single plot. Extra parameter options:
#' \describe{
#' \item{yLim}{A numerical vector of length two, used for setting the limits of
#' the y-axis of the plot. If values outside of the plotting range are given,
#' then these are ignored.}
#' }
#'
#' @section effDim plot:
#' Plots the effective dimension over time for models fitted using SpATS.
#' Extra parameter options:
#' \describe{
#' \item{whichED}{A character vector indicating which effective dimensions
#' should be plotted. This should be a subset of "colId", "rowId", "fCol",
#' "fRow", "fColRow", "colfRow", "fColfRow" and "surface". When
#' \code{useRepId = TRUE}, the effective dimensions of "colId" and "rowId"
#' become "RepId:colId" and "RepId:rowId". Default all effective dimensions
#' are plotted.}
#' \item{EDType}{A character string specifying if the effective dimension
#' ("dimension") or the ratio of effective dimensions ("ratio") should be
#' plotted. Default the dimensions are plotted.}
#' \item{yLim}{A numerical vector of length two, used for setting the limits of
#' the y-axis of the plot. If values outside of the plotting range are given,
#' then these are ignored.}
#' }
#'
#' @section variance plot:
#' Plots the residual, column and row variances over time for the fitted models.
#' Extra parameter options:
#' \describe{
#' \item{yLim}{A numerical vector of length two, used for setting the limits of
#' the y-axis of the plot. If values outside of the plotting range are given,
#' then these are ignored.}
#' }
#'
#' @section timeLapse plot:
#' Creates a time lapse of the spatial trends of models fitted using SpATS over
#' time.
#'
#' @section spatial plot:
#' Creates five plots per time point, spatial plots of the raw data,
#' fitted values, residuals and either BLUEs or BLUPs, and a histogram of the
#' BLUEs or BLUPs. When SpATS was used for modeling an extra plot with the
#' fitted spatial trend is included Extra parameter options:
#' \describe{
#' \item{spaTrend}{A character string indicating how the spatial trend should
#' be displayed. Either "raw" for raw values, or "percentage" for displaying
#' as a percentage of the original phenotypic values.}
#' }
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class fitMod.
#' @param outFile A character string indicating the .pdf file or .gif file
#' (for \code{plotType} = "timeLapse") to which the plots should be written.
#'
#' @return Depending on the plot type either a ggplot object or a list of
#' ggplot objects is invisibly returned.
#'
#' @examples
#' \donttest{
#' ## Using the first example dataset (PhenovatorDat1):
#' ## Create an object of class TP.
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
#' ## Fit a SpATS model on three points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 36))
#'
#' ## Plot the spatial trends for one time point:
#' plot(modPhenoSp,
#'      timePoints = 36,
#'      plotType = "spatial",
#'      spaTrend = "percentage")
#' }
#'
#' \dontrun{
#' ## Create a time lapse of all available time points:
#' plot(modPhenoSp,
#'      plotType = "timeLapse",
#'      outFile = "TimeLapse_modPhenoSp.gif")
#' }
#'
#' \donttest{
#' ## Plot the corrected values for a subset of four genotypes:
#' plot(modPhenoSp,
#'      plotType = "corrPred",
#'      genotypes = c("check1", "check2", "G007", "G058") )
#'
#' ## Plot the effective dimensions of all available time points in the model
#' ## for a subset of effective dimensions:
#' plot(modPhenoSp,
#'      plotType = "effDim",
#'      whichED = c("colId", "rowId", "fColRow","colfRow"),
#'      EDType = "ratio")
#' }
#'
#' @family functions for spatial modeling
#'
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
  experimentName <- attr(x = x, which = "experimentName")
  dotArgs <- list(...)
  if (!is.null(title) && (!is.character(title) || length(title) > 1)) {
    stop("title should be NULL or a character string.\n")
  }
  ## Restrict x to selected time points.
  fitMods <- x[timePoints]
  ## Get engine from fitted models.
  engine <- class(fitMods[[1]])
  if (engine == "asreml" && plotType == "effDim") {
    stop("Effective dimensions can only be plotted for models fitted ",
         "with SpATS.\n")
  }
  if (engine == "asreml" && plotType == "spatial" &&
      !attr(x = fitMods, which = "spatial")) {
    stop("spatial plots can only be made when setting spatial = TRUE when ",
         "fitting the asreml models.\n")
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
  ## Get what from fitted models.
  what <- attr(x = fitMods, which = "what")
  if (!is.null(outFile) && plotType != "timeLapse") {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    on.exit(dev.off(), add = TRUE)
    do.call(pdf, args = outFileOpts)
  }
  if (plotType == "rawPred") {
    genotypes <- dotArgs$genotypes
    plotLine <- isTRUE(dotArgs$plotLine)
    if (!is.null(genotypes) && !is.character(genotypes)) {
      stop("genotypes should be NULL or a character vector.\n")
    }
    if (is.null(title)) title <-
      paste(experimentName, "- genotypic prediction + raw data")
    if (isTRUE(dotArgs$plotChecks) && useCheck) {
      totPred <- getGenoPred(fitMods, predictChecks = TRUE)
      ## Get check predictions.
      genoPred <- totPred$genoPred
      checkPred <- totPred$checkPred
      ## Rename check column to genotype so rbinding is possible.
      colnames(checkPred)[colnames(checkPred) == "check"] <- "genotype"
      preds <- rbind(genoPred, checkPred)
    } else {
      ## Get genotypic predictions.
      preds <- getGenoPred(fitMods)$genoPred
    }
    ## Construct full raw data from models.
    if (engine == "SpATS") {
      raw <- Reduce(f = rbind, x = lapply(X = fitMods, FUN = `[[`, "data"))
    } else if (engine == "asreml") {
      raw <- Reduce(f = rbind, x = lapply(X = fitMods, FUN = function(fitMod) {
        fitMod$call$data
      }))
    }
    ## Remove observations from raw where genotype is missing.
    ## These where included when fitting spatial asreml models to create
    ## a full grid.
    raw <- raw[!is.na(raw[["genotype"]]), ]
    if (!is.null(geno.decomp) && engine == "SpATS" && !useCheck) {
      ## Genotype was converted to an interaction term of genotype and
      ## geno.decomp in the proces of fitting the model. That needs to be
      ## undone to get the genotype back in the output again.
      genoStart <- nchar(as.character(raw[["geno.decomp"]])) + 2
      raw[["genotype"]] <- as.factor(substring(raw[["genotype"]],
                                               first = genoStart))
    }
    ## Remove check genotypes from raw data.
    if (!isTRUE(dotArgs$plotChecks) && useCheck) {
      raw <- droplevels(raw[!is.na(raw[["genoCheck"]]), ])
    }
    ## Restrict genotypes.
    if (!is.null(genotypes)) {
      if (!all(genotypes %in% preds[["genotype"]])) {
        stop("All genotypes should be in ", deparse(substitute(x)), ".\n")
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
                     title = title, yLab = trait, output = output,
                     plotLine = plotLine)
  } else if (plotType == "corrPred") {
    genotypes <- dotArgs$genotypes
    plotLine <- isTRUE(dotArgs$plotLine)
    if (!is.null(genotypes) && !is.character(genotypes)) {
      stop("genotypes should be NULL or a character vector.\n")
    }
    if (is.null(title))
      title <- paste(experimentName, "- genotypic prediction + corrected data")
    if (isTRUE(dotArgs$plotChecks) && useCheck) {
      totPred <- getGenoPred(fitMods, predictChecks = TRUE)
      ## Get check predictions.
      genoPred <- totPred$genoPred
      checkPred <- totPred$checkPred
      ## Rename check column to genotype so rbinding is possible.
      colnames(checkPred)[colnames(checkPred) == "check"] <- "genotype"
      preds <- rbind(genoPred, checkPred)
    } else {
      ## Get genotypic predictions.
      preds <- getGenoPred(fitMods)$genoPred
    }
    ## Get spatial corrected values.
    corrected <- suppressWarnings(getCorrected(fitMods))
    ## Remove check genotypes from corrected data.
    if (!isTRUE(dotArgs$plotChecks) && useCheck) {
      corrected <- droplevels(corrected[corrected[["check"]] == "noCheck", ])
    }
    ## Restrict genotypes.,
    if (!is.null(genotypes)) {
      if (!all(genotypes %in% preds[["genotype"]])) {
        stop("All genotypes should be in ", deparse(substitute(x)), ".\n")
      }
      preds <- preds[preds[["genotype"]] %in% genotypes, ]
      preds <- droplevels(preds)
      corrected <- corrected[corrected[["genotype"]] %in% genotypes, ]
      corrected <- droplevels(corrected)
    }
    newTrait <- paste0(trait, "_corr")
    corrected <- corrected[c("timeNumber", "timePoint", "genotype",
                             newTrait, "plotId",
                             if (!is.null(geno.decomp)) "geno.decomp")]
    ## Add combinations missing in data to corrected.
    corrected <- addMissVals(dat = corrected, trait = newTrait)
    p <- xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = newTrait,
                     yValOverlay = "predicted.values",
                     facetVal = c("genotype",
                                  if (!is.null(geno.decomp)) "geno.decomp"),
                     title = title,
                     yLab = trait, output = output, plotLine = plotLine)
  } else if (plotType == "herit") {
    if (is.null(title)) title <- paste(experimentName, "- Heritability")
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
    p <- ggplot2::ggplot(herit,
                         ggplot2::aes(x = .data[["timePoint"]],
                                      y = .data[["h2"]],
                                      group = .data[["herit"]],
                                      color = .data[["herit"]])) +
      ggplot2::geom_point(size = 3, na.rm = TRUE) +
      plotTheme() +
      ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0,
                                                          vjust = 0.5)) +
      ggplot2::ylim(yLim) +
      ggplot2::labs(title = title)
    ## Compute the number of breaks for the time scale.
    ## If there are less than 4 time points use positions of the time points.
    ## Otherwise use 3.
    nTp <- length(unique(herit[["timePoint"]]))
    if (nTp < 5) {
      p <- p + ggplot2::scale_x_datetime(breaks = unique(herit[["timePoint"]]),
                                         labels = scales::date_format("%B %d"))
    } else {
      ## Format the time scale to Month + day.
      p <- p + ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                         labels = scales::date_format("%B %d"))
    }
    if (nTp > 1) {
      p <- p + ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE)
    }
    if (output) {
      plot(p)
    }
  } else if (plotType == "effDim") {
    useRepId <- attr(x = fitMods, which = "useRepId")
    colVarId <- ifelse(useRepId, "repId:colId", "colId")
    rowVarId <- ifelse(useRepId, "repId:rowId", "rowId")
    whichEDopts <- c(colVarId, rowVarId, "fCol", "fRow", "fColRow", "colfRow",
                     "fColfRow", "surface")
    if (is.null(dotArgs$which)) {
      whichED <- whichEDopts
    } else {
      whichED <- match.arg(dotArgs$whichED, choices = whichEDopts,
                           several.ok = TRUE)
    }
    EDType <- match.arg(dotArgs$EDType, choices = c("dimension", "ratio"))
    if (is.null(title)) title <- paste(experimentName, "- Effective dimension")
    ## Get effective dimensions.
    effDim <- getEffDims(fitMods, EDType = EDType)
    ## Convert to long format needed by ggplot.
    effDim <- reshape2::melt(effDim, measure.vars = whichED,
                             variable.name = "effDim", value.name = "ED")
    ## Manually modify limit of y-axis.
    yLim <- c(min(dotArgs$yLim[1], effDim[["ED"]]),
              max(dotArgs$yLim[2], effDim[["ED"]]))
    p <- ggplot2::ggplot(effDim, ggplot2::aes(x = .data[["timePoint"]],
                                              y = .data[["ED"]],
                                              group = .data[["effDim"]],
                                              color = .data[["effDim"]])) +
      ggplot2::geom_point(size = 3, na.rm = TRUE) +
      plotTheme() +
      ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0,
                                                          vjust = 0.5)) +
      ggplot2::ylim(yLim) +
      ggplot2::labs(title = title, color = "Effective dimension")
    ## Compute the number of breaks for the time scale.
    ## If there are less than 4 time points use positions of the time points.
    ## Otherwise use 3.
    nTp <- length(unique(effDim[["timePoint"]]))
    if (nTp < 5) {
      p <- p + ggplot2::scale_x_datetime(breaks = unique(effDim[["timePoint"]]),
                                         labels = scales::date_format("%B %d"))
    } else {
      ## Format the time scale to Month + day.
      p <- p + ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                         labels = scales::date_format("%B %d"))
    }
    if (nTp > 1) {
      p <- p + ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE)
    }
    if (output) {
      plot(p)
    }
  } else if (plotType == "variance") {
    if (is.null(title)) title <- paste(experimentName, "- Variances")
    ## Get variances.
    variance <- getVar(fitMods)
    ## Get variance columns from variance, i.e. all columns starting with var.
    varCols <- colnames(variance)[grepl(pattern = "^var",
                                        x = colnames(variance))]
    ## Construct labels for variances.
    varLabs <- c(if ("varGen" %in% varCols) "Genotypic" else
      substring(varCols[1:(length(varCols) - 3)], first = 17),
      "Residual", "Columns", "Rows")
    ## Convert to long format needed by ggplot.
    variance <- reshape2::melt(variance, measure.vars = varCols,
                               variable.name = "var")
    ## Manually modify limit of y-axis.
    yLim <- c(min(dotArgs$yLim[1], variance[["value"]]),
              max(dotArgs$yLim[2], variance[["value"]]))
    p <- ggplot2::ggplot(variance, ggplot2::aes(x = .data[["timePoint"]],
                                                y = .data[["value"]],
                                                group = .data[["var"]],
                                                color = .data[["var"]])) +
      ggplot2::geom_point(size = 3, na.rm = TRUE) +
      ggplot2::scale_color_discrete(labels = varLabs) +
      plotTheme() +
      ggplot2::theme(axis.title.y = ggplot2::element_text(angle = 0,
                                                          vjust = 0.5)) +
      ggplot2::ylim(yLim) +
      ggplot2::labs(title = title, color = "variance",
                    y = expression(sigma ^ 2))
    ## Compute the number of breaks for the time scale.
    ## If there are less than 4 time points use positions of the time points.
    ## Otherwise use 3.
    nTp <- length(unique(variance[["timePoint"]]))
    if (nTp < 5) {
      p <- p + ggplot2::scale_x_datetime(breaks = unique(variance[["timePoint"]]),
                                         labels = scales::date_format("%B %d"))
    } else {
      ## Format the time scale to Month + day.
      p <- p + ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                         labels = scales::date_format("%B %d"))
    }
    if (nTp > 1) {
      p <- p + ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE)
    }
    if (output) {
      plot(p)
    }
  } else if (plotType == "spatial") {
    p <- lapply(X = fitMods, FUN = spatPlot, trait = trait, what = what,
                geno.decomp = geno.decomp, useCheck = useCheck, engine = engine,
                output = output, ... = ..., experimentName = experimentName,
                title = title)
  } else if (plotType == "timeLapse") {
    chkFile(outFile, fileType = "gif")
    timeLapsePlot(fitMods, outFile = outFile, ...)
  }
  if (!plotType == "timeLapse") {
    invisible(p)
  }
}

#' Helper function for creating spatial plots.
#'
#' @noRd
#' @keywords internal
spatPlot <- function(fitMod,
                     trait,
                     what,
                     geno.decomp,
                     useCheck,
                     engine,
                     output = TRUE,
                     ...) {
  dotArgs <- list(...)
  ## Get plot type for spatial trend from args.
  spaTrend <- match.arg(dotArgs$spaTrend, choices = c("raw", "percentage"))
  ## Extract data from model.
  if (engine == "SpATS") {
    modDat <- fitMod$data
  } else if (engine == "asreml") {
    modDat <- fitMod$call$data
  }
  ## Check spatial information in modDat.
  if (!chkRowCol(modDat)) {
    return(NULL)
  }
  ## Extract time point from model data.
  timePoint <- modDat[["timePoint"]][1]
  ## Extract raw data.
  raw <- modDat[c("genotype", trait, "rowNum", "colNum", geno.decomp)]
  ## Extract fitted values from model.
  fitted <- fitted(fitMod)
  ## Extract predictions (BLUEs or BLUPs) from model.
  if (useCheck) {
    totPred <- predictGeno(fitMod, predictChecks = TRUE)
    ## Get check predictions.
    genoPred <- totPred$predGeno
    checkPred <- totPred$predCheck
    ## Rename check column to genotype so rbinding is possible.
    colnames(checkPred)[colnames(checkPred) == "check"] <- "genotype"
    pred <- rbind(genoPred, checkPred)
  } else {
    ## Get genotypic predictions.
    pred <- predictGeno(fitMod)$predGeno
  }
  pred <- pred[c("genotype", "predicted.values", geno.decomp)]
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
  xMin <- min(plotDat[["colNum"]])
  xMax <- max(plotDat[["colNum"]])
  ## Execute this part first since it needs plotData without missings
  ## removed.
  ## Code mimickes code from SpATS package but is adapted to create a
  ## data.frame useable by ggplot.
  if (engine == "SpATS") {
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
    for (i in seq_len(nrow(plotDat))) {
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
    plotDatSpat[["colNum"]] <- spatTr$col.p
    plotDatSpat[["rowNum"]] <- rep(x = spatTr$row.p, each = p1 * nCol)
    ## Remove missings from data.
    plotDatSpat <- ggplot2::remove_missing(plotDatSpat, na.rm = TRUE)
  }
  ## Now missing values can be removed from plotDat.
  plotDat <- ggplot2::remove_missing(plotDat, na.rm = TRUE)
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
                      paste(dotArgs$experimentName, "-", trait, "-", timePoint))
  ## Create separate plots.
  plots$p1 <- fieldPlot(plotDat = plotDat, fillVar = trait,
                        title = legends[1], colors = colors, zlim = zlim)
  plots$p2 <- fieldPlot(plotDat = plotDat, fillVar = "fitted",
                        title = legends[2], colors = colors, zlim = zlim)
  plots$p3 <- fieldPlot(plotDat = plotDat, fillVar = "resid",
                        title = legends[3], colors = colors)
  ## Spatial plot only for SpATS.
  if (engine == "SpATS") {
    ## Get tickmarks from first plot to be used as ticks.
    ## Spatial plot tends to use different tickmarks by default.
    xTicks <- ggplot2::ggplot_build(plots[[1]])$layout$panel_params[[1]]$x$breaks
    if (spaTrend == "raw") {
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
                                                          space = "Lab")(100),
                                xTicks = xTicks)
    }
  }
  plots$p5 <- fieldPlot(plotDat = plotDat, fillVar = "predicted.values",
                        title = legends[5], colors = colors)
  plots$p6 <- ggplot2::ggplot(data = plotDat) +
    ggplot2::geom_histogram(ggplot2::aes(x = .data[["predicted.values"]]),
                            fill = "white", col = "black", bins = 10,
                            boundary = 0) +
    ## Remove empty space between ticks and actual plot.
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ## No background. Center and resize title. Resize axis labels.
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.title = ggplot2::element_text(size = 9)) +
    ggplot2::labs(y = "Frequency", x = legends[5], title = legends[6])
  if (output) {
    ## do.call is needed since grid.arrange doesn't accept lists as input.
    do.call(gridExtra::grid.arrange,
            args = c(Filter(f = Negate(f = is.null), x = plots),
                     list(ncol = 3, top = plotTitle)))
  }
  return(plots)
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
  p <- ggplot2::ggplot(data = plotDat, ggplot2::aes(x = .data[["colNum"]],
                                                    y = .data[["rowNum"]],
                                                    fill = .data[[fillVar]])) +
    ggplot2::geom_tile() +
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
                   legend.text = ggplot2::element_text(size = 8,
                                                       margin = ggplot2::margin(l = 5))) +
    ggplot2::ggtitle(title)
  return(p)
}


#' Helper function for creating field plots with percentages.
#'
#' @noRd
#' @keywords internal
fieldPlotPcts <- function(plotDat,
                          fillVar,
                          title,
                          colors,
                          zlim = range(plotDat[fillVar]),
                          xTicks = ggplot2::waiver(),
                          scaleLim = Inf,
                          ...) {
  p <- ggplot2::ggplot(
    data = plotDat,
    ggplot2::aes(x = .data[["colNum"]], y = .data[["rowNum"]],
                 fill = .data[[fillVar]],
                 color = if (is.infinite(scaleLim)) NULL else "")) +
    ggplot2::geom_tile(na.rm = TRUE) +
    ## Remove empty space between ticks and actual plot.
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = xTicks) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ## Adjust plot colors.
    ggplot2::scale_fill_gradientn(limits = zlim, colors = colors, name = NULL,
                                  labels = scales::percent,
                                  breaks = seq(zlim[1], zlim[2], length.out = 5)) +
    ggplot2::scale_color_manual(values = NA) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(order = 1),
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
      for (i in seq_len(nrow(plotDat))) {
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
      zLim <- c(-1, 1) * max(c(abs(zVals), 0.1), na.rm = TRUE)
    } else {
      zLim <- c(-scaleLim, scaleLim)
    }
    ## Create a plot of the spatial trend per time point.
    for (i in seq_along(plotSpatDats)) {
      p <- fieldPlotPcts(plotDat = plotSpatDats[[i]], fillVar = "value",
                         title = names(plotSpatDats)[i],
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
#' @noRd
#' @export
`[.fitMod` <- function(x,
                       i,
                       ...) {
  timePoints <- chkTimePoints(x, i)
  timePointsX <- attr(x, which = "timePoints")
  timePointsR <- timePointsX[timePointsX[["timePoint"]] %in% timePoints, ]
  if (nrow(timePointsR) > 0) {
    class(x) <- "list"
    r <- x[timePointsR[["timePoint"]]]
    attr(r, "timePoints") <- timePointsR
    attr(r, "experimentName") <- attr(x, "experimentName")
    attr(r, "what") <- attr(x, "what")
    attr(r, "useRepId") <- attr(x, "useRepId")
    attr(r, "useCheck") <- attr(x, "useCheck")
    attr(r, "spatial") <- attr(x, "spatial")
    attr(r, "class") <- c("fitMod", "list")
    attr(r, "timestamp") <- attr(x, "timestamp")
  } else {
    r <- NULL
  }
  return(r)
}

