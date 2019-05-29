#' @keywords internal
createFitMod <- function(models) {
  fitMod <- structure(models,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

plot.fitMod <- function(x,
                        ...,
                        plotType = c("rawPred", "corrPred", "herit", "effDim",
                                     "variance", "rowPred", "colPred"),
                        timePoints = names(x),
                        traits = NULL,
                        output = TRUE) {
  ## Checks.
  if (!is.character(timePoints) || !all(hasName(x = x, name = timePoints))) {
    stop(paste0("All timePoints should be in ", deparse(substitute(x)), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (plotType == "rawPred") {
    preds <- getBLUPs(x)
    raw <- Reduce(f = rbind, x = lapply(x, `[[`, "data"))
    xyFacetPlot(baseDat = raw, overlayDat = preds, yVal = traits,
                yValOverlay = "predicted.values",
                title = "Genomic predictions + raw data",
                yLab = traits)
  } else if (plotType == "corrPred") {
    preds <- getBLUPs(x)
    corrected <- getCorrected(x)
    xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = "newTrait",
                yValOverlay = "predicted.values",
                title = "Genomic predictions + spatial corrected data",
                yLab = traits)
  } else if (plotType == "herit") {
    herit <- getHerit(x)
    ggplot2::ggplot(herit, ggplot2::aes_string(x = "timePoint", y = "h2")) +
      ggplot2::geom_line() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "Heritabilities")
  } else if (plotType == "effDim") {
    effDim <- getEffDims(x)
    effDim <- reshape2::melt(effDim, measure.vars = c("effDimSurface",
                                                      "effDimCol", "effDimRow"),
                   variable.name = "effDim", value.name = "ED")
    ggplot2::ggplot(effDim,
                    ggplot2::aes_string(x = "timePoint", y = "ED",
                                        group = "effDim", color = "effDim")) +
      ggplot2::geom_line() +
      ggplot2::scale_color_discrete(labels = c("Surface", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "Effective dimensions",
                    color = "Effective dimension")
  } else if (plotType == "variance") {
    variance <- getVar(x)
    variance <- reshape2::melt(variance, measure.vars = c("varRes", "varCol",
                                                          "varRow"),
                               variable.name = "var")
    ggplot2::ggplot(variance,
                    ggplot2::aes_string(x = "timePoint", y = "value",
                                        group = "var", color = "var")) +
      ggplot2::geom_line() +
      ggplot2::scale_color_discrete(labels = c("Residual", "Columns", "Rows")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = "Variances", color = "variance",
                    y = expression(sigma ^ 2))
  }
}
