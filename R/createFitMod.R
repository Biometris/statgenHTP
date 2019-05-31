#' @keywords internal
createFitMod <- function(models) {
  fitMod <- structure(models,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

#' @export
plot.fitMod <- function(x,
                        ...,
                        plotType = c("rawPred", "corrPred", "herit", "effDim",
                                     "variance", "rowPred", "colPred"),
                        timePoints = names(x),
                        title = NULL) {
  ## Checks.
  if (!is.character(timePoints) || !all(hasName(x = x, name = timePoints))) {
    stop(paste0("All timePoints should be in ", deparse(substitute(x)), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (plotType == "rawPred") {
    if (is.null(title)) title <- "Genomic predictions + raw data"
    trait <- x[[1]]$model$response
    preds <- getGenoPred(x)
    raw <- Reduce(f = rbind, x = lapply(x, `[[`, "data"))
    xyFacetPlot(baseDat = raw, overlayDat = preds, yVal = trait,
                yValOverlay = "predicted.values", title = title, yLab = trait)
  } else if (plotType == "corrPred") {
    if (is.null(title)) title <- "Genomic predictions + spatial corrected data"
    trait <- x[[1]]$model$response
    preds <- getGenoPred(x)
    corrected <- getCorrected(x)
    xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = "newTrait",
                yValOverlay = "predicted.values", title = title, yLab = trait)
  } else if (plotType == "herit") {
    if (is.null(title)) title <- "Heritabilities"
    herit <- getHerit(x)
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
      ggplot2::labs(title = title, color = "Effective dimension")
  } else if (plotType == "variance") {
    if (is.null(title)) title <- "Variances"
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
      ggplot2::labs(title = title, color = "variance",
                    y = expression(sigma ^ 2))
  }
}
