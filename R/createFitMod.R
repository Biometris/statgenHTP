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
                                     "rowPred", "colPred"),
                        timePoints = names(x),
                        traits = NULL,
                        output = TRUE) {
  ## Checks.
  if (!is.character(timePoints) || !all(hasName(x = x, name = timePoints))) {
    stop(paste0("All timePoints should be in ", deparse(substitute(x)), ".\n"))
  }
  plotType <- match.arg(plotType, several.ok = TRUE)
  dotArgs <- list(...)
  if ("corrPred" %in% plotType) {
    preds <- getBLUPs(x)
    preds[["time"]] <- lubridate::ymd_hms(preds[["time"]])
    corrected <- getCorrected(x)
    xyFacetPlot(baseDat = corrected, overlayDat = preds, yVal = "newTrait",
                yValOverlay = "predicted.values",
                title = "Genomic predictions",
                yLab = traits)
  }
}
