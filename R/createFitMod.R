#' @keywords internal
createFitMod <- function(models) {
  fitMod <- structure(models,
                      class = c("fitMod", "list"),
                      timestamp = Sys.time())
  return(fitMod)
}

plot.fitMod <- function(x,
                        ...,
                        plotType = c("layout", "box", "cor", "raw"),
                        timePoints = names(x),
                        traits = NULL,
                        output = TRUE) {
  ## Checks.
  if (!is.character(timePoints) || !all(hasName(x = x, name = timePoints))) {
    stop(paste0("All timePoints should be in ", deparse(substitute(x)), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
}
