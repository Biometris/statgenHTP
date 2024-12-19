#' Extract time points
#'
#' Function for extracting a data.frame with timeNumbers and timePoints from
#' an object of class TP or fitMod.
#'
#' @param x An object of class TP or fitMod
#'
#' @returns A data.frame with columns timeNumber and timePoint listing the
#' time points in x
#'
#' @examples
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
#' ## Extract the time points from the object.
#' head(getTimePoints(phenoTP))
#'
#' @family functions for data preparation
#'
#' @export
getTimePoints <- function(x) {
  ## Check input.
  if (!inherits(x, "TP") && !inherits(x, "fitMod")) {
    stop("x should be an object of class TP or fitMod")
  }
  return(attr(x = x, which = "timePoints"))
}
