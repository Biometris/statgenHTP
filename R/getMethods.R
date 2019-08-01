#' Extract genotypic predictions
#'
#' Extract predictions of the genotypic value from an object of class fitMod.
#'
#' @param fitMod An object of class fitMod.
#' @param timePoints A character or numeric vector indicating the timePoints
#' to be modeled. When using a character string to reference a timePoint, the
#' value has to be an exact match to one of the existing timePoints. When using
#' a number it will be matched by its number in the timePoints attribute of the
#' TP object.
#' @param outFile A character string indicating the .csv file to which the
#' results should be written. If \code{NULL} no file is written.
#'
#' @return A data.frame with genotypic predictions per time point.
#'
#' @export
getGenoPred <- function(fitMod,
                        timePoints = names(fitMod),
                        outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  ## Get predictions per time point.
  genoPred <- lapply(X = fitMod, FUN = predictGeno)
  ## Create one data.frame containing all time points.
  genoPred <- Reduce(f = rbind, x = genoPred)
  ## Add time numbers.
  genoPred <- addTimeNumber(fitMod, genoPred)
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "csv")
    write.csv(genoPred, file = outFile, row.names = FALSE)
  }
  return(genoPred)
}

#' Extract corrected phenotypic values
#'
#' Extract corrected phenotype from an object of class fitMod. After fitting a
#' spatial model at each time point, the raw phenotypic data is corrected by
#' subtracting the (estimated) sources of (environmental, design effect) which
#' are of no interest (nuisances). This allows keeping the data resolution at
#' the plot/plant level.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with spatially corrected values per time point.
#'
#' @export
getCorrected <- function(fitMod,
                         timePoints = names(fitMod),
                         outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  ## correctSpatial will throw warnings for every timepoint when no spatial
  ## or fixed effects are present. Catch these and only show once.
  spatCorrTP <- tryCatchExt(lapply(X = fitMod, FUN = correctSpatial))
  if (!is.null(spatCorrTP$error)) {
    stop(spatCorrTP$error)
  } else if (!is.null(spatCorrTP$warning)) {
    warning(unique(spatCorrTP$warning))
  }
  ## Create one data.frame with corrected values for all time points.
  spatCorr <- Reduce(f = rbind, x = spatCorrTP$value)
  ## Add time number.
  spatCorr <- addTimeNumber(fitMod, spatCorr)
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "csv")
    write.csv(spatCorr, file = outFile, row.names = FALSE)
  }
  return(spatCorr)
}

#' Extract variances
#'
#' Extract variances from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with variances per time point.
#'
#' @export
getVar <- function(fitMod,
                   timePoints = names(fitMod),
                   outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  useRepId <- attr(x = fitMod, which = "useRepId")
  colVarId <- ifelse(useRepId, "repId:colId", "colId")
  rowVarId <- ifelse(useRepId, "repId:rowId", "rowId")
  if (inherits(fitMod[[1]], "SpATS")) {
    varRes <- sapply(X = fitMod, FUN = function(x) x$psi[1])
    varCol <- sapply(X = fitMod, FUN = function(x) x$var.comp[colVarId])
    varRow <- sapply(X = fitMod, FUN = function(x) x$var.comp[rowVarId])
  } else if (inherits(fitMod[[1]], "asreml")) {
    varRes <- sapply(X = fitMod, FUN = function(x) x$sigma2)
    varCol <- sapply(X = fitMod, FUN = function(x) x$vparameters[colVarId])
    varRow <- sapply(X = fitMod, FUN = function(x) x$vparameters[rowVarId])
  }
  variance <- data.frame(timePoint = lubridate::as_datetime(names(varRes)),
                         varRes = varRes, varCol = varCol, varRow = varRow,
                         row.names = NULL)
  variance <- addTimeNumber(fitMod, variance)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(variance, file = outFile, row.names = FALSE)
  }
  return(variance)
}

#' Extract heritabilities
#'
#' Extract heritabilities from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with heritabilities per time point.
#'
#' @export
getHerit <- function(fitMod,
                     timePoints = names(fitMod),
                     outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (attr(x = fitMod, which = "what") == "fixed") {
    stop("Heritability can only be calculated when genotype is random.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  h2Out <- lapply(X = fitMod, FUN = heritability)
  h2Out <- dfBind(h2Out)
  h2Out <- addTimeNumber(fitMod, h2Out)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(h2Out, file = outFile, row.names = FALSE)
  }
  return(h2Out)
}

#' Extract effective dimensions
#'
#' Extract effective dimensions from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with effective dimensions per time point.
#'
#' @export
getEffDims <- function(fitMod,
                       timePoints = names(fitMod),
                       outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (!inherits(fitMod[[1]], "SpATS")) {
    stop("Models in ", fitMod, " should be fitted using SpATS.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  useRepId <- attr(x = fitMod, which = "useRepId")
  colVarId <- ifelse(useRepId, "repId:colId", "colId")
  rowVarId <- ifelse(useRepId, "repId:rowId", "rowId")
  effDimSurface <- sapply(X = fitMod, FUN = function(x) {
    sum(x$eff.dim[c("f(colNum)", "f(rowNum)", "f(colNum):rowNum",
                    "colNum:f(rowNum)","f(colNum):f(rowNum)")])
  })
  effDimCol <- sapply(X = fitMod, FUN = function(x) x$eff.dim[colVarId])
  effDimRow <- sapply(X = fitMod, FUN = function(x) x$eff.dim[rowVarId])
  effDim <- data.frame(timePoint = lubridate::as_datetime(names(effDimSurface)),
                       effDimSurface = effDimSurface,
                       effDimCol = effDimCol, effDimRow = effDimRow,
                       row.names = NULL)
  effDim <- addTimeNumber(fitMod, effDim)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(effDim, file = outFile, row.names = FALSE)
  }
  return(effDim)
}

#' Helper function for adding time numbers to data containing a time point
#' column.
#'
#' @noRd
#' @keywords internal
addTimeNumber <- function(fitMod,
                          dat) {
  ## Get data.frame containing time numbers and time points.
  timePoints <- attr(x = fitMod, which = "timePoints")
  ## Covert timePoint column to datetime format for merging.
  timePoints[["timePoint"]] <- lubridate::as_datetime(timePoints[["timePoint"]])
  ## Merge time number to data.
  dat <- merge(timePoints, dat)
  ## Put timeNumber as first column and timePoint as second.
  dat <- dat[c(2, 1, 3:ncol(dat))]
}

