#' Extract genotypic predictions
#'
#' Extract predictions of the genotypic value from an object of class fitMod.
#'
#' @param fitMod An object of class fitMod.
#' @param outFile A character string indicating the .csv file to which the
#' results should be written. If \code{NULL} no file is written.
#'
#' @return A data.frame with genotypic predictions per time point.
#'
#' @export
getGenoPred <- function(fitMod,
                        outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  genoPred <- lapply(X = fitMod, FUN = predictGeno)
  genoPred <- Reduce(f = rbind, x = genoPred)
  genoPred <- addTimeNumber(fitMod, genoPred)
  if (!is.null(outFile)) {
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
                         outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  ## correctSpatial will throw warnings for every timepoint when no spatial
  ## or fixed effects are present. Catch these and only show once.
  spatCorrTP <- tryCatchExt(lapply(X = fitMod, FUN = correctSpatial))
  if (!is.null(spatCorrTP$error)) {
    stop(spatCorrTP$error)
  } else if (!is.null(spatCorrTP$warning)) {
    warning(unique(spatCorrTP$warning))
  }
  spatCorr <- Reduce(f = rbind, x = spatCorrTP$value)
  spatCorr <- addTimeNumber(fitMod, spatCorr)
  if (!is.null(outFile)) {
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
                   outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
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
                     outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (attr(x = fitMod, which = "what") == "fixed") {
    stop("Heritability can only be calculated when genotype is random.\n")
  }
  h2 <- lapply(X = fitMod, FUN = SpATS::getHeritability)
  genoDec <- fitMod[[1]]$model$geno$geno.decomp
  h2Out <- data.frame(timePoint = lubridate::as_datetime(names(h2)),
                      row.names = NULL)
  if (!is.null(genoDec)) {
    totDat <- Reduce(f = rbind, x = lapply(fitMod, `[[`, "data"))
    h2Mat <- matrix(nrow = length(h2), ncol = nlevels(totDat[[genoDec]]),
                    dimnames = list(NULL, levels(totDat[[genoDec]])))
    for (i in seq_along(h2)) {
      h2Mat[i, match(names(h2[[i]]),
                     paste0(genoDec, colnames(h2Mat)))] <- h2[[i]]
    }
    h2Out <- cbind(h2Out, h2Mat)
  } else {
    h2Out <- cbind(h2Out, data.frame(h2 = unlist(h2), row.names = NULL))
  }
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
                       outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (!inherits(fitMod[[1]], "SpATS")) {
    stop("Models in ", fitMod, " should be fitted using SpATS.\n")
  }
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

