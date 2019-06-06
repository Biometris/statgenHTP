#' Extract genomic predictions
#'
#' Extract genomic predictions from an object of class fitMod.
#'
#' @param fitMod An object of class fitMod.
#' @param outFile A character string indicting the .csv file to which the
#' results should be written. If \code{NULL} no file is written.
#'
#' @return A data.frame with genomic predictions per time point.
#'
#' @export
getGenoPred <- function(fitMod,
                        outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  genoPred <- lapply(X = fitMod, FUN = predictGeno)
  genoPred <- Reduce(f = rbind, x = genoPred)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(genoPred, file = outFile, row.names = FALSE)
  }
  return(genoPred)
}

#' Extract column predictions
#'
#' Extract column predictions from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with column predictions per time point.
#'
#' @export
getColPred <- function(fitMod,
                       outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  colPred <- lapply(X = fitMod, FUN = predictCol)
  colPred <- Reduce(f = rbind, x = colPred)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(colPred, file = outFile, row.names = FALSE)
  }
  return(colPred)
}

#' Extract row predictions
#'
#' Extract row predictions from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with row predictions per time point.
#'
#' @export
getRowPred <- function(fitMod,
                       outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  rowPred <- lapply(X = fitMod, FUN = predictRow)
  rowPred <- Reduce(f = rbind, x = rowPred)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(rowPred, file = outFile, row.names = FALSE)
  }
  return(rowPred)
}

#' Extract genotypic BLUPs
#'
#' Extract genotypic BLUPs from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with genotypic BLUPs per time point.
#'
#' @export
getBLUPsGeno <- function(fitMod,
                         outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  BLUPsGeno <- lapply(X = fitMod, FUN = BLUPsGeno)
  BLUPsGeno <- Reduce(f = rbind, x = BLUPsGeno)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(BLUPsGeno, file = outFile, row.names = FALSE)
  }
  return(BLUPsGeno)
}

#' Extract column BLUPs
#'
#' Extract column BLUPs from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with column BLUPs per time point.
#'
#' @export
getBLUPsCol <- function(fitMod,
                        outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  BLUPsCol <- lapply(X = fitMod, FUN = BLUPsCol)
  BLUPsCol <- Reduce(f = rbind, x = BLUPsCol)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(BLUPsCol, file = outFile, row.names = FALSE)
  }
  return(BLUPsCol)
}

#' Extract row BLUPs
#'
#' Extract row BLUPs from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with row BLUPs per time point.
#'
#' @export
getBLUPsRow <- function(fitMod,
                        outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  BLUPsRow <- lapply(X = fitMod, FUN = BLUPsRow)
  BLUPsRow <- Reduce(f = rbind, x = BLUPsRow)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(BLUPsRow, file = outFile, row.names = FALSE)
  }
  return(BLUPsRow)
}

#' Extract spatially corrected values
#'
#' Extract spatially corrected values from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with spatially corrected values per time point.
#'
#' @export
getCorrected <- function(fitMod,
                         outFile = NULL) {
  ## Checks.
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  spatCorrTP <- lapply(X = fitMod, FUN = correctSpatial)
  spatCorr <- Reduce(f = rbind, x = spatCorrTP)
  if (!is.null(outFile)) {
    checkFile(outFile)
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
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (inherits(fitMod[[1]], "SpATS")) {
    varRes <- sapply(X = fitMod, FUN = function(x) x$psi[1])
    varCol <- sapply(X = fitMod, FUN = function(x) x$var.comp["colId"])
    varRow <- sapply(X = fitMod, FUN = function(x) x$var.comp["rowId"])
  } else if (inherits(fitMod[[1]], "asreml")) {
    varRes <- sapply(X = fitMod, FUN = function(x) x$sigma2)
    varCol <- sapply(X = fitMod, FUN = function(x) x$vparameters["colId"])
    varRow <- sapply(X = fitMod, FUN = function(x) x$vparameters["rowId"])
  }
  variance <- data.frame(timePoint = lubridate::as_datetime(names(varRes)),
                         varRes = varRes, varCol = varCol, varRow = varRow,
                         row.names = NULL)
  if (!is.null(outFile)) {
    checkFile(outFile)
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
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
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
  if (!is.null(outFile)) {
    checkFile(outFile)
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
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (!inherits(fitMod[[1]], "SpATS")) {
    stop("Models in ", fitMod, " should be fitted using SpATS.\n")
  }
  effDimSurface <- sapply(X = fitMod, FUN = function(x) {
    sum(x$eff.dim[c("f(colNum)", "f(rowNum)", "f(colNum):rowNum",
                    "colNum:f(rowNum)","f(colNum):f(rowNum)")])
  })
  effDimCol <- sapply(X = fitMod, FUN = function(x) x$eff.dim["colId"])
  effDimRow <- sapply(X = fitMod, FUN = function(x) x$eff.dim["rowId"])
  effDim <- data.frame(timePoint = lubridate::as_datetime(names(effDimSurface)),
                       effDimSurface = effDimSurface,
                       effDimCol = effDimCol, effDimRow = effDimRow,
                       row.names = NULL)
  if (!is.null(outFile)) {
    checkFile(outFile)
    write.csv(effDim, file = outFile, row.names = FALSE)
  }
  return(effDim)
}
