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
  if (!inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  genoPred <- lapply(X = fitMod, FUN = predictGeno)
  genoPred <- Reduce(f = rbind, x = genoPred)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(genoPred, file = outFile, row.names = FALSE)
  }
  return(genoPred)
}

#' Extract column predictions
#'
#' Extract predictions of the column values from an object of class fitMod.
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
  ## Column prediction can always be done for SpATS and
  ## for asreml if colId was used as random factor.
  if (inherits(fitMod[[1]], "SpATS") ||
      (inherits(fitMod[[1]], "asreml") &&
       any(grepl(pattern = "colId", x = attr(x = fitMod[[1]]$formulae$random,
                                             which = "term.labels"))))) {
    colPred <- lapply(X = fitMod, FUN = predictCol)
    colPred <- Reduce(f = rbind, x = colPred)
    if (!is.null(outFile)) {
      chkFile(outFile, fileType = "csv")
      write.csv(colPred, file = outFile, row.names = FALSE)
    }
    return(colPred)
  } else {
    stop("Models in ", deparse(substitute(fitMod)), " should either be fitted",
         " using SpATS or using asreml with option spatial set to TRUE.\n")
  }
}

#' Extract row predictions
#'
#' Extract predictions of the row values from an object of class fitMod.
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
  ## Row prediction can always be done for SpATS and
  ## for asreml if rowId was used as random factor.
  if (inherits(fitMod[[1]], "SpATS") ||
      (inherits(fitMod[[1]], "asreml") &&
       any(grepl(pattern = "rowId", x = attr(x = fitMod[[1]]$formulae$random,
                                             which = "term.labels"))))) {
    rowPred <- lapply(X = fitMod, FUN = predictRow)
    rowPred <- Reduce(f = rbind, x = rowPred)
    if (!is.null(outFile)) {
      chkFile(outFile, fileType = "csv")
      write.csv(rowPred, file = outFile, row.names = FALSE)
    }
    return(rowPred)
  } else {
    stop("Models in ", deparse(substitute(fitMod)), " should either be fitted",
         " using SpATS or using asreml with option spatial set to TRUE.\n")
  }
}

#' Extract genotypic BLUPs
#'
#' Extract genotypic Best Linear Unbiased Predictors (BLUPs) from an object of class fitMod.
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
    chkFile(outFile, fileType = "csv")
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
    chkFile(outFile, fileType = "csv")
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
    chkFile(outFile, fileType = "csv")
    write.csv(BLUPsRow, file = outFile, row.names = FALSE)
  }
  return(BLUPsRow)
}

#' Extract corrected phenotypic values
#'
#' Extract corrected phenotype from an object of class fitMod. After fitting a spatial model
#' at each time point, the raw phenotypic data is corrected by subtracting the (estimated)
#' sources of (environmental, design effect) which are of no interest (nuisances).
#' This allows keeping the data resolution at the plot/lant level.
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
    chkFile(outFile, fileType = "csv")
    write.csv(effDim, file = outFile, row.names = FALSE)
  }
  return(effDim)
}
