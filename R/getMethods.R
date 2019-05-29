#' @export
getBLUPs <- function(fitMod,
                     outFile = NULL) {
  BLUPsTP <- lapply(X = fitMod, FUN = method1)
  BLUPs <- Reduce(f = rbind, x = lapply(X = BLUPsTP, FUN = `[[`, "predGeno"))
  if (!is.null(outFile)) {
    write.csv(BLUPs, file = outFile, row.names = FALSE)
  }
  return(BLUPs)
}

#' @export
getCorrected <- function(fitMod,
                         outFile = NULL) {
  spatCorrTP <- lapply(X = fitMod, FUN = method2)
  spatCorr <- dataCorrSpat <- Reduce(f = rbind, x = spatCorrTP)
  if (!is.null(outFile)) {
    write.csv(spatCorr, file = outFile, row.names = FALSE)
  }
  return(spatCorr)
}

#' @export
getVar <- function(fitMod,
                   outFile = NULL) {
  varRes <- sapply(X = fitMod, FUN = function(x) x$psi[1])
  varCol <- sapply(X = fitMods, FUN = function(x) x$var.comp["colId"])
  varRow <- sapply(X = fitMods, FUN = function(x) x$var.comp["rowId"])
  variance <- data.frame(timePoint = lubridate::as_datetime(names(varRes)),
                         varRes = varRes, varCol = varCol, varRow = varRow,
                         row.names = NULL)
  if (!is.null(outFile)) {
    write.csv(variance, file = outFile, row.names = FALSE)
  }
  return(variance)
}

#' @export
getHerit <- function(fitMod,
                     outFile = NULL) {
  h2 <- sapply(X = fitMod, FUN = SpATS::getHeritability)
  h2 <- data.frame(timePoint = lubridate::as_datetime(names(h2)),
                   h2 = h2, row.names = NULL)
  if (!is.null(outFile)) {
    write.csv(h2, file = outFile, row.names = FALSE)
  }
  return(h2)
}

#' @export
getEffDims <- function(fitMod,
                       outFile = NULL) {
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
    write.csv(effDim, file = outFile, row.names = FALSE)
  }
  return(effDim)
}






