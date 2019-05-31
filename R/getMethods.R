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
  varCol <- sapply(X = fitMod, FUN = function(x) x$var.comp["colId"])
  varRow <- sapply(X = fitMod, FUN = function(x) x$var.comp["rowId"])
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
    h2Out <- cbind(h2Out, h2)
  }
  if (!is.null(outFile)) {
    write.csv(h2, file = outFile, row.names = FALSE)
  }
  return(h2Out)
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






