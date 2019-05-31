#' @export
getGenoPred <- function(fitMod,
                     outFile = NULL) {
  genoPred <- lapply(X = fitMod, FUN = predictGeno)
  genoPred <- Reduce(f = rbind, x = genoPred)
  if (!is.null(outFile)) {
    write.csv(genoPred, file = outFile, row.names = FALSE)
  }
  return(genoPred)
}

#' @export
getColPred <- function(fitMod,
                       outFile = NULL) {
  colPred <- lapply(X = fitMod, FUN = predictCol)
  colPred <- Reduce(f = rbind, x = colPred)
  if (!is.null(outFile)) {
    write.csv(colPred, file = outFile, row.names = FALSE)
  }
  return(colPred)
}

#' @export
getRowPred <- function(fitMod,
                        outFile = NULL) {
  rowPred <- lapply(X = fitMod, FUN = predictRow)
  rowPred <- Reduce(f = rbind, x = rowPred)
  if (!is.null(outFile)) {
    write.csv(rowPred, file = outFile, row.names = FALSE)
  }
  return(rowPred)
}

#' @export
getBLUPsGeno <- function(fitMod,
                         outFile = NULL) {
  BLUPsGeno <- lapply(X = fitMod, FUN = BLUPsGeno)
  BLUPsGeno <- Reduce(f = rbind, x = BLUPsGeno)
  if (!is.null(outFile)) {
    write.csv(BLUPsGeno, file = outFile, row.names = FALSE)
  }
  return(BLUPsGeno)
}

#' @export
getBLUPsCol <- function(fitMod,
                         outFile = NULL) {
  BLUPsCol <- lapply(X = fitMod, FUN = BLUPsCol)
  BLUPsCol <- Reduce(f = rbind, x = BLUPsCol)
  if (!is.null(outFile)) {
    write.csv(BLUPsCol, file = outFile, row.names = FALSE)
  }
  return(BLUPsCol)
}

#' @export
getBLUPsRow <- function(fitMod,
                         outFile = NULL) {
  BLUPsRow <- lapply(X = fitMod, FUN = BLUPsRow)
  BLUPsRow <- Reduce(f = rbind, x = BLUPsRow)
  if (!is.null(outFile)) {
    write.csv(BLUPsRow, file = outFile, row.names = FALSE)
  }
  return(BLUPsRow)
}

#' @export
getCorrected <- function(fitMod,
                         outFile = NULL) {
  spatCorrTP <- lapply(X = fitMod, FUN = correctSpatial)
  spatCorr <- Reduce(f = rbind, x = spatCorrTP)
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




# # Add results
# genoPred <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "predGeno"))
# colPred <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "predCol"))
# rowPred <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "predRow"))
# genoBLUPs <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "BLUPsGeno"))
# colBLUPs <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "BLUPsCol"))
# rowBLUPs <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "BLUPsRow"))





