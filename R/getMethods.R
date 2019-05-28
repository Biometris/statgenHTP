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
getSigma2 <- function(fitMod,
                      outFile = NULL) {
  sigma2 <- sapply(X = fitMod, FUN = function(x) x$psi[1])
  if (!is.null(outFile)) {
    write.csv(sigma2, file = outFile, row.names = FALSE)
  }
  return(sigma2)
}

#' @export
getHerit <- function(fitMod,
                     outFile = NULL) {
  h2 <- sapply(X = fitMods, FUN = SpATS::getHeritability)
  if (!is.null(outFile)) {
    write.csv(h2, file = outFile, row.names = FALSE)
  }
  return(h2)
}

