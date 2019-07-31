#' @noRd
#' @keywords internal
heritability <- function(fitMod) {
  ## Get engine from fitted model.
  engine <- class(fitMod)
  ## All steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (engine == "SpATS") {
    h2out <- heritabilitySpATS(fitMod)
  } else {
    h2out <- heritabilityAsreml(fitMod)
  }
  ## return results.
  return(h2out)
}

#' @noRd
#' @keywords internal
heritabilitySpATS <- function(fitMod) {
  h2 <- SpATS::getHeritability(fitMod)
  genoDec <- fitMod[[1]]$model$geno$geno.decomp
  h2Out <- data.frame(timePoint = fitMod$data[["timePoint"]][1],
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
  return(h2Out)
}

#' Code copied and modified from
#' https://github.com/PaulSchmidtGit/Heritability/tree/master/Alternative%20Heritability%20Measures
#'
#' @noRd
#' @keywords internal
heritabilityAsreml <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$formulae$fixed))
  ## Get name of genotype column used.
  genoCol <- if (useCheck) "genoCheck" else "genotype"
  ## Get genetic variance.
  varGen <- fitMod$vparameters[genoCol] * fitMod$sigma2
  ## Obtain squared s.e.d. matrix.
  vdBLUP <- predict(fitMod, classify = genoCol, only = genoCol,
                    sed = TRUE)$sed^2
  ## Compute mean variance of a difference of two genotypic BLUPs.
  vdBLUPMean <- mean(as.numeric(vdBLUP), na.rm = TRUE)
  ## Compute heritability.
  h2Out <- data.frame(timePoint = fitMod$call$data[["timePoint"]][1],
                      h2 = 1 - (vdBLUPMean / 2 / varGen),
                      row.names = NULL)
  return(h2Out)
}


