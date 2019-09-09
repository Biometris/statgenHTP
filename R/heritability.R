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
  ## SpATS heritabily function returns a numeric value if geno.decomp was
  ## not used and a vector otherwise.
  h2 <- SpATS::getHeritability(fitMod)
  ## Get variable used as geno.decomp.
  geno.decomp <- fitMod$model$geno$geno.decomp
  ## Create a base data.frame to which the heritability can be merged.
  h2Out <- data.frame(timePoint = fitMod$data[["timePoint"]][1],
                      row.names = NULL)
  if (!is.null(geno.decomp)) {
    names(h2) <- gsub(pattern = "geno.decomp", replacement = "", x = names(h2))
    h2Mat <- matrix(h2, nrow = 1, dimnames = list(NULL, names(h2)))
    h2Out <- cbind(h2Out, h2Mat)
  } else {
    h2Out <- cbind(h2Out, data.frame(h2 = h2, row.names = NULL))
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
  ## Create a base data.frame to which the heritability can be merged.
  h2Out <- data.frame(timePoint = fitMod$call$data[["timePoint"]][1],
                      row.names = NULL)
  if (!"geno.decomp" %in% all.vars(fitMod$formulae$random)) {
    ## Get genetic variance.
    varGen <- fitMod$vparameters[genoCol] * fitMod$sigma2
    ## Obtain squared s.e.d. matrix.
    vdBLUP <- predictAsreml(fitMod, classify = genoCol, vcov = FALSE,
                            only = genoCol, sed = TRUE)$sed^2
    ## Compute mean variance of a difference of two genotypic BLUPs.
    vdBLUPMean <- mean(as.numeric(vdBLUP), na.rm = TRUE)
    ## Compute heritability.
    h2Out <- cbind(h2Out, data.frame(h2 = 1 - (vdBLUPMean / 2 / varGen),
                                     row.names = NULL))
  } else {
    decompLabs <- levels(droplevels(fitMod$call$data[["geno.decomp"]]))
    h2 <- sapply(X = seq_along(decompLabs), FUN = function(i) {
      varName <- paste0("at(geno.decomp, ", decompLabs[i], "):", genoCol)
      ## Get genetic variance.
      levVarGen <- fitMod$vparameters[varName] * fitMod$sigma2
      ## Get predictions for genotype on current level of geno.decomp.
      levPred <- predictAsreml(fitMod,
                               classify = paste0("geno.decomp:", genoCol),
                               vcov = FALSE,
                               present = c("geno.decomp", "genotype"),
                               levels = list(geno.decomp = i), sed = TRUE)
      ## Compute mean variance of a difference of two genotypic BLUPs.
      levVdBLUP <- levPred$sed^2
      ## Only use values that were actually estimated and not aliased.
      levVdBLUP <- levVdBLUP[levPred$pvals$status == "Estimable",
                             levPred$pvals$status == "Estimable"]
      ## Compute mean variance of a difference of two genotypic BLUPs.
      levVdBLUPMean <- mean(as.numeric(levVdBLUP), na.rm = TRUE)
      ## Compute heritability.
      1 - (levVdBLUPMean / 2 / levVarGen)
    })
    h2Mat <- matrix(h2, nrow = 1, dimnames = list(NULL, decompLabs))
    h2Out <- cbind(h2Out, h2Mat)
  }
  return(h2Out)
}
