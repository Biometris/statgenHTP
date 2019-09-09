#' @keywords internal
predictGeno <- function(fitMod) {
  ## Get engine from fitted model.
  engine <- class(fitMod)
  ## Most steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (engine == "SpATS") {
    pred <- predictGenoSpATS(fitMod)
  } else if (engine == "asreml") {
    pred <- predictGenoAsreml(fitMod)
  }
  ## return results.
  return(pred)
}

#' @importFrom SpATS predict.SpATS
#' @keywords internal
predictGenoSpATS <- function(fitMod) {
  ## Get name of genotype column used.
  genoCol <- fitMod$model$geno$genotype
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
  useGenoDecomp <- !is.null(fitMod$model$geno$geno.decomp)
  ## Genotype prediction (including the effect of geno.decomp as well as
  ## the intercept).
  predGeno <- predict(fitMod, which = c(genoCol,
                                        if (useGenoDecomp) "geno.decomp",
                                        if (useCheck) "check"))
  ## Repeat for the check genotypes.
  if (useCheck) {
    ## Predict check genotypes.
    predCheck <- predict(fitMod, which = "check")
    ## Rename check to genotype for merging with genotype predictions.
    predCheck[["genotype"]] <- predCheck[["check"]]
    ## Remove noCheck
    predCheck <- predCheck[predCheck[["genotype"]] != "noCheck", ]
    ## Rename genoCol to genotype for merging with check predictions.
    predGeno[["genotype"]] <- predGeno[[genoCol]]
    predGeno <- rbind(predGeno, predCheck)
  }
  if (useGenoDecomp) {
    if (!hasName(x = predGeno, name = "geno.decomp")) {
      genoGenoDecomp <- unique(fitMod$data[c(genoCol, "geno.decomp")])
      predGeno <- merge(predGeno, genoGenoDecomp,
                        by.x = "genotype", by.y = genoCol)
    }
    ## Genotype was converted to an interaction term of genotype and
    ## geno.decomp in the proces of fitting the model. That needs to be
    ## undone to get the genotype back in the output again.
    genoStart <- nchar(as.character(predGeno[["geno.decomp"]])) + 2
    predGeno[["genotype"]] <- as.factor(substring(predGeno[["genotype"]],
                                                  first = genoStart))
  }
  ## Temporary fix for difference between SpATS and asreml predictions.
  ## asreml predicts marginal means whereas SpATS predicts conditional means.
  ## By adding the means of the fixed effects to the conditional means the
  ## marginal means are calculated.
  ## Note that this means the standard errors are no longer correct.
  corVars <- setdiff(all.vars(fitMod$model$fixed),
                     c(genoCol, "geno.decomp", "check"))
  if (length(corVars) > 0) {
    ## Order in descreasing order so variables that are substrings of other
    ## variables are treated correctly.
    corVars <- corVars[order(nchar(corVars), decreasing = TRUE)]
    ## Get coefficients for fixed variables.
    coeffs <- fitMod$coeff[!attr(fitMod$coeff, "random")]
    ## Loop over corVars and adjust predicted value by mean of fixed effects
    ## for corVar. Then remove it from coeff so it isn't used again by a
    ## shorter variable, i.e. repId1 and repId
    for (corVar in corVars) {
      corMean <- mean(c(0, coeffs[grepl(corVar, names(coeffs))]))
      predGeno[["predicted.values"]] <-
        predGeno[["predicted.values"]] + corMean
      coeffs <- coeffs[!grepl(corVar, names(coeffs))]
    }
  }
  ## Include time point.
  predGeno[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  ## Select the variables needed for subsequent analyses.
  predGeno <- predGeno[c("timePoint", if (useGenoDecomp) "geno.decomp",
                         "genotype", "predicted.values", "standard.errors")]
  return(predGeno)
}

#' @keywords internal
predictGenoAsreml <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$formulae$fixed))
  ## Get name of genotype column used.
  genoCol <- if (useCheck) "genoCheck" else "genotype"
  useGenoDecomp <- "geno.decomp" %in% all.vars(fitMod$formulae$random) |
    "geno.decomp" %in% all.vars(fitMod$formulae$fixed)
  ## Genotype prediction (including the effect of geno.decomp as well as
  ## the intercept).
  classForm <- paste0(if (useGenoDecomp) "geno.decomp:", genoCol)
  predGeno <- predictAsreml(fitMod, classify = classForm,
                            present = c(genoCol, if (useCheck) "check"),
                            vcov = FALSE)$pvals
  genoGenoDecomp <- unique(fitMod$call$data[c(genoCol,
                                              if (useGenoDecomp) "geno.decomp")])
  predGeno <- merge(predGeno, genoGenoDecomp)
  ## Repeat for the check genotypes.
  if (useCheck) {
    ## Predict check genotypes.
    classFormChk <- paste0(if (useGenoDecomp) "geno.decomp:", "check")
    predCheck <- predictAsreml(fitMod, classify = classFormChk,
                               vcov = FALSE)$pvals
    ## Rename check to genotype for merging with genotype predictions.
    predCheck[["genotype"]] <- predCheck[["check"]]
    ## Remove noCheck
    predCheck <- predCheck[predCheck[["genotype"]] != "noCheck", ]
    ## Rename genoCol to genotype for merging with check predictions.
    predGeno[["genotype"]] <- predGeno[[genoCol]]
    if (useGenoDecomp) {
      predGeno <- rbind(predGeno[-2], predCheck[-2])
    } else {
      predGeno <- rbind(predGeno[-1], predCheck[-1])
    }
  }
  ## Rename columns to match those from SpATS predictions.
  colnames(predGeno)[colnames(predGeno) == "predicted.value"] <-
    "predicted.values"
  colnames(predGeno)[colnames(predGeno) == "std.error"] <-
    "standard.errors"
  ## Include time point.
  predGeno[["timePoint"]] <- fitMod$call$data[["timePoint"]][1]
  ## Select the variables needed for subsequent analyses.
  predGeno <- predGeno[c("timePoint", if (useGenoDecomp) "geno.decomp",
                         "genotype", "predicted.values", "standard.errors")]
  return(predGeno)
}
