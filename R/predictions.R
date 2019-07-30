#' @importFrom SpATS predict.SpATS
#' @keywords internal
predictGeno <- function(fitMod) {
  if (inherits(fitMod, "SpATS")) {
    ## Get name of genotype column used.
    genoCol <- fitMod$model$geno$genotype
    ## Check if check was used when fitting model.
    useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))

    geno.decomp <- fitMod$model$geno$geno.decomp
    ## Genotype prediction (including the effect of geno.decomp as well as
    ## the intercept).
    predGeno <- predict(fitMod,
                        which = c(genoCol, geno.decomp,
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
    ## Temporary fix for difference between SpATS and asreml predictions.
    ## asreml predicts marginal means whereas SpATS predicts conditional means.
    ## By adding the means of the fixed effects to the conditional means the
    ## marginal means are calculated.
    ## Note that this means the standard errors are no longer correct.
    corVars <- setdiff(all.vars(fitMod$model$fixed), c(geno.decomp, "check"))
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
  } else if (inherits(fitMod, "asreml")) {
    ## Check if check was used when fitting model.
    useCheck <- grepl(pattern = "check", x = deparse(fitMod$formulae$fixed))
    ## Get name of genotype column used.
    genoCol <- if (useCheck) "genoCheck" else "genotype"
    useGenoDecomp <- "geno.decomp" %in% all.vars(fitMod$formulae$random)
    ## Genotype prediction (including the effect of geno.decomp as well as
    ## the intercept).
    classForm <- paste0(if (useGenoDecomp) "geno.decomp:", genoCol)
    predGeno <- predictAsreml(fitMod, classify = classForm,
                              present = c(genoCol,
                          #                if (useGenoDecomp) "geno.decomp",
                                          if (useCheck) "check"),
                              vcov = FALSE)$pvals
    if (useGenoDecomp) {
      predGeno <- predGeno[startsWith(as.character(predGeno[[genoCol]]),
                                      as.character(predGeno[["geno.decomp"]])), ]
    }
    predGeno <- predGeno[predGeno[["status"]] == "Estimable", ]
    ## Repeat for the check genotypes.
    if (useCheck) {
      ## Predict check genotypes.
      predCheck <- predict(fitMod, classify = "check", vcov = FALSE)$pvals
      ## Rename check to genotype for merging with genotype predictions.
      predCheck[["genotype"]] <- predCheck[["check"]]
      ## Remove noCheck
      predCheck <- predCheck[predCheck[["genotype"]] != "noCheck", ]
      ## Rename genoCol to genotype for merging with check predictions.
      predGeno[["genotype"]] <- predGeno[[genoCol]]
      if (useGenoDecomp) {
        predGeno <- rbind(predGeno[-c(1, 2)], predCheck[-1])
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
  }
  ## Select the variables needed for subsequent analyses.
  predGeno <- predGeno[c("timePoint", "genotype", "predicted.values",
                         "standard.errors")]
  return(predGeno)
}

#' @keywords internal
predictCol <- function(fitMod) {
  if (inherits(fitMod, "SpATS")) {
    ## Col prediction (including intercept).
    predCol <- predict(fitMod, which = "colId")
    ## Include time point.
    predCol[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  } else if (inherits(fitMod, "asreml")) {
    ## Col prediction (including intercept).
    predCol <- predict(fitMod, classify = "colId")$pvals
    ## Rename columns to match those from SpATS predictions.
    colnames(predCol)[colnames(predCol) == "predicted.value"] <-
      "predicted.values"
    colnames(predCol)[colnames(predCol) == "std.error"] <-
      "standard.errors"
    ## Include time point.
    predCol[["timePoint"]] <- fitMod$call$data[["timePoint"]][1]
  }
  ## Select the variables needed for subsequent analyses.
  predCol <- predCol[c("timePoint", "colId", "predicted.values",
                       "standard.errors")]
  return(predCol)
}

#' @keywords internal
predictRow <- function(fitMod) {
  if (inherits(fitMod, "SpATS")) {
    ## Row prediction (including intercept).
    predRow <- predict(fitMod, which = "rowId")
    ## Include time point.
    predRow[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  } else if (inherits(fitMod, "asreml")) {
    ## Col prediction (including intercept).
    predRow <- predict(fitMod, classify = "rowId")$pvals
    ## Rename columns to match those from SpATS predictions.
    colnames(predRow)[colnames(predRow) == "predicted.value"] <-
      "predicted.values"
    colnames(predRow)[colnames(predRow) == "std.error"] <-
      "standard.errors"
    ## Include time point.
    predRow[["timePoint"]] <- fitMod$call$data[["timePoint"]][1]
  }
  ## Select the variables needed for subsequent analyses.
  predRow <- predRow[c("timePoint", "rowId", "predicted.values",
                       "standard.errors")]
  return(predRow)
}
