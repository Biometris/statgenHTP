#' @keywords internal
predictGeno <- function(fitMod,
                        predictChecks = FALSE) {
  ## Most steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (inherits(fitMod, "SpATS")) {
    pred <- predictGenoSpATS(fitMod, predictChecks = predictChecks)
  } else if (inherits(fitMod, "asreml")) {
    pred <- predictGenoAsreml(fitMod, predictChecks = predictChecks)
  }
  ## return results.
  return(pred)
}

#' @importFrom SpATS predict.SpATS
#' @keywords internal
predictGenoSpATS <- function(fitMod,
                             predictChecks) {
  ## Get name of genotype column used.
  genoCol <- fitMod$model$geno$genotype
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
  if (predictChecks & !useCheck) {
    warning("check was not used when fitting the model.\n",
            "predictChecks set to FALSE", call. = FALSE)
    predictChecks <- FALSE
  }
  useGenoDecomp <- !is.null(fitMod$model$geno$geno.decomp)
  ## Genotype prediction (including the effect of geno.decomp as well as
  ## the intercept).
  predTot <- predict(fitMod, which = c(genoCol,
                                       if (useGenoDecomp) "geno.decomp",
                                       if (useCheck) "check"),
                     predFixed = "marginal")
  predGeno <- droplevels(predTot[!is.na(predTot[[genoCol]]),
                                 colnames(predTot) != "geno.decomp"])
  if (predictChecks) {
    predCheck <- droplevels(predTot[is.na(predTot[[genoCol]]), ])
  } else {
    predCheck <- NULL
  }
  if (useGenoDecomp) {
    ## Merge geno.decomp to predGeno.
    genoGenoDecomp <- unique(fitMod$data[c(genoCol, "geno.decomp")])
    predGeno <- merge(predGeno, genoGenoDecomp, by = genoCol)
    ## Genotype was converted to an interaction term of genotype and
    ## geno.decomp in the proces of fitting the model. That needs to be
    ## undone to get the genotype back in the output again.
    genoStart <- nchar(as.character(predGeno[["geno.decomp"]])) + 2
    predGeno[[genoCol]] <- as.factor(substring(predGeno[[genoCol]],
                                               first = genoStart))
  }
  ## Rename genoCol to genotype for consistency with models without check.
  predGeno[["genotype"]] <- predGeno[[genoCol]]
  ## Include time point.
  predGeno[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  ## Select the variables needed for subsequent analyses.
  predGeno <- predGeno[c("timePoint", if (useGenoDecomp) "geno.decomp",
                         "genotype", "predicted.values", "standard.errors")]
  if (predictChecks) {
    ## Include time point.
    predCheck[["timePoint"]] <- fitMod$data[["timePoint"]][1]
    ## Select the variables needed for subsequent analyses.
    predCheck <- predCheck[c("timePoint", if (useGenoDecomp) "geno.decomp",
                             "check", "predicted.values", "standard.errors")]
  }
  return(list(predGeno = predGeno, predCheck = predCheck))
}

#' @keywords internal
predictGenoAsreml <- function(fitMod,
                              predictChecks) {
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
                            present = c(genoCol,
                                        if (useGenoDecomp) "geno.decomp",
                                        if (useCheck) "check"),
                            vcov = FALSE)$pvals
  genoGenoDecomp <- unique(fitMod$call$data[c(genoCol,
                                              if (useGenoDecomp) "geno.decomp")])
  predGeno <- merge(predGeno, genoGenoDecomp)
  ## Rename columns to match those from SpATS predictions.
  colnames(predGeno)[colnames(predGeno) == "predicted.value"] <-
    "predicted.values"
  colnames(predGeno)[colnames(predGeno) == "std.error"] <-
    "standard.errors"
  ## Rename genoCol to genotype.
  predGeno[["genotype"]] <- predGeno[[genoCol]]
  ## Include time point.
  predGeno[["timePoint"]] <- fitMod$call$data[["timePoint"]][1]
  ## Select the variables needed for subsequent analyses.
  predGeno <- predGeno[c("timePoint", if (useGenoDecomp) "geno.decomp",
                         "genotype", "predicted.values", "standard.errors")]
  ## Repeat for the check genotypes.
  if (predictChecks) {
    ## Predict check genotypes.
    classFormChk <- paste0(if (useGenoDecomp) "geno.decomp:", "check")
    predCheck <- predictAsreml(fitMod, classify = classFormChk,
                               vcov = FALSE)$pvals
    ## Remove aliased observations
    predCheck <- predCheck[predCheck[["status"]] != "Aliased", ]
    ## Remove noCheck
    predCheck <- predCheck[predCheck[["check"]] != "noCheck", ]
    ## Rename columns to match those from SpATS predictions.
    colnames(predCheck)[colnames(predCheck) == "predicted.value"] <-
      "predicted.values"
    colnames(predCheck)[colnames(predCheck) == "std.error"] <-
      "standard.errors"
    ## Include time point.
    predCheck[["timePoint"]] <- fitMod$call$data[["timePoint"]][1]
    ## Select the variables needed for subsequent analyses.
    predCheck <- predCheck[c("timePoint", if (useGenoDecomp) "geno.decomp",
                           "check", "predicted.values", "standard.errors")]
  } else {
    predCheck <- NULL
  }
  return(list(predGeno = predGeno, predCheck = predCheck))
}
