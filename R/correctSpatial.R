#' Correct for spatial effects and other unnecessary factors.
#'
#' @noRd
#' @keywords internal
correctSpatial <- function(fitMod) {
  ## All steps are different for SpATS and asreml so just move them to
  ## separate functions.
  if (inherits(fitMod, "SpATS")) {
    pred <- correctSpatialSpATS(fitMod)
  } else if (inherits(fitMod, "asreml")) {
    pred <- correctSpatialAsreml(fitMod)
  }
  ## return results.
  return(pred)
}

#' @noRd
#' @keywords internal
correctSpatialSpATS <- function(fitMod) {
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
  ## Get trait from fitted model.
  trait <- fitMod$model$response
  ## Set name for new trait.
  newTrait <- paste0(trait, "_corr")
  ## Get geno.decomp from fitted models.
  geno.decomp <- fitMod$model$geno$geno.decomp
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing.
  fixVars <- attr(terms(fitMod$model$fixed), "term.labels")
  ## Subset on order = 1 to remove interaction terms.
  fixVars <- fixVars[attr(terms(fitMod$model$fixed), "order") == 1]
  fixVars <- setdiff(fixVars, c("genotype", "genoCheck", "check", geno.decomp))

  randVars <- attr(terms(fitMod$model$random), "term.labels")
  ## Subset on order = 1 to remove interaction terms.
  randVars <- randVars[attr(terms(fitMod$model$random), "order") == 1]
  randVars <- setdiff(randVars, c("genotype", "genoCheck", geno.decomp,
                                  "rowNum", "colNum"))
  ## Get name of genotype column used.
  genoCol <- fitMod$model$geno$genotype
  useGenoDecomp <- !is.null(fitMod$model$geno$geno.decomp)
  ## Genotype prediction (including the effect of geno.decomp as well as
  ## the intercept).
  predVars <- c(genoCol,
                if (useGenoDecomp) "geno.decomp",
                if (useCheck) "check")
  pred <- predict(fitMod, which = predVars, predFixed = "marginal",
                  return.vcov.matrix = TRUE)
  ## Compute Weights based on the inverse of the var-cov (vcov) matrix.
  ## Get the vcov matrix from pred.
  vcov <- attr(x = pred, which = "vcov")
  ## Add the residual error to the diagonal of the vcov matrix
  ## Required since the residuals are added to the predictions.
  vcovComb <- as.matrix(vcov) + fitMod$psi[1] * diag(nrow(vcov))
  pred[["wt"]] <- diag(solve(vcovComb))
  ## Remove redundant columns since these are added from the data
  ## used for fitting the model.
  pred <- pred[, !colnames(pred) %in% c(fixVars, randVars)]
  ## Merge genotype and timepoint to data.
  modDat <- fitMod$data[c("genotype", "plotId", "timePoint", trait,
                          geno.decomp, fixVars, randVars,
                          if (useCheck) c("check", "genoCheck"))]
  modDat[["resid"]] <- fitMod$residuals
  pred <- merge(modDat, pred, by = predVars, all.x = TRUE)
  if (!is.null(geno.decomp) && !useCheck) {
    ## Genotype was converted to an interaction term of genotype and
    ## geno.decomp in the process of fitting the model. That needs to be
    ## undone to get the genotype back in the output again.
    genoStart <- nchar(as.character(pred[["geno.decomp"]])) + 2
    pred[["genotype"]] <- as.factor(substring(pred[["genotype"]],
                                                first = genoStart))
  }
  ## Obtain the corrected trait.
  pred[[newTrait]] <- pred[["predicted.values"]] + pred[["resid"]]
  ## Select the variables needed for subsequent analyses.
  pred <- pred[c(newTrait, trait, "wt", "genotype", if (useCheck) "check",
                 geno.decomp, fixVars, randVars, "plotId", "timePoint")]
}

#' @noRd
#' @keywords internal
correctSpatialAsreml <- function(fitMod) {
  ## Get trait from fitted model.
  trait <- all.vars(update(fitMod$formulae$fixed, .~0))
  ## Set name for new trait.
  newTrait <- paste0(trait, "_corr")
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
                                        if (useCheck) "check"))
  genoGenoDecomp <- unique(fitMod$call$data[c(genoCol,
                                              if (useGenoDecomp) "geno.decomp")])
  vcovGeno <- predGeno$vcov
  vcovGeno <- vcovGeno[predGeno$pvals[["status"]] == "Estimable",
                       predGeno$pvals[["status"]] == "Estimable"]
  ## Add the residual error to the diagonal of the vcov matrix
  ## Required since the residuals are added to the predictions.
  vcovComb <- as.matrix(vcovGeno) + fitMod$sigma2 * diag(nrow(vcovGeno))
  predGeno <- predGeno$pvals[predGeno$pvals[["status"]] == "Estimable", ]
  predGeno[["wt"]] <- diag(solve(vcovComb))
  ## Rename genoCol to genotype.
  colnames(predGeno)[colnames(predGeno) == genoCol] <- "genotype"
  ## Repeat for the check genotypes.
  if (useCheck) {
    ## Predict check genotypes.
    classFormChk <- paste0(if (useGenoDecomp) "geno.decomp:", "check")
    predCheck <- predictAsreml(fitMod, classify = classFormChk)
    vcovCheck <- predCheck$vcov
    vcovCheck <- vcovCheck[predCheck$pvals[["status"]] == "Estimable",
                         predCheck$pvals[["status"]] == "Estimable"]
    ## Add the residual error to the diagonal of the vcov matrix
    ## Required since the residuals are added to the predictions.
    vcovComb <- as.matrix(vcovCheck) + fitMod$sigma2 * diag(nrow(vcovCheck))
    predCheck <- predCheck$pvals
    ## Remove aliased observations
    predCheck <- predCheck[predCheck[["status"]] != "Aliased", ]
    ## Remove noCheck
    predCheck[["wt"]] <- sqrt(diag(solve(vcovCheck)))
    predCheck <- predCheck[predCheck[["check"]] != "noCheck", ]
    ## Rename genoCol to genotype.
    predCheck[["genotype"]] <- predCheck[["check"]]
    ## Add check column to predGeno.
    predGeno[["check"]] <- "noCheck"
  } else {
    predCheck <- NULL
  }
  predTot <- rbind(predGeno, predCheck)
  ## Include in the prediction the factors (variables) whose effect we are
  ## interested in removing.
  fixVars <- all.vars(update(fitMod$formulae$fixed, 0~.))
  fixVars <- setdiff(fixVars, c("genotype", "genoCheck",
                                if (useGenoDecomp) "geno.decomp"))
  randVars <- all.vars(fitMod$formulae$random)
  randVars <- setdiff(randVars, c("genotype", "genoCheck", "units",
                                  if (useGenoDecomp) "geno.decomp"))
  predVars <- c(fixVars, randVars)
  modDat <- fitMod$call$data[union(c("genotype", if (useCheck) "check",
                                     "plotId", "timePoint", trait,
                                     if (useGenoDecomp) "geno.decomp"),
                                   predVars)]
  ## Extract residuals.
  modDat[["resid"]] <- residuals(fitMod)
  ## Construct residuals from 'full' residuals and residuals with nugget.
  # res1 <- residuals(fitMod)
  # res2 <- residuals(fitMod, spatial = "plot") # include nugget
  # res <- res2 - res1
  # modDat[["resid"]] <- res
  predTot <- merge(modDat, predTot, all.x = TRUE)
  predTot[[newTrait]] <- predTot[["predicted.value"]] + predTot[["resid"]]
  if (length(predVars) == 0) {
    warning("No spatial or fixed effects to correct for. Returning raw data.\n")
  }
  ## Remove row/col combinations added when fitting models.
  predTot <- predTot[!is.na(predTot[["plotId"]]), ]
  ## Select the variables needed for subsequent analyses.
  predTot <- predTot[c(newTrait, trait, "wt", "genotype",
                       if (useGenoDecomp) "geno.decomp", fixVars, randVars,
                       "plotId", "timePoint")]
}
