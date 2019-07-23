#' @importFrom SpATS predict.SpATS
#' @keywords internal
predictGeno <- function(fitMod) {
  if (inherits(fitMod, "SpATS")) {
    ## Get name of genotype column used.
    genoCol <- fitMod$model$geno$genotype
    ## Check if check was used when fitting model.
    useCheck <- grepl(pattern = "check", x = deparse(fitMod$model$fixed))
    ## Genotype prediction (including the effect of geno.decomp as well as
    ## the intercept).
    predGeno <- predict(fitMod,
                        which = c(genoCol, fitMod$model$geno$geno.decomp))
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
    ## Include time point.
    predGeno[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  } else if (inherits(fitMod, "asreml")) {
    ## Check if check was used when fitting model.
    useCheck <- grepl(pattern = "check", x = deparse(fitMod$formulae$fixed))
    ## Get name of genotype column used.
    genoCol <- if (useCheck) "genoCheck" else "genotype"
    ## Genotype prediction (including the effect of geno.decomp as well as
    ## the intercept).
    predGeno <- predictAsreml(fitMod, classify = "geno.decomp+genoCheck",
                               vcov = FALSE)$pvals

    predGeno <- predGeno[substring(predGeno$genoCheck, 1,1) ==
                           as.character(predGeno$geno.decomp), ]

    ## Repeat for the check genotypes.
    if (useCheck) {
      ## Predict check genotypes.
      predCheck <- predict(fitMod, classify = "geno.decomp+check", vcov = FALSE)$pvals
      ## Rename check to genotype for merging with genotype predictions.
      predCheck[["genotype"]] <- predCheck[["check"]]
      ## Remove noCheck
      predCheck <- predCheck[predCheck[["genotype"]] != "noCheck", ]


      predCheck <- predCheck[substring(predCheck$check, 1,1) ==
                             as.character(predCheck$geno.decomp), ]

      ## Rename genoCol to genotype for merging with check predictions.
      predGeno[["genotype"]] <- predGeno[[genoCol]]
      predGeno <- rbind(predGeno[-2], predCheck[-2])
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
  ## Col prediction (including intercept)
  predCol <- predict(fitMod, which = "colId")
  ## Include time point.
  predCol[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  } else if (inherits(fitMod, "asreml")) {
    ## Col prediction (including intercept)
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
    ## Row prediction (including intercept)
    predRow <- predict(fitMod, which = "rowId")
    ## Include time point.
    predRow[["timePoint"]] <- fitMod$data[["timePoint"]][1]
  } else if (inherits(fitMod, "asreml")) {
    ## Col prediction (including intercept)
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

#' @keywords internal
BLUPsGeno <- function(fitMod) {
  ## BLUPs.
  genotypes <- fitMod$terms$geno$geno_names
  BLUPsGeno <- fitMod$coeff[genotypes]
  ## Standard error.
  seBLUPsGeno <- sqrt(diag(fitMod$vcov$C11_inv[genotypes, genotypes]))
  ## Create data.frame with variables needed for subsequent analyses.
  BLUPsGeno <- data.frame(genotype = genotypes, predicted.values = BLUPsGeno,
                          standard.errors = seBLUPsGeno,
                          timePoint = fitMod$data[["timePoint"]][1],
                          row.names = NULL)
  return(BLUPsGeno)
}

#' @keywords internal
BLUPsCol <- function(fitMod) {
  ## Col BLUPs.
  columns <- paste0("colId", levels(fitMod$data[["colId"]]))
  BLUPsCol <- fitMod$coeff[columns]
  ## Standard error.
  seBLUPsCol <- sqrt(diag(fitMod$vcov$C22_inv[columns, columns]))
  ## Create data.frame variables needed for subsequent analyses.
  BLUPsCol <- data.frame(colId = columns, predicted.values = BLUPsCol,
                         standard.errors = seBLUPsCol,
                         timePoint = fitMod$data[["timePoint"]][1],
                         row.names = NULL)
  return(BLUPsCol)
}

#' @keywords internal
BLUPsRow <- function(fitMod) {
  ## Col BLUPs.
  rows <- paste0("colId", levels(fitMod$data[["rowId"]]))
  BLUPsRow <- fitMod$coeff[rows]
  ## Standard error.
  seBLUPsRow <- sqrt(diag(fitMod$vcov$C22_inv[rows, rows]))
  ## Create data.frame variables needed for subsequent analyses.
  BLUPsRow <- data.frame(colId = rows, predicted.values = BLUPsRow,
                         standard.errors = seBLUPsRow,
                         timePoint = fitMod$data[["timePoint"]][1],
                         row.names = NULL)
  return(BLUPsRow)
}



