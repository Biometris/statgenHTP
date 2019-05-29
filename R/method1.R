#' @importFrom SpATS predict.SpATS
method1 <- function(fit.SpATS) {
  ##############################################################################
  # Approach 1: obtain the genotypic predictions
  ##############################################################################

  ## Get name of genotype column used.
  genoCol <- fit.SpATS$model$geno$genotype
  ## Check if check was used when fitting model.
  useCheck <- grepl(pattern = "check", x = deparse(fit.SpATS$model$fixed))
  checkCol <- if (useCheck) "check" else NULL
  ##############################
  # Genotype predictions
  ##############################
  ## Genotype prediction (including the effect of TrtPop, as well as the intercept)
  predGeno <- predict(fit.SpATS,
                      which = c(genoCol, fit.SpATS$model$geno$geno.decomp))
  ## Repeat for the check genotypes.
  if (useCheck) {
    ## Predict check genotypes.
    predCheck <- predict(fit.SpATS, which = checkCol)
    ## Rename check to genotype for merging with genotype predictions.
    predCheck[["genotype"]] <- predCheck[[checkCol]]
    ## Remove noCheck
    predCheck <- predCheck[predCheck[["genotype"]] != "noCheck", ]
    ## Rename genoCol to genotype for merging with check predictions.
    predGeno[["genotype"]] <- predGeno[[genoCol]]
    predGeno <- rbind(predGeno, predCheck)
  }
  ## Include time point.
  predGeno[["timePoint"]] <- fit.SpATS$data[["timePoint"]][1]
  ## Rename genoCol to genotype.
  ## Select the variables needed for subsequent analyses.
  predGeno <- predGeno[c("timePoint", "genotype", "predicted.values",
                         "standard.errors")]
  ##############################
  # Col predictions
  ##############################
  # Col prediction (including intercept)
  predCol <- predict(fit.SpATS, which = "colId")
  ## Include time point
  predCol[["timePoint"]] <- fit.SpATS$data[["timePoint"]][1]
  # Select the needed variables for subsequent analyses
  predCol <- predCol[c("timePoint", "colId", "predicted.values",
                       "standard.errors")]
  ##############################
  # Row predictions
  ##############################
  # Row prediction (including intercept)
  predRow <- predict(fit.SpATS, which = "rowId")
  ## Include time point
  predRow[["timePoint"]] <- fit.SpATS$data[["timePoint"]][1]
  # Select the needed variables for subsequent analyses
  predRow <- predRow[c("timePoint", "rowId", "predicted.values",
                       "standard.errors")]
  ##############################
  # Genotype BLUPs
  ##############################
  # BLUPs
  genotypes <- fit.SpATS$terms$geno$geno_names
  BLUPsGeno <- fit.SpATS$coeff[genotypes]
  # Standard error
  seBLUPsGeno <- sqrt(diag(fit.SpATS$vcov$C11_inv[genotypes, genotypes]))

  # Create data.frame the needed variables for subsequent analyses
  BLUPsGeno <- data.frame(genotype = genotypes, predicted.values = BLUPsGeno,
                          standard.errors = seBLUPsGeno,
                          timePoint = fit.SpATS$data[["timePoint"]][1])
  ##############################
  # Col BLUPs
  ##############################
  columns <- paste0("colId", levels(fit.SpATS$data[["colId"]]))
  # BLUPs
  BLUPsCol <- fit.SpATS$coeff[columns]
  # Standard error
  seBLUPsCol <- sqrt(diag(fit.SpATS$vcov$C22_inv[columns, columns]))
  # Create data.frame the needed variables for subsequent analyses
  BLUPsCol <- data.frame(colId = columns, predicted.values = BLUPsCol,
                         standard.errors = seBLUPsCol,
                         timePoint = fit.SpATS$data[["timePoint"]][1])
  ##############################
  # Row BLUPs
  ##############################
  rows <- paste0("rowId", levels(fit.SpATS$data[["rowId"]]))
  # BLUPs
  BLUPsRow <- fit.SpATS$coeff[rows]
  # Standard error
  seBLUPsRow <- sqrt(diag(fit.SpATS$vcov$C22_inv[rows, rows]))
  # Create data.frame the needed variables for subsequent analyses
  BLUPsRow <- data.frame(rowId = rows, predicted.values = BLUPsRow,
                         standard.errors = seBLUPsRow,
                         timePoint = fit.SpATS$data[["timePoint"]][1])
  return(list(predGeno = predGeno, predCol = predCol, predRow = predRow,
              BLUPsGeno = BLUPsGeno, BLUPsCol = BLUPsCol, BLUPsRow = BLUPsRow))
}
