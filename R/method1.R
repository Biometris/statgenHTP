#' @importFrom SpATS predict.SpATS
method1 <- function(fit.SpATS) {
  ##############################################################################
  # Approach 1: obtain the genotypic predictions
  ##############################################################################

  ##############################
  # Genotype predictions
  ##############################
  # Genotype prediction (including the effect of TrtPop, as well as the intercept)
  predGeno <- predict(fit.SpATS, which = "genotype") # ,"Check"
  ## Include time point
  predGeno[["time"]] <- fit.SpATS$data[["time"]][1]
  # Select the needed variables for subsequent analyses
  predGeno <- predGeno[c("time", "genotype", "predicted.values", "standard.errors")]

  ##############################
  # Col predictions
  ##############################
  # Col prediction (including intercept)
  predCol <- predict(fit.SpATS, which = "colId")
  ## Include time point
  predCol[["time"]] <- fit.SpATS$data[["time"]][1]
  # Select the needed variables for subsequent analyses
  predCol <- predCol[c("time", "colId", "predicted.values", "standard.errors")]

  ##############################
  # Row predictions
  ##############################
  # Row prediction (including intercept)
  predRow <- predict(fit.SpATS, which = "rowId")
  ## Include time point
  predRow[["time"]] <- fit.SpATS$data[["time"]][1]
  # Select the needed variables for subsequent analyses
  predRow <- predRow[c("time", "rowId", "predicted.values", "standard.errors")]

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
                          time = fit.SpATS$data[["time"]][1])

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
                         time = fit.SpATS$data[["time"]][1])

  ##############################
  # Row BLUPs
  ##############################
  # BLUPs
  rows <- paste0("rowId", levels(fit.SpATS$data[["rowId"]]))
  # BLUPs
  BLUPsRow <- fit.SpATS$coeff[rows]
  # Standard error
  seBLUPsRow <- sqrt(diag(fit.SpATS$vcov$C22_inv[rows, rows]))
  # Create data.frame the needed variables for subsequent analyses
  BLUPsRow <- data.frame(rowId = rows, predicted.values = BLUPsRow,
                         standard.errors = seBLUPsRow,
                         time = fit.SpATS$data[["time"]][1])

  return(list(predGeno = predGeno, predCol = predCol, predRow = predRow,
              BLUPsGeno = BLUPsGeno, BLUPsCol = BLUPsCol, BLUPsRow = BLUPsRow))
}
