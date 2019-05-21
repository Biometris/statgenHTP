method1 <- function(fit.SpATS) {
  ##############################################################################
  # Approach 1: obtain the genotypic predictions
  ##############################################################################

  ##############################
  # Genotype predictions
  ##############################
  # Genotype prediction (including the effect of TrtPop, as well as the intercept)
  predGeno <- predict(fit.SpATS, which = c("Geno")) # ,"Check"
  ## Include time point
  predGeno[["Time"]] <- fit.SpATS$data[["Time"]][1]
  # Select the needed variables for subsequent analyses
  predGeno <- predGeno[c("Time", "Geno", "predicted.values", "standard.errors")]

  ##############################
  # Col predictions
  ##############################
  # Col prediction (including intercept)
  predCol <- predict(fit.SpATS, which = "Col")
  ## Include time point
  predCol[["Time"]] <- fit.SpATS$data[["Time"]][1]
  # Select the needed variables for subsequent analyses
  predCol <- predCol[c("Time", "Col", "predicted.values", "standard.errors")]

  ##############################
  # Row predictions
  ##############################
  # Row prediction (including intercept)
  predRow <- predict(fit.SpATS, which = "Row")
  ## Include time point
  predRow[["Time"]] <- fit.SpATS$data[["Time"]][1]
  # Select the needed variables for subsequent analyses
  predRow <- predRow[c("Time", "Row", "predicted.values", "standard.errors")]

  ##############################
  # Genotype BLUPs
  ##############################
  # BLUPs
  genotypes <- fit.SpATS$terms$geno$geno_names
  BLUPsGeno <- fit.SpATS$coeff[genotypes]
  # Standard error
  seBLUPsGeno <- sqrt(diag(fit.SpATS$vcov$C11_inv[genotypes, genotypes]))

  # Create data.frame the needed variables for subsequent analyses
  BLUPsGeno <- data.frame(TrtGeno = names(BLUPsGeno),
                          predicted.values = BLUPsGeno,
                          standard.errors = seBLUPsGeno,
                          Time = fit.SpATS$data[["Time"]][1])

  ##############################
  # Col BLUPs
  ##############################

  columns <- paste0("Col", levels(fit.SpATS$data[["Col"]]))
  # BLUPs
  BLUPsCol <- fit.SpATS$coeff[columns]
  # Standard error
  seBLUPsCol <- sqrt(diag(fit.SpATS$vcov$C22_inv[columns, columns]))
  # Create data.frame the needed variables for subsequent analyses
  BLUPsCol <- data.frame(Col = names(BLUPsCol),
                         predicted.values = BLUPsCol,
                         standard.errors = seBLUPsCol,
                         Time = fit.SpATS$data[["Time"]][1])

  ##############################
  # Row BLUPs
  ##############################
  # BLUPs
  rows <- paste0("Row", levels(fit.SpATS$data[["Row"]]))
  # BLUPs
  BLUPsRow <- fit.SpATS$coeff[rows]
  # Standard error
  seBLUPsRow <- sqrt(diag(fit.SpATS$vcov$C22_inv[rows, rows]))
  # Create data.frame the needed variables for subsequent analyses
  BLUPsRow <- data.frame(Row = names(BLUPsRow),
                         predicted.values = BLUPsRow,
                         standard.errors = seBLUPsRow,
                         Time = fit.SpATS$data[["Time"]][1])

  return(list(predGeno = predGeno, predCol = predCol, predRow = predRow,
              BLUPsGeno = BLUPsGeno, BLUPsCol = BLUPsCol, BLUPsRow = BLUPsRow))
}
