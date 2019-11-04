#' Extract genotypic predicted values
#'
#' Extract predictions of the genotypic value from an object of class fitMod.
#'
#' @param fitMod An object of class fitMod.
#' @param timePoints A character or numeric vector indicating the time point(s)
#' to be modeled. When using a character string to reference a time point, the
#' value has to be an exact match to one of the existing time points. When using
#' a number it will be matched by its number ("timeNumber") in the timePoints
#' attribute of the TP object.
#' @param outFile A character string indicating the .csv file to which the
#' results should be written. If \code{NULL} no file is written.
#'
#' @return A data.frame with genotypic predicted values per time point.
#'
#' @examples
#' ## Using the first example dataset (PhenovatorDat1):
#' data("PhenovatorDat1")
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2", "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = seq(1,73,by=5))
#'
#' ## Extract the genotypic predictions for one time point:
#' genoPredSp <- getGenoPred(modPhenoSp, timePoints = 6)
#'
#' @export
getGenoPred <- function(fitMod,
                        timePoints = names(fitMod),
                        outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  ## Get predictions per time point.
  genoPred <- lapply(X = fitMod, FUN = predictGeno)
  ## Create one data.frame containing all time points.
  genoPred <- do.call(what = rbind, args = genoPred)
  ## Create a data.frame with combinations of genotype and geno.decomp
  ## Only combinations that are present in at least one of the timePoints are
  ## included.
  full <- unique(genoPred[colnames(genoPred) %in% c("genotype", "geno.decomp")])
  full <- merge(unique(genoPred["timePoint"]), full)
  ## Merge to the predictions to get NA predictions for genotypes that are
  ## completely missing at a certain timepoint.
  genoPred <- merge(full, genoPred, all.x = TRUE)
  ## Add time numbers.
  genoPred <- addTimeNumber(fitMod, genoPred)
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "csv")
    write.csv(genoPred, file = outFile, row.names = FALSE)
  }
  return(genoPred)
}

#' Extract corrected phenotypic values
#'
#' Extract corrected phenotype from an object of class fitMod. After fitting a
#' spatial model at each time point, the raw phenotypic data is corrected by
#' subtracting the (estimated) sources of (environmental, design effect) which
#' are of no interest (nuisances). This allows keeping the data resolution at
#' the plot/plant level.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with spatially corrected values per time point.
#'
#' @examples
#' ## Using the first example dataset (PhenovatorDat1):
#' data("PhenovatorDat1")
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2", "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = seq(1,73,by=5))
#'
#' ## Extract the corrected values for one time point:
#' spatCorrSp <- getCorrected(modPhenoSp, timePoints = 6)
#'
#' @export
getCorrected <- function(fitMod,
                         timePoints = names(fitMod),
                         outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  ## correctSpatial will throw warnings for every timepoint when no spatial
  ## or fixed effects are present. Catch these and only show once.
  spatCorrTP <- tryCatchExt(lapply(X = fitMod, FUN = correctSpatial))
  if (!is.null(spatCorrTP$error)) {
    stop(spatCorrTP$error)
  } else if (!is.null(spatCorrTP$warning)) {
    warning(unique(spatCorrTP$warning))
  }
  ## Create one data.frame with corrected values for all time points.
  spatCorr <- Reduce(f = rbind, x = spatCorrTP$value)
  ## Add time number.
  spatCorr <- addTimeNumber(fitMod, spatCorr)
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "csv")
    write.csv(spatCorr, file = outFile, row.names = FALSE)
  }
  return(spatCorr)
}

#' Extract variances
#'
#' Extract variances from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with variances per time point.
#'
#' @examples
#' ## Using the first example dataset (PhenovatorDat1):
#' data("PhenovatorDat1")
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2", "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = seq(1,73,by=5))
#'
#' ## Extract the variances for all available time points in the model object:
#' varianceSp <- getVar(modPhenoSp)
#'
#' @export
getVar <- function(fitMod,
                   timePoints = names(fitMod),
                   outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  useRepId <- attr(x = fitMod, which = "useRepId")
  colVarId <- ifelse(useRepId, "repId:colId", "colId")
  rowVarId <- ifelse(useRepId, "repId:rowId", "rowId")
  if (inherits(fitMod[[1]], "SpATS")) {
    varGen <- sapply(X = fitMod, FUN = function(x) x$var.comp["genotype"])
    varRes <- sapply(X = fitMod, FUN = function(x) x$psi[1])
    varCol <- sapply(X = fitMod, FUN = function(x) x$var.comp[colVarId])
    varRow <- sapply(X = fitMod, FUN = function(x) x$var.comp[rowVarId])
  } else if (inherits(fitMod[[1]], "asreml")) {
    varGen <- sapply(X = fitMod, FUN = function(x) {
      x$vparameters["genotype"] * x$sigma2
    })
    varRes <- sapply(X = fitMod, FUN = function(x) x$sigma2)
    varCol <- sapply(X = fitMod, FUN = function(x) {
      x$vparameters[colVarId] * x$sigma2
    })
    varRow <- sapply(X = fitMod, FUN = function(x) {
      x$vparameters[rowVarId] * x$sigma2
    })
  }
  variance <- data.frame(timePoint = lubridate::as_datetime(names(varRes)),
                         varGen = varGen, varRes = varRes, varCol = varCol,
                         varRow = varRow, row.names = NULL)
  variance <- addTimeNumber(fitMod, variance)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(variance, file = outFile, row.names = FALSE)
  }
  return(variance)
}

#' Extract heritabilities
#'
#' Extract heritabilities from an object of class fitMod. When
#' \code{geno.decomp} is used, the heritabilities of each level
#' of geno.decomp are stored in separate columns.
#'
#' @inheritParams getGenoPred
#'
#' @return A data.frame with heritabilities per time point.
#'
#' @examples
#' ## Using the first example dataset (PhenovatorDat1):
#' data("PhenovatorDat1")
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2", "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = seq(1,73,by=5))
#'
#' ## Extract the heritabilities for all available time points in the model object:
#' heritSp    <- getHerit(modPhenoSp)
#'
#' @export
getHerit <- function(fitMod,
                     timePoints = names(fitMod),
                     outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  if (attr(x = fitMod, which = "what") == "fixed") {
    stop("Heritability can only be calculated when genotype is random.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  h2Out <- lapply(X = fitMod, FUN = heritability)
  h2Out <- dfBind(h2Out)
  h2Out <- addTimeNumber(fitMod, h2Out)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(h2Out, file = outFile, row.names = FALSE)
  }
  return(h2Out)
}

#' Extract effective dimensions
#'
#' Extract effective dimensions from an object of class fitMod.
#'
#' @inheritParams getGenoPred
#'
#' @param EDType A character string specifying if the effective dimension
#' ("dimension") or the ratio of effective dimensions ("ratio") should be
#' returned.
#'
#' @return A data.frame with effective dimensions per time point.
#'
#' @examples
#' ## Using the first example dataset (PhenovatorDat1):
#' data("PhenovatorDat1")
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2", "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = seq(1,73,by=5))
#'
#' ## Extract the effective dimensions for all available time points in the model object:
#' effDimSp <- getEffDims(modPhenoSp)
#'
#' @export
getEffDims <- function(fitMod,
                       timePoints = names(fitMod),
                       EDType = c("dimension", "ratio"),
                       outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop("fitMod should be an object of class fitMod.\n")
  }
  if (!inherits(fitMod[[1]], "SpATS")) {
    stop("Models in ", deparse(substitute(fitMod)), " should be fitted using",
         " SpATS.\n")
  }
  EDType <- match.arg(EDType)
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  useRepId <- attr(x = fitMod, which = "useRepId")
  colVarId <- ifelse(useRepId, "repId:colId", "colId")
  rowVarId <- ifelse(useRepId, "repId:rowId", "rowId")
  effDimOut <- data.frame(timePoint = lubridate::as_datetime(names(fitMod)),
                          row.names = NULL)
  effDimNames <- c(colVarId, rowVarId, "f(colNum)", "f(rowNum)",
                   "f(colNum):rowNum", "colNum:f(rowNum)","f(colNum):f(rowNum)")
  ## Get effective dimensions for spatial terms.
  effDims <- sapply(X = fitMod, FUN = function(x) {
    x$eff.dim[effDimNames]
  })
  ## Transpose to get dimensions in columns.
  effDims <- t(effDims)
  ## Add effective dimension for surface.
  effDims <- cbind(effDims, rowSums(effDims[, 3:7, drop = FALSE]))
  if (EDType == "ratio") {
    ## Get nominal values for spatial terms.
    effDimNom <- sapply(X = fitMod, FUN = function(x) {
      x$dim.nom[effDimNames]
    })
    ## Transpose to get dimensions in columns.
    effDimNom <- t(effDimNom)
    ## Add nominal for surface.
    effDimNom <- cbind(effDimNom, rowSums(effDimNom[, 3:7, drop = FALSE]))
    effDims <- effDims / effDimNom
  }
  ## Rename columns to more readable format.
  colnames(effDims) <- c(colVarId, rowVarId, "fCol", "fRow", "fColRow",
                         "colfRow", "fColfRow", "surface")
  effDimOut <- cbind(effDimOut, effDims)
  effDimOut <- addTimeNumber(fitMod, effDimOut)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "csv")
    write.csv(effDimOut, file = outFile, row.names = FALSE)
  }
  return(effDimOut)
}

#' Helper function for adding time numbers to data containing a time point
#' column.
#'
#' @noRd
#' @keywords internal
addTimeNumber <- function(fitMod,
                          dat) {
  ## Get data.frame containing time numbers and time points.
  timePoints <- attr(x = fitMod, which = "timePoints")
  ## Covert timePoint column to datetime format for merging.
  timePoints[["timePoint"]] <- lubridate::as_datetime(timePoints[["timePoint"]])
  ## Merge time number to data.
  dat <- merge(timePoints, dat)
  ## Put timeNumber as first column and timePoint as second.
  dat <- dat[c(2, 1, 3:ncol(dat))]
}

