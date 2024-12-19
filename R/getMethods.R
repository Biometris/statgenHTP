#' Extract predicted genotypic values
#'
#' Extract predictions of the genotypic value from an object of class
#' \code{fitMod}.
#'
#' @param fitMod An object of class \code{fitMod}.
#' @param timePoints A character or numeric vector indicating the time point(s)
#' for which the predictions should be extracted. When using a character string
#' to reference a time point, the value has to be an exact match to one of the
#' existing time points. When using a number it will be matched by its number
#' ("timeNumber") in the timePoints attribute of the TP object.
#' @param predictChecks Should predictions of the check genotypes be included
#' in the ouptut. If \code{TRUE} a list of two \code{data.frames} is returned
#' from the function, one with the predictions for the regular genotypes and
#' one with the predictions for the checks.
#' @param outFile A character string indicating the .csv file to which the
#' results should be written. If \code{NULL} no file is written.
#'
#' @returns A list of two data.frames with predicted genotypic values per time
#' point. \code{genoPred} with the predicted values for the genotypes and
#' \code{checkPred} with the predicted values for the checks. If
#' \code{predictChecks = FALSE} the latter will be \code{NULL}.
#'
#' @examples
#' ## Using the first example dataset (PhenovatorDat1).
#' \donttest{
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points.
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 20))
#'
#' ## Extract the genotypic predictions for one time point:
#' genoPredSp <- getGenoPred(modPhenoSp,
#'                           timePoints = 6)
#' head(genoPredSp)
#' }
#'
#' @family functions for spatial modeling
#'
#' @export
getGenoPred <- function(fitMod,
                        timePoints = names(fitMod),
                        predictChecks = FALSE,
                        outFile = NULL) {
  ## Checks.
  if (missing(fitMod) || !inherits(fitMod, "fitMod")) {
    stop(fitMod, " should be an object of class fitMod.\n")
  }
  timePoints <- chkTimePoints(fitMod, timePoints)
  ## Restrict fitMod to selected timePoints.
  fitMod <- fitMod[timePoints]
  ## Get predictions per time point.
  ## predictGeno will throw warnings for every timepoint when no check
  ## was used in the model but predictChecks is TRUE.
  ## Remove the duplicate warnings.
  totPred <- tryCatchExt(lapply(X = fitMod, FUN = predictGeno,
                                predictChecks = predictChecks))
  if (!is.null(totPred$error)) {
    stop(totPred$error)
  } else if (!is.null(totPred$warning)) {
    warning(unique(totPred$warning))
  }
  totPred <- totPred$value
  ## Create one data.frame containing all genotypes for all time points.
  genoPred <- do.call(what = rbind, args = lapply(totPred, `[[`, "predGeno"))
  ## Create one data.frame containing all checks for all time points.
  checkPred <- do.call(what = rbind, args = lapply(totPred, `[[`, "predCheck"))
  ## Create a data.frame with combinations of genotype and geno.decomp.
  ## Only combinations that are present in at least one of the timePoints are
  ## included.
  genoFull <- unique(genoPred[colnames(genoPred) %in%
                                c("genotype", "geno.decomp")])
  genoFull <- merge(unique(genoPred["timePoint"]), genoFull)
  ## Merge to the predictions to get NA predictions for genotypes and checks
  ## that are completely missing at a certain timepoint.
  genoPred <- merge(genoFull, genoPred, all.x = TRUE)
  ## Add time numbers.
  genoPred <- addTimeNumber(fitMod, genoPred)
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "csv")
    write.csv(genoPred, file = outFile, row.names = FALSE)
  }
  if (!is.null(checkPred)) {
    ## Repeat actions for checkPred.
    checkFull <- unique(checkPred[colnames(checkPred) %in%
                                    c("check", "geno.decomp")])
    checkFull <- merge(unique(checkPred["timePoint"]), checkFull)
    checkPred <- merge(checkFull, checkPred, all.x = TRUE)
    checkPred <- addTimeNumber(fitMod, checkPred)
    if (!is.null(outFile)) {
      ## Construct name for outfile.
      outFileCheck <- paste0(substring(outFile, first = 1,
                                       last = nchar(outFile) - 4),
                             "Check.csv")
      ## Check if file exists and is writable.
      chkFile(outFileCheck, fileType = "csv")
      write.csv(checkPred, file = outFileCheck, row.names = FALSE)
    }
  } else {
    checkPred <- NULL
  }
  return(list(genoPred = genoPred, checkPred = checkPred))
}

#' Extract corrected phenotypic values
#'
#' Extract corrected phenotypic values from an object of class fitMod. After
#' fitting a spatial model at each time point, the raw phenotypic data is
#' corrected by subtracting the (estimated) sources of variation (environmental,
#' design  effect) that are of no interest (nuisances). This allows keeping
#' the data resolution at the plot/plant level.
#'
#' @inheritParams getGenoPred
#'
#' @param timePoints A character or numeric vector indicating the time point(s)
#' for which the corrected values should be extracted. When using a character
#' string to reference a time point, the value has to be an exact match to one
#' of the existing time points. When using a number it will be matched by its
#' number ("timeNumber") in the timePoints attribute of the TP object.
#'
#' @returns A data.frame with spatially corrected values per time point.
#'
#' @examples
#' \donttest{
#' ## Using the first example dataset (PhenovatorDat1).
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points.
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 20))
#'
#' ## Extract the corrected values for one time point:
#' spatCorrSp <- getCorrected(modPhenoSp,
#'                            timePoints = 6)
#' head(spatCorrSp)
#' }
#'
#' @family functions for spatial modeling
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
#' @param timePoints A character or numeric vector indicating the time point(s)
#' for which the variances should be extracted. When using a character
#' string to reference a time point, the value has to be an exact match to one
#' of the existing time points. When using a number it will be matched by its
#' number ("timeNumber") in the timePoints attribute of the TP object.
#'
#' @returns A data.frame with variances per time point.
#'
#' @examples
#' \donttest{
#' ## Using the first example dataset (PhenovatorDat1):
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points.
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 20))
#'
#' ## Extract the variances for all available time points.
#' getVar(modPhenoSp)
#' }
#'
#' @family functions for spatial modeling
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
  useCheck <- attr(x = fitMod, which = "useCheck")
  colVarId <- ifelse(useRepId, "repId:colId", "colId")
  rowVarId <- ifelse(useRepId, "repId:rowId", "rowId")
  genoCol <- if (useCheck) "genoCheck" else "genotype"
  if (inherits(fitMod[[1]], "SpATS")) {
    geno.decomp <- fitMod[[1]]$model$geno$geno.decomp
    if (!is.null(geno.decomp)) {
      varGenDf <- lapply(X = fitMod, FUN = function(x) {
        levGD <- levels(x$data[["geno.decomp"]])
        varGD <- x$var.comp[paste0("geno.decomp", levGD)]
        names(varGD) <- paste0("var_geno.decomp_", levGD)
        return(as.data.frame(t(varGD)))
      })
      varGen <- dfBind(varGenDf)
    } else {
      varGen <- sapply(X = fitMod, FUN = function(x) {
        x$var.comp[genoCol]
      })
      varGen <- matrix(varGen, dimnames = list(NULL,"varGen"))
    }
    varRes <- sapply(X = fitMod, FUN = function(x) x$psi[1])
    varCol <- sapply(X = fitMod, FUN = function(x) x$var.comp[colVarId])
    varRow <- sapply(X = fitMod, FUN = function(x) x$var.comp[rowVarId])
  } else if (inherits(fitMod[[1]], "asreml")) {
    ## Get geno.decomp from fitted models.
    if ("geno.decomp" %in% all.vars(fitMod[[1]]$formulae$random) ||
        "geno.decomp" %in% all.vars(fitMod[[1]]$formulae$fixed)) {
      geno.decomp <- "geno.decomp"
    } else {
      geno.decomp <- NULL
    }
    if (!is.null(geno.decomp)) {
      varGenDf <- lapply(X = fitMod, FUN = function(x) {
        levGD <- levels(x$call$data[["geno.decomp"]])
        varGD <- x$vparameters[paste0("at(geno.decomp, ", levGD, "):",
                                      genoCol)] * x$sigma2
        names(varGD) <- paste0("var_geno.decomp_", levGD)
        return(as.data.frame(t(varGD)))
      })
      varGen <- dfBind(varGenDf)
    } else {
      varGen <- sapply(X = fitMod, FUN = function(x) {
        x$vparameters[genoCol] * x$sigma2
      })
      varGen <- matrix(varGen, dimnames = list(NULL,"varGen"))
    }
    varRes <- sapply(X = fitMod, FUN = function(x) x$sigma2)
    varCol <- sapply(X = fitMod, FUN = function(x) {
      x$vparameters[colVarId] * x$sigma2
    })
    varRow <- sapply(X = fitMod, FUN = function(x) {
      x$vparameters[rowVarId] * x$sigma2
    })
  }
  variance <- data.frame(timePoint = lubridate::as_datetime(names(varRes)),
                         varGen, varRes = varRes, varCol = varCol,
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
#' @param timePoints A character or numeric vector indicating the time point(s)
#' for which the heritabilities should be extracted. When using a character
#' string to reference a time point, the value has to be an exact match to one
#' of the existing time points. When using a number it will be matched by its
#' number ("timeNumber") in the timePoints attribute of the TP object.
#'
#' @returns A data.frame with heritabilities per time point.
#'
#' @examples
#' \donttest{
#' ## Using the first example dataset (PhenovatorDat1):
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points.
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 20))
#'
#' ## Extract the heritabilities for all available time points.#'
#' getHerit(modPhenoSp)
#' }
#'
#' @family functions for spatial modeling
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
#' @description Extract effective dimensions from an object of class fitMod.
#' The table below gives an overview of the effective dimensions and an
#' explanation of their meaning.
#'
#' | Effective Dimension | Explanation |
#' | :---------------------------- |:----------------------------------|
#' | colId | Linear trend along columns |
#' | rowId | Linear trend along rows |
#' | fCol | Smooth trend along columns |
#' | fRow | Smooth trend along rows |
#' | fColRow | Linear trend in rows changing smoothly along cols |
#' | colfRow | Linear trend in cols changing smoothly along rows |
#' | fColfRow | Smooth-by-smooth interaction trend over rows and cols |
#' | surface | Sum of smooth trends |
#'
#' @inheritParams getGenoPred
#'
#' @param timePoints A character or numeric vector indicating the time point(s)
#' for which the effective dimension should be extracted. When using a character
#' string to reference a time point, the value has to be an exact match to one
#' of the existing time points. When using a number it will be matched by its
#' number ("timeNumber") in the timePoints attribute of the TP object.
#'
#'
#' @param EDType A character string specifying if the effective dimension
#' ("dimension") or the ratio of effective dimensions ("ratio") should be
#' returned.
#'
#' @returns A data.frame with effective dimensions per time point.
#'
#' @examples
#' \donttest{
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
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#'
#' ## Fit a SpATS model on few time points:
#' modPhenoSp <- fitModels(TP = phenoTP,
#'                         trait = "EffpsII",
#'                         timePoints = c(1, 6, 20))
#'
#' ## Extract the effective dimensions for all available time points in the
#' ## model object:
#' effDimSp <- getEffDims(modPhenoSp)
#' }
#'
#' @family functions for spatial modeling
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

