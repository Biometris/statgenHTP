#' Fit spatial models per time point
#'
#' Fit spatial models per time point in an object of class TP.
#'
#' @param TP An object of class TP.
#' @param trait A character string indicating the trait used as response
#' variable in the model.
#' @param covariates A character vector indicating the variables used as
#' fixed effects in the model.
#' @param geno.decomp A character vector indicating the variables used to
#' group the genotypes in the model.
#' @param useCheck Should check genotypes be used as an extra factor in the
#' model?
#'
#' @return An object of class fitMod, a list of fitted spatial models.
#'
#' @export
fitModels <- function(TP,
                      trait,
                      covariates = NULL,
                      geno.decomp = NULL,
                      useCheck = FALSE,
                      engine = c("SpATS", "asreml")) {
  ## Checks.
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  if (!all(sapply(X = TP, FUN = hasName, name = trait))) {
    stop(trait, " should be a column in TP for all timePoints.\n")
  }
  if (!is.null(covariates)) {
    for (covar in covariates) {
      if (!all(sapply(X = TP, FUN = hasName, name = covar))) {
        stop(covar, " should be a column in TP for all timePoints.\n")
      }
    }
  }
  if (!is.null(geno.decomp)) {
    for (gd in geno.decomp) {
      if (!all(sapply(X = TP, FUN = hasName, name = gd))) {
        stop(gd, " should be a column for all timePoints.\n")
      }
    }
  }
  if (useCheck) {
    if (!all(sapply(X = TP, FUN = hasName, name = "check"))) {
      stop("check should be a column in TP for all timePoints.\n")
    }
  }
  engine <- match.arg(engine)
  ## If geno.decomp is used genotype and covariates have to be replaced by
  ## an interaction of genotype and covariates with the geno.decomp variables.
  ## Construct an interaction of all variables in geno.decomp.
  if (length(geno.decomp) > 1) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      timePoint[["geno.decomp"]] <- interaction(timePoint[geno.decomp],
                                                sep = "_")
      return(timePoint)
    })
    ## Set geno.decomp to newly constructed variable.
    geno.decomp <- "geno.decomp"
  }
  ## Replace genotype and covariates by their interaction with geno.decomp.
  if (!is.null(geno.decomp)) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      timePoint[["genotype"]] <- interaction(timePoint[[geno.decomp]],
                                             timePoint[["genotype"]], sep = "_")
      for (covar in covariates) {
        timePoint[[covar]] <- interaction(timePoint[[geno.decomp]],
                                          timePoint[[covar]], sep = "_")
      }
      return(timePoint)
    })
  }
  ## All covariates should be factors. If not convert them to factor.
  for (covar in covariates) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      if (!is.factor(timePoint[[covar]])) {
        timePoint[[covar]] <- as.factor(timePoint[[covar]])
      }
      return(timePoint)
    })
  }
  ## Fixed part consists of covariates, geno.decomp and check.
  if (!is.null(c(covariates, geno.decomp))) {
    fixedForm <- formula(paste("~", paste(c(covariates, geno.decomp),
                                          collapse = "+"),
                               if (useCheck) "+ check"))
  } else {
    fixedForm <- if (useCheck) formula("~ check") else NULL
  }
  if (engine == "SpATS") {
    ## Loop on timepoint to run SpATS.
    fitMods <- lapply(X = TP, function(timePoint) {
      message(timePoint[["timePoint"]][1])
      ## Only keep columns needed for analysis.
      modCols <- c("timePoint", "plotId", "genotype", "genoCheck", "check",
                   "colId", "rowId", "colNum", "rowNum", covariates,
                   geno.decomp, trait)
      modDat <- timePoint[colnames(timePoint) %in% modCols]
      modDat <- droplevels(modDat)
      ## number of segments for SpATS.
      nseg = c(nlevels(modDat[["colId"]]), nlevels(modDat[["rowId"]]))
      ## Fit and return the model.
      SpATS::SpATS(response = trait, fixed = fixedForm,
                   random = ~ colId + rowId,
                   spatial = ~ SpATS::PSANOVA(colNum, rowNum, nseg = nseg,
                                              nest.div = c(2, 2)),
                   genotype = if (useCheck) "genoCheck" else "genotype",
                   genotype.as.random = TRUE, geno.decomp = geno.decomp,
                   data = modDat,
                   control = list(maxit = 50, tolerance = 1e-03,
                                  monitoring = 0))
    })
  } else if (engine == "asreml") {
    ## Loop on timepoint to run asreml.
    fitMods <- lapply(X = TP, function(timePoint) {
      message(timePoint[["timePoint"]][1])
      ## Only keep columns needed for analysis.
      modDat <- droplevels(timePoint)
      ## Add empty observations.
      TDTab <- as.data.frame(table(modDat$colId, modDat$rowId))
      TDTab <- TDTab[TDTab$Freq == 0, , drop = FALSE]
      if (nrow(TDTab) > 0) {
        extObs <- setNames(as.data.frame(matrix(nrow = nrow(TDTab),
                                                ncol = ncol(modDat))),
                           colnames(modDat))
        extObs$trial <- modDat$trial[1]
        extObs[, c("colId", "rowId")] <- TDTab[, c("Var1", "Var2")]
        extObs[, c("colNum", "rowNum")] <-
          c(as.numeric(levels(TDTab[, "Var1"]))[TDTab[, "Var1"]],
            as.numeric(levels(TDTab[, "Var2"]))[TDTab[, "Var2"]])
        modDat <- rbind(modDat, extObs)
      }
      ## TP needs to be sorted by row and column to prevent asreml from crashing.
      modDat <- modDat[order(modDat$rowId, modDat$colId), ]
      ## Run model.
      ## na.method for x set to include to deal with missing rows/columns.
      asreml::asreml(fixed = formula(paste(trait, "~ ", geno.decomp)),
                     random = formula(paste("~ rowId + colId + at(",
                                            geno.decomp, "):genotype")),
                     residual = ~ar1(rowId):ar1(colId),
                     data = modDat, trace = FALSE, maxiter = 200,
                     na.action = asreml::na.method(x = "include"))
    })
  }
  return(createFitMod(fitMods))
}



