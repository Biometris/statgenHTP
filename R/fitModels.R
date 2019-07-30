#' Fit spatial models per time point
#'
#' Perform REML analysis given a specific experimental design using either
#' SpATS or asreml. SpATS is used as a default method. See details for the exact
#' models fitted.
#'
#' The actual model fitted depends on the function parameters secified. The
#' basic model is the following:\cr
#' trait = \strong{genotype} + e\cr
#' In case \code{useCheck = TRUE}, instead of genotype genoCheck is used and
#' and check is used as an extra fixed effect. So then the model becomes:\cr
#' trait = \emph{check} + \strong{genoCheck} + e\cr
#' Variables in \code{covariates} are fitted as extra fixed effects.\cr\cr
#' When \code{SpATS} is used for modeling, an extra spatial term is included
#' in the model. This term is constructed using the function
#' \code{\link[SpATS]{PSANOVA}} from the SpATS package as\cr
#' \code{PSANOVA(colNum, rowNum, nseg = nSeg, nest.div = 2)}
#' where\cr \code{nSeg = (number of columns, number of rows)}.\cr\cr
#' When \code{asreml} is used for modeling and \code{spatial = TRUE}
#' four models are fitted with different random terms and covariance structure.
#' The best model is determined based on a goodness-of-fit criterion, either
#' AIC or BIC. This can be set using the control parameter \code{criterion},
#' default is AIC.
#' The following combinations of random and spatial terms are fitted
#' \itemize{
#' \item{random = repId:rowId, spatial = NULL}
#' \item{random = repId:rowId, spatial = ar1(rowId):colId}
#' \item{random = repId:colId, spatial = rowId:ar1(colId)}
#' \item{random = repId:rowId + repId:colId, spatial = ar1(rowId):ar1(colId)}
#' }
#' If there are no replicates in the model, repId is left out from the random
#' parts above.
#'
#' @param TP An object of class TP.
#' @param trait A character string indicating the trait used as response
#' variable in the model.
#' @param timePoints A character or numeric vector indicating the timePoints
#' to be modeled. When using a character string to reference a timePoint, the
#' value has to be an exact match to one of the existing timePoints. When using
#' a number it will be matched by its number in the timePoints attribute of the
#' TP object.
#' @param covariates A character vector indicating the variables used as
#' fixed effects in the model.
#' @param geno.decomp A character vector indicating the variables used to
#' group the genotypes in the model.
#' @param what A character vector specifying whether "genotype" should
#' be fitted as "random" or "fixed" effect. Note that when using SpATS and
#' geno.decomp fitting a model with genotype as "fixed" effect is not possible.
#' @param useCheck Should check genotypes be used as an extra factor in the
#' model?
#' @param useRepId Should repId be used as a fixed effect in the model? When
#' fitting a spatial model repId is also added as an interaction term with
#' rowId and colId in the random part of the model.
#' @param engine A character string indicating the engine used to fit the
#' models.
#' @param spatial Should a spatial model be fitted for asreml?
#'
#' @return An object of class fitMod, a list of fitted models.
#'
#' @references
#' Maria Xose Rodriguez-Alvarez, Martin P. Boer, Fred A. van Eeuwijk, Paul H.C.
#' Eilers (2017). Correcting for spatial heterogeneity in plant breeding
#' experiments with P-splines. Spatial Statistics
#' \url{https://doi.org/10.1016/j.spasta.2017.10.003}
#' @references
#' Butler, D. G., et al. (2018).ASReml-R Reference Manual Version 4. VSN
#' International Ltd, http://asreml.org
#'
#' @export
fitModels <- function(TP,
                      trait,
                      timePoints = names(TP),
                      covariates = NULL,
                      geno.decomp = NULL,
                      what = c("random", "fixed"),
                      useCheck = FALSE,
                      useRepId = FALSE,
                      engine = c("SpATS", "asreml"),
                      spatial = FALSE) {
  ## Checks.
  if (missing(TP) || !inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  ## Check time points and convert numerical input to corresponding
  ## character values.
  timePoints <- chkTimePoints(TP, timePoints)
  ## Restrict TP to selected time points.
  TP <- TP[timePoints]
  ## Check trait.
  if (missing(trait) || !is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string.\n")
  }
  if (!all(sapply(X = TP, FUN = hasName, name = trait))) {
    stop(trait, " should be a column in TP for all timePoints.\n")
  }
  ## Check covariates.
  if (useRepId) {
    ## Add repId to covariates so it is used as fixed effect.
    covariates <- c(covariates, "repId")
  }
  if (!is.null(covariates)) {
    if (!is.character(covariates)) {
      stop("covariates should be a character vector.\n")
    }
    for (covar in covariates) {
      if (!all(sapply(X = TP, FUN = hasName, name = covar))) {
        stop(covar, " should be columns in TP for all timePoints.\n")
      }
    }
  }
  ## Check geno.decomp.
  if (!is.null(geno.decomp)) {
    if (!is.character(geno.decomp)) {
      stop("geno.decomp should be a character vector.\n")
    }
    for (gd in geno.decomp) {
      if (!all(sapply(X = TP, FUN = hasName, name = gd))) {
        stop(gd, " should be a column for all timePoints.\n")
      }
    }
  }
  what <- match.arg(what)
  genoRand <- what == "random"
  if (useCheck) {
    if (!all(sapply(X = TP, FUN = hasName, name = "check"))) {
      stop("check should be a column in TP for all timePoints.\n")
    }
  }
  engine <- match.arg(engine)
  ## For spatial models spatial columns are required.
  if (engine == "SpATS" || (engine == "asreml" && spatial)) {
    spatCols <- c("rowId", "colId", "rowNum", "colNum")
    if (!all(sapply(X = TP, FUN = hasName, name = trait))) {
      stop(spatCols, " should be a columns in TP for all timePoints when ",
           "fitting spatial models.\n")
    }
  }
  ## Extract timepoints attribute for re-adding in the end.
  timePoints <- attr(TP, which = "timePoints")
  ## If geno.decomp is used genotype and covariates have to be replaced by
  ## an interaction of genotype and covariates with the geno.decomp variables.
  ## Construct an interaction of all variables in geno.decomp.
  if (length(geno.decomp) > 0) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      timePoint[["geno.decomp"]] <- interaction(timePoint[geno.decomp],
                                                sep = "_")
      return(timePoint)
    })
    ## Set geno.decomp to newly constructed variable.
    geno.decomp <- "geno.decomp"
  }
  ## Get column containing genotype.
  genoCol <- if (useCheck) "genoCheck" else "genotype"
  ## Replace genotype and covariates by their interaction with geno.decomp.
  if (!is.null(geno.decomp)) {
    TP <- lapply(X = TP, FUN = function(timePoint) {
      timePoint[[genoCol]] <- interaction(timePoint[[geno.decomp]],
                                          timePoint[[genoCol]], sep = "_")
      for (covar in covariates) {
        timePoint[[covar]] <- interaction(timePoint[[geno.decomp]],
                                          timePoint[[covar]], sep = "_")
      }
      if (useCheck) {
        timePoint[["check"]] <- interaction(timePoint[[geno.decomp]],
                                            timePoint[["check"]], sep = "_")
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
  fixedForm <- formula("~ 1")
  if (!is.null(covariates)) {
    fixedForm <- update(fixedForm,
                        paste("~ . +" , paste(c(covariates), collapse = "+")))
  }
  if (useCheck) {
    fixedForm <- update(fixedForm, "~ . + check")
  }
  if (!is.null(geno.decomp) && is.null(covariates) && !useCheck) {
    fixedForm <- update(fixedForm, "~ . + geno.decomp")
  }
  if (engine == "SpATS") {
    if (useRepId) {
      randForm <- formula("~ repId:rowId + repId:colId")
    } else {
      randForm <- formula("~ rowId + colId")
    }
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
      nseg = c(nlevels(modDat[["colId"]]), nlevels(modDat[["rowId"]])) / 2
      ## Fit and return the model.
      SpATS::SpATS(response = trait, fixed = fixedForm,
                   random = randForm,
                   spatial = ~ SpATS::PSANOVA(colNum, rowNum, nseg = nseg,
                                              nest.div = c(2, 2)),
                   genotype = genoCol, genotype.as.random = genoRand,
                   geno.decomp = geno.decomp, data = modDat,
                   control = list(maxit = 50, tolerance = 1e-03,
                                  monitoring = 0))
    })
  } else if (engine == "asreml") {
    ## fixed in asreml needs response variable on lhs of formula.
    fixedForm <- update(fixedForm, paste(trait, "~ ."))
    ## Genotype should be specified in either the fixed or the random part
    ## of the model.
    if (genoRand) {
      ## Construct formula for random part of the model.
      randForm <- formula(paste("~ ", if (is.null(geno.decomp)) genoCol else
        paste0("at(", geno.decomp, "):", genoCol)))
    } else {
      ## For genotype fixed the base random formula is empty.
      ## Genotype is added to the fixedForm.
      fixedForm <- update(fixedForm,
                          paste("~ . + ", if (is.null(geno.decomp)) genoCol else
                            paste0("at(", geno.decomp, "):", genoCol)))
      randForm <- ~NULL
    }
    if (!spatial) {
      ## Loop on timepoint to run asreml.
      fitMods <- lapply(X = TP, function(timePoint) {
        message(timePoint[["timePoint"]][1])
        ## Only keep columns needed for analysis.
        modDat <- droplevels(timePoint)
        ## Run model.
        asrFit <- asreml::asreml(fixed = fixedForm, random = randForm,
                                 data = modDat, trace = FALSE, maxiter = 200,
                                 na.action = asreml::na.method(x = "include"))
        ## evaluate call terms so predict can be run.
        ## The first (unnamed) item in call contains the full asreml function.
        ## This is replaced by a named reference to drastically reduce output
        ## size.
        asrFit$call[[1]] <- quote(asreml::asreml)
        asrFit$call$fixed <- eval(asrFit$call$fixed)
        asrFit$call$random <- eval(asrFit$call$random)
        asrFit$call$data <- substitute(modDat)
        return(asrFit)
      })
    } else {
      ## First find the best spatial model over all timePoints by fitting
      ## all models for 20% of the time points (with a minimum of 10)
      bestNum <- min(max(length(TP) / 5, 10), length(TP))
      bestMod <- character()
      for (i in round(seq(1, length(TP), length.out = bestNum))) {
        timePoint <- TP[[i]]
        ## Only keep columns needed for analysis.
        modDat <- droplevels(timePoint)
        asrFitSpat <- bestSpatMod(modDat = modDat, traits = trait,
                                  fixedForm = fixedForm, randomForm = randForm)
        bestMod <- c(bestMod, attr(asrFitSpat[["sumTab"]][[trait]], "chosen"))
      }
      bestMod <- names(which.max(table(bestMod)))
      fitMods <- setNames(vector(mode = "list", length = length(TP)), names(TP))
      for (i in seq_along(TP)) {
        timePoint <- TP[[i]]
        message(timePoint[["timePoint"]][1])
        ## Only keep columns needed for analysis.
        modDat <- droplevels(timePoint)
        asrFitSpat <- bestSpatMod(modDat = modDat, traits = trait,
                                  fixedForm = fixedForm, randomForm = randForm,
                                  useRepId = useRepId, spatTerms = bestMod)
        asrFit <- asrFitSpat[["fitMods"]][[trait]]
        attr(x = asrFit, which = "sumTab") <- asrFitSpat[["sumTab"]]
        ## evaluate call terms so predict can be run.
        ## The first (unnamed) item in call contains the full asreml function.
        ## This is replaced by a named reference to drastically reduce output
        ## size.
        asrFit$call[[1]] <- quote(asreml::asreml)
        asrFit$call$fixed <- eval(asrFit$call$fixed)
        asrFit$call$random <- eval(asrFit$call$random)
        asrFit$call$data <- eval(asrFit$call$data)
        fitMods[[i]] <- asrFit
      }
    }
  }
  return(createFitMod(fitMods,
                      what = what,
                      useRepId = useRepId,
                      timePoints = timePoints))
}

#' Helper function for calculating best spatial model using asreml.
#' @noRd
#' @keywords internal
bestSpatMod <- function(modDat,
                        traits,
                        criterion = "AIC",
                        fixedForm,
                        randomForm,
                        useRepId = FALSE,
                        spatTerms = c("none", "AR1(x)id", "id(x)AR1",
                                      "AR1(x)AR1"),
                        ...) {
  dotArgs <- list(...)
  spatTerms <- match.arg(spatTerms, several.ok = TRUE)
  ## Increase max number of iterations for asreml.
  maxIter <- 200
  ## Add empty observations.
  TPTab <- as.data.frame(table(modDat[["colId"]], modDat[["rowId"]]))
  TPTab <- TPTab[TPTab$Freq == 0, , drop = FALSE]
  if (nrow(TPTab) > 0) {
    extObs <- setNames(as.data.frame(matrix(nrow = nrow(TPTab),
                                            ncol = ncol(modDat))),
                       colnames(modDat))
    extObs[["timePoint"]] <- modDat[["timePoint"]][1]
    extObs[, c("colId", "rowId")] <- TPTab[, c("Var1", "Var2")]
    extObs[, c("colNum", "rowNum")] <-
      c(as.numeric(levels(TPTab[, "Var1"]))[TPTab[, "Var1"]],
        as.numeric(levels(TPTab[, "Var2"]))[TPTab[, "Var2"]])
    modDat <- rbind(modDat, extObs)
  }
  ## modDat needs to be sorted by row and column to prevent asreml from crashing.
  modDat <- modDat[order(modDat[["rowId"]], modDat[["colId"]]), ]
  spatCh <- c("none", "AR1(x)id", "id(x)AR1", "AR1(x)AR1")
  spatSel <- which(spatCh %in% spatTerms)
  spatCh <- spatCh[spatSel]
  spatTerm <- c(NA, paste("~", c("ar1(rowId):colId", "rowId:ar1(colId)",
                                 "ar1(rowId):ar1(colId)")))[spatSel]
  ## Create empty base lists.
  fitMods <- spatial <- sumTab <- setNames(vector(mode = "list",
                                                  length = length(traits)),
                                           traits)
  btCols <- c("spatial", "AIC", "BIC", "row", "col", "error", "converge")
  for (trait in traits) {
    ## Reset criterion to Inf.
    criterionBest <- Inf
    ## Create data.frame for storing summary for current trait.
    modSum <- as.data.frame(matrix(nrow = length(spatCh), ncol = length(btCols),
                                   dimnames = list(NULL, btCols)))
    ## Fit model with genotype random for all different random/spatial terms.
    for (i in seq_along(spatTerm)) {
      ## Add extra random term to random part.
      if (useRepId) {
        randForm <- update(randomForm, "~ . + repId:rowId + repId:colId")
      } else {
        randForm <- update(randomForm, "~ . + rowId + colId")
      }
      asrArgs <- c(list(fixed = fixedForm, random = randForm, aom = TRUE,
                        data = modDat, maxiter = maxIter, trace = FALSE,
                        na.action = asreml::na.method(x = "include")),
                   dotArgs)
      if (!is.na(spatTerm[i])) {
        asrArgs[["residual"]] <- formula(spatTerm[i])
      }
      capture.output(fitMod <- tryCatchExt(do.call(what = asreml::asreml,
                                                   args = asrArgs)),
                     file = tempfile())
      if (!is.null(fitMod$warning)) {
        fitMod <- chkLastIter(fitMod)
        fitMod <- wrnToErr(fitMod)
      }
      if (length(fitMod$warning) != 0) {
        warning(paste0("Warning in asreml for model ", spatCh[i],
                       " genotype random, trait ", trait, " in timePoint ",
                       modDat[["timePoint"]][1], ":\n", fitMod$warning,
                       "\n"), call. = FALSE)
      }
      if (is.null(fitMod$error)) {
        fitMod <- fitMod$value
      } else {
        warning(paste0("Error in asreml for model ", spatCh[i],
                       " genotype random, trait ", trait, " in timePoint ",
                       modDat[["timePoint"]][1], ":\n", fitMod$error,
                       "\n"), call. = FALSE)
        fitMod <- NULL
      }
      ## Fill model summary table.
      modSum[i, "spatial"] <- spatCh[i]
      modSum[i, "converge"] <- isTRUE(!is.null(fitMod) & fitMod$converge)
      if (!is.null(fitMod)) {
        summ <- summary(fitMod)$varcomp["component"]
        modSum[i, "AIC"] <- -2 * fitMod$loglik + 2 * nrow(summ)
        modSum[i, "BIC"] <- -2 * fitMod$loglik +
          log(length(fitted(fitMod))) * nrow(summ)
        ## Row and column output differs for regular/non-regular.
        ## Always max. one of the possibilities is in summary so rowVal and
        ## colVal are always a single value.
        rowVal <- summ[rownames(summ) == "rowId:colId!rowId!cor", ]
        modSum[i, "row"] <- ifelse(length(rowVal) == 0, NA, rowVal)
        colVal <- summ[rownames(summ) == "rowId:colId!colId!cor", ]
        modSum[i, "col"] <- ifelse(length(colVal) == 0, NA, colVal)
        modSum[i, "error"] <- summ[rownames(summ) %in%
                                     c("units!R", "rowId:colId!R"), ]
        ## If current model is better than best so far based on chosen criterion
        ## define best model as current model.
        if (fitMod$converge) {
          if (criterion == "AIC") {
            criterionCur <- modSum[i, "AIC"]
          } else {
            criterionCur <- modSum[i, "BIC"]
          }
        }
        if (criterionCur < criterionBest) {
          bestModTr <- fitMod
          ## Evaluate call terms in bestModTr so predict can be run.
          ## The first (unnamed) item in call contains the full asreml function.
          ## This is replaced by a named reference to drastically reduce output
          ## size.
          ## Needs to be called in every iteration to prevent final result
          ## from always having the values of the last iteration.
          bestModTr$call[[1]] <- quote(asreml::asreml)
          bestModTr$call$fixed <- eval(bestModTr$call$fixed)
          bestModTr$call$random <- eval(bestModTr$call$random)
          bestModTr$call$residual <- eval(bestModTr$call$residual)
          bestModTr$call$data <- substitute(modDat)
          criterionBest <- criterionCur
          bestMod <- spatCh[i]
        }
      }
    }
    fitMods[[trait]] <- bestModTr
    spatial[[trait]] <- spatCh[bestMod]
    attr(x = modSum, which = "chosen") <- bestMod
    sumTab[[trait]] <- modSum
  } # End for traits.
  return(list(fitMods = fitMods, spatial = spatial, sumTab = sumTab))
}

