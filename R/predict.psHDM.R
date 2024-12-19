#' Predict the P-Splines Hierarchical Curve Data Model
#'
#' Function that predicts the P-spline Hierarchical Curve Data Model (see
#' \code{\link{fitSplineHDM}}) on a dense grid. It provides standard errors
#' for curves at each level of the hierarchy. User has to be aware that
#' standard errors at the plot level demand large memory. We suggest set
#' that option at the \code{FALSE} level
#'
#' @param object An object of class "psHDM" as obtained after fitting
#' (\code{\link{fitSplineHDM}}) the P-spline Hierarchical Curve Data Model
#' @param newtimes A numeric vector with timepoints at which predictions are
#' desired
#' @param pred A list that controls the hierarchical levels at which
#' predictions are desired (population/genotypes/plots).  The default is
#' \code{TRUE}.
#' @param se A list that controls the hierarchical levels at which standard
#' errors are desired (population/genotypes/plots).  The default is
#' \code{TRUE} except at the plot level.
#' @param ... Not used.
#' @param trace An optional value that controls the function trace.
#' The default is \code{TRUE}.
#'
#' @returns An object of class \code{psHDM}, a list with the following outputs:
#' predict.psHDM
#' \code{newtimes} A numeric vector with the timepoints at which predictions
#' and/or standard errors have been obtained.
#' \code{popLevel} A data.frame with the estimated population trajectories
#' and first and second order derivatives, and if required their respective
#' standard errors, at the \code{newtimes}.
#' \code{genoLevel} A data.frame with the estimated genotype-specific
#' deviations and trajectories and their respective first and second order
#' derivatives, and if required their respective standard errors,
#' at the \code{newtimes}.
#' \code{plotLevel} A data.frame with the estimated plot-specific
#' deviations and trajectories and their respective first and second order
#' derivatives, and if required their respective standard errors,
#' at the \code{newtimes}.
#' \code{plotObs} A data.frame with the raw data at the original timepoints.
#'
#' @examples
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## We need to specify the genotype-by-treatment interaction.
#' ## Treatment: water regime (WW, WD).
#' spatCorrectedArch[["treat"]] <- substr(spatCorrectedArch[["geno.decomp"]],
#'                                       start = 1, stop = 2)
#' spatCorrectedArch[["genoTreat"]] <-
#'   interaction(spatCorrectedArch[["genotype"]],
#'              spatCorrectedArch[["treat"]], sep = "_")
#'
#' ## Fit P-Splines Hierarchical Curve Data Model for selection of genotypes.
#' fit.psHDM  <- fitSplineHDM(inDat = spatCorrectedArch,
#'                           trait = "LeafArea_corr",
#'                           genotypes = c("GenoA14_WD", "GenoA51_WD",
#'                                        "GenoB11_WW", "GenoB02_WD",
#'                                        "GenoB02_WW"),
#'                           time = "timeNumber",
#'                           pop = "geno.decomp",
#'                           genotype = "genoTreat",
#'                           plotId = "plotId",
#'                           difVar = list(geno = FALSE, plot = FALSE),
#'                           smoothPop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothGeno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothPlot = list(nseg = 4, bdeg = 3, pord = 2),
#'                           weights = "wt",
#'                           trace = FALSE)
#'
#' ## Predict the P-Splines Hierarchical Curve Data Model on a dense grid
#' ## with standard errors at the population and genotype levels
#' pred.psHDM <- predict(object = fit.psHDM,
#'                      newtimes = seq(min(fit.psHDM$time[["timeNumber"]]),
#'                                    max(fit.psHDM$time[["timeNumber"]]),
#'                                    length.out = 100),
#'                      pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
#'                      se = list(pop = TRUE, geno = TRUE, plot = FALSE))
#'
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#'
#' ## Plots at population level.
#' plot(pred.psHDM,
#'     plotType = "popTra")
#'
#' ## Plots at genotype level.
#' plot(pred.psHDM,
#'     plotType = "popGenoTra")
#'
#' ## Plots of derivatives at genotype level.
#' plot(pred.psHDM,
#'     plotType = "popGenoDeriv")
#'
#' ## Plots of deviations at genotype level.
#' plot(pred.psHDM,
#'     plotType = "genoDev")
#'
#' ## Plots at plot level.
#' plot(pred.psHDM,
#'     plotType = "genoPlotTra")
#'
#' @references Pérez-Valencia, D.M., Rodríguez-Álvarez, M.X., Boer, M.P. et al.
#' A two-stage approach for the spatio-temporal analysis of high-throughput
#' phenotyping data. Sci Rep 12, 3177 (2022). \doi{10.1038/s41598-022-06935-9}
#'
#' @family functions for fitting hierarchical curve data models
#'
#' @export
predict.psHDM <- function(object,
                          newtimes,
                          pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                          se = list(pop = TRUE, geno = TRUE, plot = FALSE),
                          trace = TRUE,
                          ...) {
  ## Checks.
  if (missing(newtimes)) {
    xp <- object$time[["timeNumber"]]
  } else {
    if (!is.vector(newtimes) || !is.numeric(newtimes)) {
      stop("newtimes should be a numerical vector.\n")
    }
    xp <- newtimes
  }
  if (!is.list(pred) || length(pred) != 3 ||
      !(setequal(names(pred), c("pop", "geno", "plot")))) {
        stop("pred should be a named list of length 3.\n")
  }
  if (!is.list(se) || length(se) != 3 ||
      !(setequal(names(se), c("pop", "geno", "plot")))) {
    stop("se should be a named list of length 3.\n")
  }
  if (isTRUE(pred$plot) && !(isTRUE(pred$geno) && isTRUE(pred$pop))) {
    stop("Predictions at plot level can only be made if predictions are ",
         "also made at geno and pop level.\n")
  }
  if (isTRUE(pred$geno) && !isTRUE(pred$pop)) {
    stop("Predictions at geno level can only be made if predictions are ",
         "also made at pop level.\n")
  }
  if (isTRUE(se$plot) && !(isTRUE(se$geno) && isTRUE(se$pop))) {
    stop("Standard errors at plot level can only be computed if standard ",
         "errors are also computed at geno and pop level.\n")
  }
  if (isTRUE(se$geno) && !isTRUE(se$pop)) {
    stop("Standard errors at geno level can only be computed if standard ",
         "errors are also computed at pop level.\n")
  }
  if (isTRUE(se$pop) && !isTRUE(pred$pop)) {
    stop("Standard errors at population level can only be computed ",
         "if predictions are also made at population level.\n")
  }
  if (isTRUE(se$geno) && !isTRUE(pred$geno)) {
    stop("Standard errors at genotype level can only be computed ",
         "if predictions are also made at genotype level.\n")
  }
  if (isTRUE(se$plot) && !isTRUE(pred$plot)) {
    stop("Standard errors at plot level can only be computed ",
         "if predictions are also made at plot level.\n")
  }
  ## Output data.
  res <- list(newtimes = xp)
  ## Normalize time
  xp <- xp - min(xp) + 1
  ## Number of parameters: fixed and random (for each component)
  np <- object$dim
  npComp <- c(np[1] + np[2], np[3] + np[4], np[5] + np[6])
  names(npComp) <- c("pop", "geno", "plot")
  npE <- cumsum(npComp)
  npS <- npE - npComp + 1
  ## Predictions.
  ## Functions at population level.
  if (isTRUE(pred$pop)) {
    popLevel <- mixmodToBsplinePred(what = "pop", object = object,
                                    npS, npE, xp, dev = FALSE)
    if (trace) {
      print("Population-specific growth curves OK")
    }
    if (isTRUE(se$pop)) {
      popLevel$pred <-
        append(popLevel$pred,
               standardErrors(Tm = popLevel$Tm, B = popLevel$B,
                              Bd1 = popLevel$Bd1, Bd2 = popLevel$Bd2,
                              what = "pop", dev = FALSE, object = object,
                              npS, npE))
      if (trace) {
        print("Standard errors for population-specific growth curves OK")
      }
    }
    ## Data frame with all the information at population level
    res <- append(res, listToDf(object1 = popLevel$pred, object2 = object,
                                what = "pop", xp = res$newtimes))
  }
  ## Functions at genotype level.
  if (isTRUE(pred$geno)) {
    ## Genotype-specific deviations and first- and second-order derivatives
    genoDev <- mixmodToBsplinePred(what = "geno", object = object,
                                   npS, npE, xp, dev = TRUE)
    genoLevel <- list(genoDev = genoDev$pred)
    if (trace) {
      print("Genotype-specific deviations OK")
    }
    ## Genotype-specific growth curves and first- and second-order derivatives
    ## Contrast matrix: Assign genotypes to populations
    if(length(object$popLevs) == 1) {
      mmGenoPop <- matrix(1, ncol = 1, nrow = length(object$genoLevs))
    } else {
      genoPop <- rep(as.factor(1:length(object$popLevs)), object$nGenoPop)
      mmGenoPop <- Matrix::sparse.model.matrix(~ 0 + genoPop) # The contrast matrix changes here!!!!!!
    }
    TGeno <- Matrix::bdiag(popLevel$Tm, genoDev$Tm)
    genoTra <- mixmodToBsplinePred(what = "geno", object = object, npS,
                                   npE, xp, Tmat = TGeno,
                                   modMat = mmGenoPop, dev = FALSE,
                                   Bbasis = list(B = genoDev$B,
                                                 Bd1 = genoDev$Bd1,
                                                 Bd2 = genoDev$Bd2))
    genoLevel$genoTra <- genoTra$pred
    if (trace) {
      print ("Genotype-specific growth curves OK")
    }
    if (isTRUE(se$geno)) {
      ## Genotype-specific deviations
      genoLevel$genoDev <-
        append(genoLevel$genoDev,
               standardErrors(Tm = genoDev$Tm, B = genoDev$B,
                              Bd1 = genoDev$Bd1, Bd2 = genoDev$Bd2,
                              what = "geno", dev = TRUE, object = object,
                              npS, npE))
      if (trace) {
        print("Standard errors for genotype-specific deviations OK")
      }
      ## Genotype-specific growth curves
      genoLevel$genoTra <-
        append(genoLevel$genoTra,
               standardErrors(Tm = TGeno, B = genoTra$B,
                              Bd1 = genoTra$Bd1, Bd2 = genoTra$Bd2,
                              what = "geno", dev = FALSE, object = object,
                              npS, npE))
      if (trace) {
        print("Standard errors for genotype-specific growth curves OK")
      }
    }
    ## Data frame with all the information at genotype level.
    res <- append(res, listToDf(object1 = genoLevel, object2 = object,
                                what = "geno", xp = res$newtimes))
  }
  ## Functions at plot level.
  if (isTRUE(pred$plot)) {
    # Plot-specific deviations and first- and second-order derivatives
    plotDev <- mixmodToBsplinePred(what = "plot", object = object, npS,
                                   npE, xp, dev = TRUE)
    plotLevel <- list(plotDev = plotDev$pred)
    if (trace) {
      print("Plot-specific deviations OK")
    }
    ## Plot-specific growth curves.
    ## Contrast matrix: Assign plots to populations.
    if (length(object$popLevs) == 1) {
      mmPlotPop <- matrix(1, ncol = 1, nrow = object$nPlotPop)
    } else {
      plotPop <- rep(as.factor(1:length(object$popLevs)), object$nPlotPop)
      mmPlotPop <- Matrix::sparse.model.matrix(~ 0 + plotPop) # The contrast matrix changes here!!!!!!
    }
    ## Contrast matrix: Assign plots to genotypes.
    if (length(object$genoLevs) == 1) {
      mmPlotGeno <- matrix(1, ncol = 1, nrow = length(object$genoLevs))
    } else {
      plotGeno <- rep(as.factor(1:length(object$genoLevs)), object$nPlotGeno)
      mmPlotGeno <- Matrix::sparse.model.matrix(~ 0 + plotGeno) # The contrast matrix changes here!!!!!!
    }
    TPlot <- Matrix::bdiag(popLevel$Tm, genoDev$Tm, plotDev$Tm)
    plotTra <- mixmodToBsplinePred(what = "plot", object = object, npS,
                                   npE, xp, Tmat = TPlot,
                                   modMat = list(mmPlotPop, mmPlotGeno),
                                   dev = FALSE,
                                   Bbasis = list(B = plotDev$B,
                                                 Bd1 = plotDev$Bd1,
                                                 Bd2 = plotDev$Bd2))
    plotLevel$plotTra <- plotTra$pred
    if (trace) {
      print("Plot-specific growth curves OK")
    }
    if (isTRUE(se$plot)) {
      ## Plot-specific deviations.
      plotLevel$plotDev <-
        append(plotLevel$plotDev,
               standardErrors(Tm = plotDev$Tm, B = plotDev$B,
                              Bd1 = plotDev$Bd1, Bd2 = plotDev$Bd2,
                              what = "plot", dev = TRUE, object = object,
                              npS, npE))
      if (trace) {
        print("Standard errors for plot-specific deviations OK")
      }
      ## Standard errors for plot deviations + geno deviations + population effects.
      plotLevel$plotTra <-
        append(plotLevel$plotTra,
               standardErrors(Tm = TPlot, B = plotTra$B,
                              Bd1 = plotTra$Bd1, Bd2 = plotTra$Bd2,
                              what = "plot", dev = FALSE, object = object,
                              npS, npE))
      if (trace) {
        print("Standard errors for plot-specific growth curves OK")
      }
    }
    ## Data frame with all the information at genotype level.
    res <- append(res,
                  listToDf(object1 = plotLevel, object2 = object,
                           what = "plot", xp = res$newtimes))
  }
  class(res) <- c("psHDM", "list")
  attr(res, which = "trait") <- attributes(object)[["trait"]] #object$trait
  return(res)
}


### Help functions

#' mixmodToBsplinePred
#'
#' @noRd
#' @keywords internal
mixmodToBsplinePred <- function(what = c("pop", "geno", "plot"),
                                object,
                                npS,
                                npE,
                                xp,
                                Tmat = NULL,
                                modMat = NULL,
                                dev = TRUE,
                                Bbasis = NULL){
  ## Predictions.
  ## Transformation matrix Tm, Mixed model coefficients MMCoeff,
  ## theta and B, and fitted/predicted values f.
  whatS <- whatE <- what
  if (isFALSE(dev) && what != "pop") {
    if (what == "geno") {
      whatS <- "pop"
      whatE <- "geno"
    } else {
      whatS <- "pop"
      whatE <- "plot"
    }
  }
  lW <- length(object[[paste0(what, "Levs")]])
  MMW <- paste0("MM", tools::toTitleCase(what))
  if (is.null(Tmat)) {
    Tm <- cbind(Matrix::kronecker(Matrix::Diagonal(lW), object$MM[[MMW]]$U.X),
                Matrix::kronecker(Matrix::Diagonal(lW), object$MM[[MMW]]$U.Z))
  } else{
    Tm <- Tmat
  }
  if (is.null(modMat)) {
    mmat <- Matrix::Diagonal(lW)
  } else if (is.list(modMat)) {
    mmat <- list(mmat1 = modMat[[1]], mmat2 = modMat[[2]])
  } else{
    mmat <- modMat
  }
  MMCoeff <- Matrix::Matrix(object$coeff[npS[whatS]:npE[whatE]],
                            ncol = 1)
  theta <- Tm %*% MMCoeff
  BFull <- function(mmat, what, deriv) {
    MMW <- paste0("MM", tools::toTitleCase(what))
    smW <- paste0("smooth", tools::toTitleCase(what))
    Matrix::kronecker(mmat,
                      spline.bbase(knots = object$MM[[MMW]]$knots,
                                   X. = xp, BDEG. = object$smooth[[smW]]$bdeg,
                                   deriv = deriv))
  }
  if (is.null(Bbasis)) {
    B <- BFull(mmat, what, deriv = 0)
    Bd1 <- BFull(mmat, what, deriv = 1)
    Bd2 <- BFull(mmat, what, deriv = 2)
  } else {
    if (what == "geno"){
      B <- cbind(BFull(mmat, what = "pop", deriv = 0), Bbasis$B)
      Bd1 <- cbind(BFull(mmat, what = "pop", deriv = 1), Bbasis$Bd1)
      Bd2 <- cbind(BFull(mmat, what = "pop", deriv = 2), Bbasis$Bd2)
    } else if(what == "plot") {
      B <- cbind(BFull(mmat$mmat1, what = "pop", deriv = 0),
                 BFull(mmat$mmat2, what = "geno", deriv = 0),
                 Bbasis$B)
      Bd1 <- cbind(BFull(mmat$mmat1, what = "pop", deriv = 1),
                   BFull(mmat$mmat2, what = "geno", deriv = 1),
                   Bbasis$Bd1)
      Bd2 <- cbind(BFull(mmat$mmat1, what = "pop", deriv = 2),
                   BFull(mmat$mmat2, what = "geno", deriv = 2),
                   Bbasis$Bd2)
    }
  }
  f <- matrix(B %*% theta, ncol = lW) # Note that f == eta_what
  fd1 <- matrix(Bd1 %*% theta, ncol = lW)
  fd2 <- matrix(Bd2 %*% theta, ncol = lW)
  aux <- lapply(list(f = f, fd1 = fd1, fd2 = fd2), function(x) {
    colnames(x) <- object[[paste0(what, "Levs")]]
    return(x)
  })
  ## Object to be returned.
  res <- list(level = what,
              Tm = Tm,
              B = B,
              Bd1 = Bd1,
              Bd2 = Bd2,
              pred = aux)
  return(res)
}

#' standardErrors
#'
#' @noRd
#' @keywords internal
standardErrors <- function(Tm,
                           B,
                           Bd1,
                           Bd2,
                           what = c("pop", "geno", "plot"),
                           dev = TRUE,
                           object,
                           npS,
                           npE) {
  whatS <- whatE <- what
  if (isFALSE(dev) && what != "pop"){
    if (what == "geno") {
      whatS <- "pop"
      whatE <- "geno"
    } else {
      whatS <- "pop"
      whatE <- "plot"
    }
  }
  Vp <- spam::chol2inv(object$cholHn)
  lW <- length(object[[paste0(what, "Levs")]])
  seTheta <- Matrix::Matrix(Tm) %*%
    Matrix::Matrix(Vp[npS[whatS]:npE[whatE],
                      npS[whatS]:npE[whatE]]) %*% Matrix::t(Tm)
  sef <- matrix(sqrt(Matrix::colSums(
    Matrix::t(B) * Matrix::tcrossprod(seTheta, B))), ncol = lW)
  sefd1 <- matrix(sqrt(Matrix::colSums(
    Matrix::t(Bd1) * Matrix::tcrossprod(seTheta, Bd1))), ncol = lW)
  sefd2 <- matrix(sqrt(Matrix::colSums(
    Matrix::t(Bd2) * Matrix::tcrossprod(seTheta, Bd2))), ncol = lW)
  res <- list(sef = sef,
              sefd1 = sefd1,
              sefd2 = sefd2)
  res <- lapply(res, function(x){
    colnames(x) <- object[[paste0(what, "Levs")]]
    return(x)
  })
  return(res)
}

#' listToDf
#'
#' @noRd
#' @keywords internal
listToDf <- function(object1,
                     object2,
                     what,
                     xp) {
  if(!inherits(object2, "psHDM")) {
    stop("The object class is not correct")
  }
  res <- list()
  if (hasName(x = object2$time, name = "timePoint")) {
    ## Create time point range.
    ## Has to take into account irregular grids for xp.
    timePointRange <- min(object2$time[["timePoint"]]) +
      (xp - min(object2$time[["timeNumber"]])) /
      (diff(range(object2$time[["timeNumber"]]))) *
      diff(range(object2$time[["timePoint"]]))
  }
  ## Population-specific growth curves.
  if (what == "pop") {
    dfPopTra <- data.frame(timeNumber = rep(xp, length(object2$popLevs)),
                           pop = rep(object2$popLevs, each = length(xp)),
                           fPop = c(object1$f),
                           fPopDeriv1 = c(object1$fd1),
                           fPopDeriv2 = c(object1$fd2))
    if (!is.null(object1$sef)) {
      dfPopTra$sePop <- c(object1$sef)
      dfPopTra$sePopDeriv1 <- c(object1$sefd1)
      dfPopTra$sePopDeriv2 <- c(object1$sefd2)
    }
    if (hasName(x = object2$time, name = "timePoint")) {
      dfPopTra[["timePoint"]] <- rep(timePointRange, length(object2$popLevs))
      dfPopTra <- dfPopTra[c("timeNumber", "timePoint",
                             setdiff(colnames(dfPopTra),
                                     c("timeNumber", "timePoint")))]
    }
    res$popLevel <- dfPopTra
  }
  ## Genotypic-specific growth curves and deviations
  if (what == "geno") {
    dfGenoTra <- data.frame(timeNumber = rep(xp, length(object2$genoLevs)),
                            pop = rep(object2$popLevs,
                                      object2$nGenoPop * length(xp)),
                            genotype = rep(object2$genoLevs, each = length(xp)),
                            fGeno = as.vector(object1$genoTra$f),
                            fGenoDeriv1 = as.vector(object1$genoTra$fd1),
                            fGenoDeriv2 = as.vector(object1$genoTra$fd2),
                            fGenoDev = as.vector(object1$genoDev$f),
                            fGenoDevDeriv1 = as.vector(object1$genoDev$fd1),
                            fGenoDevDeriv2 = as.vector(object1$genoDev$fd2))
    if (!is.null(object1$genoDev$sef)) {
      dfGenoTra$seGeno <- as.vector(object1$genoTra$sef)
      dfGenoTra$seGenoDeriv1 <- as.vector(object1$genoTra$sefd1)
      dfGenoTra$seGenoDeriv2 <- as.vector(object1$genoTra$sefd2)
      dfGenoTra$seGenoDev <- as.vector(object1$genoDev$sef)
      dfGenoTra$seGenoDevDeriv1 <- as.vector(object1$genoDev$sefd1)
      dfGenoTra$seGenoDevDeriv2 <- as.vector(object1$genoDev$sefd2)
    }
    if (hasName(x = object2$time, name = "timePoint")) {
      dfGenoTra[["timePoint"]] <- rep(timePointRange, length(object2$genoLevs))
      dfGenoTra <- dfGenoTra[c("timeNumber", "timePoint",
                               setdiff(colnames(dfGenoTra),
                                       c("timeNumber", "timePoint")))]
    }
    res$genoLevel <- dfGenoTra
  }
  ## Plot-specific growth curves and deviations
  if (what == "plot") {
    dfPlotTra <- data.frame(timeNumber = rep(xp, sum(object2$nPlotGeno)),
                            pop = rep(object2$popLevs,
                                      object2$nPlotPop * length(xp)),
                            genotype = rep(object2$genoLevs,
                                           object2$nPlotGeno * length(xp)),
                            plotId = rep(object2$plotLevs, each = length(xp)),
                            fPlot = as.vector(object1$plotTra$f),
                            fPlotDeriv1 = as.vector(object1$plotTra$fd1),
                            fPlotDeriv2 = as.vector(object1$plotTra$fd2),
                            fPlotDev = as.vector(object1$plotDev$f),
                            fPlotDevDeriv1 = as.vector(object1$plotDev$fd1),
                            fPlotDevDeriv2 = as.vector(object1$plotDev$fd2))
    if (!is.null(object1$plotDev$sef)){
      dfPlotTra$sePlot <- as.vector(object1$plotTra$sef)
      dfPlotTra$sePlotDeriv1 <- as.vector(object1$plotTra$sefd1)
      dfPlotTra$sePlotDeriv2 <- as.vector(object1$plotTra$sefd2)
      dfPlotTra$sePlotDev <- as.vector(object1$plotDev$sef)
      dfPlotTra$sePlotDevDeriv1 <- as.vector(object1$plotDev$sefd1)
      dfPlotTra$sePlotDevDeriv2 <- as.vector(object1$plotDev$sefd2)
    }
    if (hasName(x = object2$time, name = "timePoint")) {
      dfPlotTra[["timePoint"]] <- rep(timePointRange, sum(object2$nPlotGeno))
      dfPlotTra <- dfPlotTra[c("timeNumber", "timePoint",
                               setdiff(colnames(dfPlotTra),
                                       c("timeNumber", "timePoint")))]
    }
    res$plotLevel <- dfPlotTra
    if (isTRUE(setequal(xp, object2$time[["timeNumber"]]))) {
      res$plotLevel$obsPlot <- c(do.call("cbind", object2$y))
    } else {
      ## Raw plot growth curves.
      ## (Observed data: in the raw time not in newtimes (xp)).
      dfPlotObs <-
        data.frame(timeNumber = rep(object2$time[["timeNumber"]],
                                    sum(object2$nPlotGeno)),
                   pop = rep(object2$popLevs,
                             object2$nPlotPop * nrow(object2$time)),
                   genotype = rep(object2$genoLevs,
                                  object2$nPlotGeno * nrow(object2$time)),
                   plotId = rep(object2$plotLevs, each = nrow(object2$time)),
                   obsPlot = c(do.call("cbind", object2$y)))
      if (hasName(x = object2$time, name = "timePoint")) {
        dfPlotObs[["timePoint"]] <- rep(object2$time[["timePoint"]],
                                        sum(object2$nPlotGeno))
        dfPlotObs <- dfPlotObs[c("timeNumber", "timePoint",
                                 setdiff(colnames(dfPlotObs),
                                         c("timeNumber", "timePoint")))]
      }
      res$plotObs <- dfPlotObs
    }
  }
  return(res)
}

