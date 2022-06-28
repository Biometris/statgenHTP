# Function ----------------------------------------------------------------

#' fitSplineHDM
#'
#' Fit the P-spline Hierarchical Curve Data Model used in the second stage of
#' the two-stage approach proposed by Pérez-Valencia et al. (2022). This model
#' assumes a three-level hierarchical structure in the data, with plants nested
#' in genotypes, genotypes nested in populations. The input for this function
#' is the spatially corrected data, as obtained from the first stage of the
#' approach (see \code{\link{fitModels}} and \code{\link{getCorrected}}).
#' The number of segments is chosen by the user, as well as the B-spline degree,
#' and the penalty order for the three-levels of the hierarchy. The user can
#' also decide if different variances for random effects at genotype (separately
#' for each population) and plant (separately for each genotype) levels are
#' desired. The function outputs are fitted values (time series of trajectories
#' and deviations) and their first and second derivatives for the three-levels
#' of the hierarchy. The outputs can then be used to estimate relevant parameters
#' from the curves for further analysis (see \code{\link{estimateSplineParameters}}).
#'
#' @param inDat A data.frame with corrected spatial data.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param time A character string indicating the timepoints at which the trait
#' is measured for each plant.
#' @param pop A character string indicating the the populations to which each
#' genotype/variety belongs. This variable must be a factor in the data frame.
#' @param genotype A character string indicating the populations to which each
#' genotype/variety belongs. This variable must be a factor in the data frame.
#' @param plotId A character string indicating the genotypes/varieties to which
#' each plant/plot/individual belongs. This variable must be a factor in the
#' data frame.
#' @param difVar Should different variances for random effects at genotype
#' (separately for each population) and plant level (separately for each
#' genotype) be considered?.
#' @param smoothPop A list specifying the P-Spline model at the population
#' level (nseg: number of segments; bdeg: degree of the B-spline basis; pord:
#' penalty order).
#' @param smoothGeno A list specifying the P-Spline model at the genotype
#' level.
#' @param smoothPlot A list specifying the P-Spline model at the plant level.
#' @param offset An optional numerical vector containing an a priori known
#' component to be included in the linear predictor during fitting. The default
#' is \code{NULL}.
#' @param weights A character string indicating the weights to be used in the
#' fitting process (for error propagation from first stage to second stage).
#' By default, the weights are considered to be one. The default is \code{NULL}.
#' @param family An object of class \code{family} specifying the distribution
#' and link function. The default is \code{gaussian()}.
#' @param maxit An optional value that controls the maximum number of iterations
#' of the algorithm. The default is 200.
#' @param trace An optional value that controls the function trace.
#' The default is \code{TRUE}.
#' @param thr An optional value that controls the convergence threshold of the
#' algorithm. The default is 1.e-03.
#'
#' @return An object of class \code{psHDM}, a list with the following outputs:
#' \code{time}, a numeric vector with the timepoints.
#' \code{genoLevs}, a vector with the names of the genotypes.
#' \code{popLevs}, a vector with the names of the populations
#' \code{plotLevs}, a vector with the names of the plants
#' \code{nPlotPop}, a numeric vector with the number of plants per
#' population.
#' \code{nGenoPop}, a numeric vector with the number of genotypes per
#' population.
#' \code{nPlotGeno}, a numeric vector with the number of plants per
#' genotype.
#' \code{MM}, a list with the design matrices at plant, genotype and
#' population levels.
#' \code{ed}, a numeric vector with the estimated effective dimension
#' (or effective degrees of freedom) for each random component of the
#' model (intercept, slope and non-linear trend) at each level of the
#' hierarchy (population, genotype and plant)
#' \code{tot_ed}, a numeric value with the sum of the effective
#' dimensions for all components of the model.
#' \code{vc}, a numeric vector with the (REML) variance component
#' estimates for each random component of the model (intercept,
#' slope and non-linear trend) at each level of the hierarchy
#' (population, genotype and plant)
#' \code{phi}, a numeric value with the error variance estimate.
#' \code{coeff}, a numeric vector with the estimated fixed and random
#' effect coefficients.
#' \code{popLevel}, a data.frame with the estimated population trajectories
#' and first and second order derivatives.
#' \code{genoLevel}, a data.frame with the estimated genotype-specific
#' deviations and trajectories, and their respective first and second
#' order derivatives.
#' \code{plotLevel}, a data.frame with the estimated plant-specific
#' deviations and trajectories, and their respective first and second
#' order derivatives.
#' \code{deviance}, the (REML) deviance at convergence.
#' \code{convergence}, a logical value indicating whether the algorithm
#' managed to converge before the given number of iterations.
#' \code{dim}, a numeric vector with the (model) dimension of each
#' model component (fixed and/or random) at each level of the
#' hierarchy (population, genotype, and plant).
#' These values correspond to the number of parameters to be estimated.
#' \code{family}, an object of class family specifying the distribution
#' and link function.
#' \code{Vp}, the variance-covariance matrix for the coefficients.
#' \code{smooth}, a list with the information about number of segments
#' (nseg), degree of the B-spline basis (bdeg) and penalty order (pord)
#' used for the three levels of the hierarchy.
#'
#' @examples
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#' head(spatCorrectedArch)
#' ggplot2::ggplot(data = spatCorrectedArch,
#'                 ggplot2::aes(x= timeNumber, y = LeafArea_corr, group = plotId)) +
#'   ggplot2::geom_line() +
#'   ggplot2::facet_grid(~geno.decomp)
#'
#' ## We need to specify the genotype-by-treatment interaction
#' ## Treatment: water regime (WW, WD)
#' spatCorrectedArch$treat <- factor(spatCorrectedArch$geno.decomp,
#'                                   labels = substr(levels(spatCorrectedArch$geno.decomp), 1, 2))
#' spatCorrectedArch$genobytreat <- paste0(spatCorrectedArch$genotype,"_",spatCorrectedArch$treat)
#'
#' ## Fit P-Splines Hierarchical Curve Data Model
#' fit.psHDM  <- fitSplineHDM(inDat = spatCorrectedArch,
#'                           trait = "LeafArea_corr",
#'                           time = "timeNumber",
#'                           pop = "geno.decomp",
#'                           genotype = "genobytreat",
#'                           plotId = "plotId",
#'                           difVar = list(geno = FALSE, plant = FALSE),
#'                           smoothPop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothGeno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothPlot = list(nseg = 4, bdeg = 3, pord = 2),
#'                           weights = "wt")
#'
#' ## Visualize the data.frames with predicted values at the three levels of the hierarchy.
#' ## Population level
#' head(fit.psHDM$popLevel)
#' ## Genotype level
#' head(fit.psHDM$genoLevel)
#' ## Plot level
#' head(fit.psHDM$plotLevel)
#'
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#' ## Plots at plant level for some genotypes (as illustration)
#' plot(x = fit.psHDM,
#'     genotypes = c("GenoA14_WD", "GenoA51_WD", "GenoB11_WW",
#'                  "GenoB02_WD","GenoB02_WW"),
#'     themeHDM = themeHDM())
#'
#' @references Pérez-Valencia, D.M., Rodríguez-Álvarez, M.X., Boer, M.P. et al.
#' A two-stage approach for the spatio-temporal analysis of high-throughput
#' phenotyping data. Sci Rep 12, 3177 (2022). \doi{10.1038/s41598-022-06935-9}
#'
#' @export
fitSplineHDM <- function(inDat,
                         trait,
                         time,
                         pop,
                         genotype,
                         plotId,
                         difVar = list(geno = FALSE, plot = FALSE),
                         smoothPop = list(nseg = 10, bdeg = 3, pord = 2),
                         smoothGeno = list(nseg = 10, bdeg = 3, pord = 2),
                         smoothPlot = list(nseg = 10, bdeg = 3, pord = 2),
                         offset = NULL,
                         weights = NULL,
                         family = gaussian(),
                         maxit = 200,
                         trace = TRUE,
                         thr = 1e-03) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(inDat, "data.frame")) {
    stop("inDat should be a data.frame.\n")
  }
  corrCols <- c(genotype, trait, time, pop, plotId)
  if (!all(hasName(x = inDat, name = corrCols))) {
    stop("inDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.\n")
  }
  if (!is.numeric(thr) || length(thr) > 1 || thr < 0) {
    stop("thr should be a positive numerical value.\n")
  }
  ## Unused levels might cause strange behaviour.
  inDat <- droplevels(inDat)
  ## Create a full data set of observations for all combinations of
  ## timepoints, row and column.
  timeDat <- data.frame(time = sort(unique(inDat[[time]])))
  colnames(timeDat) <- time
  fullGrid <- merge(unique(inDat[c(pop, genotype, plotId, "colId", "rowId")]),
                    timeDat)
  inDat <- merge(fullGrid, inDat, all.x = TRUE)
  ## Order the data
  inDat <- inDat[order(inDat[, pop], inDat[, genotype], inDat[, plotId],
                       inDat[, time]), ]
  ## Normalize time.
  raw.time <- timeDat[[time]]
  inDat[[time]] <- inDat[[time]] - min(inDat[[time]]) + 1
  ## Define offset.
  if (is.null(offset)) {
    inDat$offset <- 0
  }
  ## Specify weights.
  if (is.null(weights)) {
    weights <- rep(1, nrow(inDat))
  } else{
    weights <- inDat[[weights]] ^ 2
  }
  ## Convert to factors - only if not already factors.
  for (facVar in c(pop, genotype, plotId)) {
    if (!is.factor(inDat[[facVar]])) {
      inDat[[facVar]] <- factor(inDat[[facVar]])
    }
  }
  ## Elements and number of elements by level of the hierarchy.
  popLevs <- unique(inDat[[pop]])
  genoLevs <- unique(inDat[[genotype]])
  plotLevs <- unique(inDat[[plotId]])
  nPop <- nlevels(inDat[[pop]])
  nPlotPop <- apply(X = table(inDat[[pop]], inDat[[plotId]]),
                    MARGIN = 1, FUN = function(x) { sum(x!=0) })
  nGenoPop <- apply(X = table(inDat[[pop]], inDat[[genotype]]),
                    MARGIN = 1, FUN = function(x) { sum(x!=0) })
  nGeno <- nlevels(inDat[[genotype]])
  nPlotGeno <- apply(table(inDat[[genotype]], inDat[[plotId]]),
                     MARGIN = 1, FUN = function(x) { sum(x!=0) })[genoLevs]
  nTot <- nlevels(inDat[[plotId]])
  # Time interval
  x <- sort(unique(inDat[[time]]))
  ## Construct design matrices: data in an array
  ## Design matrices for the population curves.
  MMPop <- MM.basis(x = x, ndx = smoothPop$nseg, bdeg = smoothPop$bdeg,
                    pord = smoothPop$pord)
  XPop <- MMPop$X
  ZPop <- MMPop$Z
  #$ Design matrices for the genotype curves.
  MMGeno <- MM.basis(x = x, ndx = smoothGeno$nseg, bdeg = smoothGeno$bdeg,
                     pord = smoothGeno$pord)
  XGeno <- MMGeno$X
  ZGeno <- MMGeno$Z
  ## Design matrices for the individual curves.
  MMPlot <- MM.basis(x = x, ndx = smoothPlot$nseg, bdeg = smoothPlot$bdeg,
                     pord = smoothPlot$pord)
  XPlot <- MMPlot$X
  ZPlot <- MMPlot$Z
  ## Response, offset and weights
  y <- inDat[[trait]]
  offset.f <- inDat[["offset"]]
  weights.f <- weights
  ## Set missing values and corresponding weights to 0.
  nas <- is.na(y)
  y[nas] <- 0
  weights.f[nas] <- 0
  ## Construct matrices assigning plants to populations and genotypes
  ## Matrix to assign plants to populations
  plotPlotPop <- rep(1:nPop, nPlotPop)
  if (nPop == 1) {
    CPop <- matrix(rep(1, nTot), ncol = 1)
  } else {
    xxt <- data.frame(plotPop = as.factor(plotPlotPop))
    mft <- model.frame(~ plotPop - 1, data = xxt, drop.unused.levels = TRUE)
    mtt <- terms(mft)
    fTermst <- attr(mtt, "term.labels")[attr(mtt,"dataClasses") == "factor"]
    CPop <- Matrix::sparse.model.matrix(mtt, data = mft,
                                        contrasts.arg = lapply(X = mft[, fTermst, drop = FALSE],
                                                               FUN = contrasts,
                                                               contrasts = FALSE))
    attr(mtt, "xlev") <- .getXlevels(mtt, mft)
    attr(mtt, "contrast") <- attr(CPop, "contrast")
  }
  ## Matrix to assign plants to genotypes.
  plotPlotGeno <- rep(1:nGeno, nPlotGeno)
  if (nGeno == 1) {
    CGeno <- matrix(rep(1, nTot), ncol = 1)
  } else {
    xxg <- data.frame(plotGeno = as.factor(plotPlotGeno))
    mfg <- model.frame(~ plotGeno - 1, xxg, drop.unused.levels = TRUE)
    mtg <- terms(mfg)
    fTermsg <- attr(mtg, "term.labels")[attr(mtg,"dataClasses") == "factor"]
    CGeno <- Matrix::sparse.model.matrix(mtg, data = mfg,
                                         contrasts.arg = lapply(X = mfg[, fTermsg, drop = FALSE],
                                                                FUN = contrasts,
                                                                contrasts = FALSE))
    attr(mtg, "xlev") <- .getXlevels(mtg, mfg)
    attr(mtg, "contrast") <- attr(CGeno,"contrast")
  }
  ## GLAM matrices
  GLAM <- TRUE
  X <- list(X1 = CPop,
            X2 = XPop) # Parametric population effect
  Z <- list(Z1 = list(Z11 = CPop,
                      Z12 = ZPop),
            Z2 = list(Z21 = CGeno,
                      Z22 = XGeno), # Random intercept and slope
            Z3 = list(Z31 = CGeno,
                      Z32 = ZGeno), # Genotype specific curves
            Z4 = list(Z41 = Matrix::Diagonal(n = nTot),
                      Z42 = XPlot), # Random intercept and slopes (plot)
            Z5 = list(Z51 = Matrix::Diagonal(n = nTot),
                      Z52 = ZPlot)) # Genotype specific curves (plot)
  ## Number of parameters: fixed and random (for each component)
  np <- c(ncol(XPop) * nPop, ncol(ZPop) * nPop,
          ncol(XGeno) * nGeno, ncol(ZGeno) * nGeno,
          ncol(XPlot) * nTot, ncol(ZPlot) * nTot)
  names(np) <- c("pop.fixed", "pop.random", "geno.x.random", "geno.z.random",
                 "plot.x.random", "plot.z.random")
  ## Construct precision matrix
  ## Smooth main effect: one per population
  g <- rep(x = list(MMPop$d), times = nPop)
  if (isTRUE(difVar$geno)) {
    ## Random intercepts and slopes (genotype).
    for (i in 1:nPop) {
      g[[nPop + i]] <- lapply(X = 1:smoothGeno$pord, FUN = function(j) {
        rep(x = sapply(X = 1:smoothGeno$pord, FUN = function(k) {
          ifelse(k == j, 1, 0)
        }), times = nGenoPop[i])
      })
    }
    ## Smooth effects (genotype).
    for (i in 1:nPop) {
      g[[nPop * 2 + i]] <- rep(x = MMGeno$d, times = nGenoPop[i])
    }
    if (isTRUE(difVar$plot)) {
      ## Random intercepts and slopes (individual).
      for (i in 1:nGeno) {
        g[[nPop * 3 + i]] <- lapply(X = 1:smoothPlot$pord, FUN = function(j) {
          rep(x = sapply(X = 1:smoothPlot$pord, FUN = function(k) {
            ifelse(k == j, 1, 0)
          }), times = nPlotGeno[i])
        })
      }
      ## Smooth effects (individual).
      for (i in 1:nGeno) {
        g[[nPop * 3 + nGeno + i]] <- rep(x = MMPlot$d, times = nPlotGeno[i])
      }
    } else {
      ## Random intercepts and slopes (individual).
      g[[nPop * 3 + 1]] <- lapply(X = 1:smoothPlot$pord, FUN = function(j) {
        rep(x = sapply(X = 1:smoothPlot$pord, FUN = function(k) {
          ifelse(k == j, 1, 0)
        }), times = nTot)
      })
      ## Smooth effects (individual).
      g[[nPop * 3 + 2]] <- list(rep(x = MMPlot$d, times = nTot))
    }
  } else if (isFALSE(difVar$geno)) {
    ## Random intercepts and slopes (genotype).
    g[[nPop + 1]] <- lapply(X = 1:smoothGeno$pord, FUN = function(j) {
      rep(x = sapply(X = 1:smoothGeno$pord, FUN = function(k) {
        ifelse(k == j, 1, 0)
      }), times = nGeno)
    })
    # Smooth effects (genotype)
    g[[nPop + 2]] <- list(rep(x = MMGeno$d, times = nGeno))
    if (isTRUE(difVar$plot)) {
      ## Random intercepts and slopes (individual).
      for (i in 1:nGeno) {
        g[[nPop + 2 + i]] <- lapply(X = 1:smoothPlot$pord, FUN = function(j) {
          rep(x = sapply(X = 1:smoothPlot$pord, FUN = function(k) {
            ifelse(k == j, 1, 0)
          }), times = nPlotGeno[i])
        })
      }
      ## Smooth effects (individual).
      for (i in 1:nGeno) {
        g[[nPop + 2 + nGeno + i]] <- rep(x = MMPlot$d, times = nPlotGeno[i])
      }
    } else {
      ## Random intercepts and slopes (individual).
      g[[nPop + 3]] <- lapply(X = 1:smoothPlot$pord, FUN = function(j) {
        rep(x = sapply(X = 1:smoothPlot$pord, FUN = function(k) {
          ifelse(k == j, 1, 0)
        }), times = nTot)
      })
      ## Smooth effects (individual).
      g[[nPop + 4]] <- list(rep(x = MMPlot$d, times = nTot))
    }
  }
  ## Construct the components of the precision matrix (as needed by the algorithm).
  g <- constructCapitalLambda(g)
  ## Initialise the parameters.
  la = rep(x = 1, l = length(g) + 1)
  devold = 1e10
  mustart <- etastart <- NULL
  nobs <- length(y)
  eval(family$initialize)
  mu <- mustart
  eta <- family$linkfun(mustart)
  ## Iteration process to estimate coefficients and variance components.
  for (iter in 1:maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offset.f) + (y - mu) / deriv
    w <- as.vector(deriv ^ 2 / family$variance(mu))
    w <- w * weights.f
    mat <- construct.matrices(X = X, Z = Z, z = z, w = w, GLAM = GLAM)
    V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., mat$ZtZ.) # Does not change from iteration to iteration
    ## V is an object of class Matrix, transform it to spam.
    V <- spam::as.spam.dgCMatrix(V)
    ## number of coef for each random component.
    gExt <- lapply(X = g, FUN = function(x) { c(rep(0, np[1]), x) })
    lGExt <- lapply(X = gExt, FUN = spam::diag.spam)
    lP <- append(list(V), lGExt)
    ADcholC <- LMMsolver:::ADchol(lP)
    EDmax <- sapply(g, function(x) sum(abs(x) > 1.0e-10))
    ## First iteration (cholesky decomposition).
    Ginv <- vector(length = length(g[[1]]))
    for (i in 1:length(g)) {
      Ginv <- Ginv + g[[i]]
    }
    G <- 1 / Ginv
    D <- spam::diag.spam(c(rep(x = 0, times = np[1]), Ginv))
    cholHn <- chol(V + D) # The sparsness structure does not change (Ginv is diagonal)
    for (it in 1:maxit) {
      Ginv <- vector(length = length(g[[1]]))
      for (i in 1:length(g)) {
        Ginv <- Ginv + (1 / la[i + 1]) * g[[i]]
      }
      G <- 1 / Ginv
      D <- spam::diag.spam(c(rep(x = 0, times = np[1]), Ginv))
      Hn <- (1 / la[1]) * V + D
      cholHn <- update(cholHn, Hn)
      b <- spam::backsolve(cholHn,
                           spam::forwardsolve(cholHn, (1 / la[1]) * mat$u))
      theta <- 1 / la
      EDc <- theta * LMMsolver:::dlogdet(ADcholC, theta)
      ED <- EDmax - EDc[-1]
      ## Fixed and random coefficients.
      b.fixed  <- b[1:np[1]]
      b.random <- b[-(1:np[1])]
      ## Variance components as in SOP.
      ssv <- ed <- tau <- vector(mode = "list", length = length(g))
      for (i in 1:length(g)) {
        ed[[i]] = ED[i]
        ed[[i]] <- ifelse(ed[[i]] <= 1e-10, 1e-10, ed[[i]])
        ssv[[i]] <- sum(b.random ^ 2 * g[[i]])
        tau[[i]] <- ssv[[i]] / ed[[i]]
        tau[[i]] <- ifelse(tau[[i]] <= 1e-10, 1e-10, tau[[i]])
      }
      ssr <- mat$yty. - t(c(b.fixed, b.random)) %*% (2 * mat$u - V %*% b)
      dev <- deviance_spam(C = cholHn, G = spam::diag.spam(G), w = w[w != 0],
                           sigma2 = la[1], ssr = ssr,
                           edf = sum(b.random ^ 2 * Ginv))[1]
      if (family$family %in% c("gaussian", "quasipoisson")) {
        phi <- as.numeric((ssr / (length(z[w != 0]) - sum(unlist(ed)) - np[1])))
      } else {
        phi <- 1
      }
      ## New variance components and convergence check.
      lanew <- c(phi, unlist(tau))
      dla <- abs(devold - dev)
      if (trace) {
        if(it == 1){
          ## Print header.
          ed.names <- c("", "Deviance", paste0("ed.p", 1:nPop),
                        "ed.g.int", "ed.g.slp", "ed.g.smooth",
                        "ed.i.int", "ed.i.slp", "ed.i.smooth")
          cat(sprintf("%12.12s", ed.names), sep = "")
          cat('\n')
        }
        cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
        cat(sprintf("%12.3f", unlist(ed)), sep = "")
        cat('\n')
      }
      if (dla < thr) break
      la <- lanew
      devold <- dev
    }
    etaOld <- eta
    eta <- Xtheta(X, b.fixed) + Ztheta(Z, b.random, np[-1]) + offset.f
    mu <- family$linkinv(eta)
    ## Convergence criterion: linear predictor.
    tol <- sum((eta - etaOld) ^ 2) / sum(eta ^ 2)
    if (tol < 1e-6 || family$family == "gaussian") break
  }
  coeff <- c(b.fixed, b.random)
  ## Plant raw data.
  aux <- matrix(inDat[, trait], ncol = nTot)
  obsPlot <- lapply(X = split(aux, rep(plotPlotGeno, each = nrow(aux))),
                    FUN = function(x, nobs) {
                      matrix(x, nrow = nobs)
                    }, nobs = length(x))
  names(obsPlot) <- genoLevs
  ## Object to be returned
  res <- structure(
    list(y = obsPlot,
         time = raw.time,
         genoLevs = genoLevs,
         popLevs = popLevs,
         plotLevs = plotLevs,
         nPlotPop = nPlotPop,
         nGenoPop = nGenoPop,
         nPlotGeno = nPlotGeno,
         MM = list(MMPop = MMPop, MMGeno = MMGeno, MMPlot = MMPlot),
         ed = unlist(ed),
         tot_ed = sum(unlist(ed)),
         vc = unlist(tau),
         phi = phi,
         coeff = coeff,
         deviance = dev,
         convergence = it < maxit,
         dim = np,
         family = family,
         Vp = spam::chol2inv(cholHn),
         smooth = list(smoothPop = smoothPop, smoothGeno = smoothGeno,
                       smoothPlot = smoothPlot)),
    class = c("psHDM", "list"),
    trait = trait
  )
  ## Make predictions on original data points.
  preds <- predict.psHDM(res,
                         se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                         trace = FALSE)
  ## Add predictions to results.
  res <- c(res, preds[c("popLevel", "genoLevel", "plotLevel")])
  class(res) <- c("psHDM", "list")
  attr(res, which = "trait") <- trait
  return(res)
}

