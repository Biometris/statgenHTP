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
#' @param geno A character string indicating the populations to which each
#' genotype/variety belongs. This variable must be a factor in the data frame.
#' @param plant A character string indicating the genotypes/varieties to which
#' each plant/plot/individual belongs. This variable must be a factor in the
#' data frame.
#' @param dif.var Should different variances for random effects at genotype
#' (separately for each population) and plant level (separately for each
#' genotype) be considered?.
#' @param smooth.pop A list specifying the P-Spline model at the population
#' level (nseg: number of segments; bdeg: degree of the B-spline basis; pord:
#' penalty order).
#' @param smooth.geno A list specifying the P-Spline model at the genotype
#' level.
#' @param smooth.plant A list specifying the P-Spline model at the plant level.
#' @param offset An optional numerical vector containing an a priori known
#' component to be included in the linear predictor during fitting. The default
#' is \code{NULL}.
#' @param weights A character string indicating the weights to be used in the
#' fitting process (for error propagation from first stage to second stage).
#' By default, the weights are considered to be one. The default is \code{NULL}.
#' @param family An object of class \code{family} specifying the distribution
#' and link function.  The default is \code{gaussian()}.
#' @param maxit An optional value that controls the maximum number of iterations
#' of the algorithm.  The default is 200.
#' @param trace An optional value that controls the function trace.
#' The default is \code{TRUE}.
#' @param thr An optional value that controls the convergence threshold of the
#' algorithm. The default is 1.e-03.
#'
#' @return An object of class \code{psHDM}, a list with the following outputs:
#' \code{time}, a numeric vector with the timepoints.
#' \code{l.geno}, a vector with the names of the genotypes.
#' \code{l.pop}, a vector with the names of the populations
#' \code{l.plant}, a vector with the names of the plants
#' \code{n.plants_p_pop}, a numeric vector with the number of plants per
#' population.
#' \code{n.geno_p_pop}, a numeric vector with the number of genotypes per
#' population.
#' \code{n.plants_p_geno}, a numeric vector with the number of plants per
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
#' \code{pop.level}, a data.frame with the estimated population trajectories
#' and first and second order derivatives.
#' \code{geno.level}, a data.frame with the estimated genotype-specific
#' deviations and trajectories, and their respective first and second
#' order derivatives.
#' \code{plant.level}, a data.frame with the estimated plant-specific
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
#' ggplot(data = spatCorrectedArch, aes(x= timeNumber, y = LeafArea_corr, group = plotId))+
#'  geom_line() + facet_grid(~geno.decomp)
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
#'                           geno = "genobytreat",
#'                           plant = "plotId",
#'                           dif.var = list(geno = FALSE, plant = FALSE),
#'                           smooth.pop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smooth.geno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smooth.plant = list(nseg = 4, bdeg = 3, pord = 2),
#'                           weights = "wt")
#'
#' ## Visualize the data.frames with predicted values at the three levels of the hierarchy.
#' ## Population level
#'   head(fit.psHDM$pop.level)
#' ## Genotype level
#'   head(fit.psHDM$geno.level)
#' ## Plant level
#'   head(fit.psHDM$plant.level)
#'
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#' ## Plots at plant level for some genotypes (as illustration)
#' plot(object = fit.psHDM,
#'     geno.sub = c("GenoA14_WD","GenoA51_WD","GenoB11_WW","GenoB02_WD","GenoB02_WW"),
#'     my.theme = my.theme())
#'
#' @export
fitSplineHDM <- function(inDat,
                         trait,
                         time,
                         pop,
                         geno,
                         plant,
                         dif.var = list(geno = FALSE, plant = FALSE),
                         smooth.pop = list(nseg = 10, bdeg = 3, pord = 2),
                         smooth.geno = list(nseg = 10, bdeg = 3, pord = 2),
                         smooth.plant = list(nseg = 10, bdeg = 3, pord = 2),
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
  corrCols <- c(geno, trait, time, pop, plant)
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
  fullGrid <- merge(unique(inDat[c(pop, geno, plant, "colId", "rowId")]),
                    timeDat)
  inDat <- merge(fullGrid, inDat, all.x = TRUE)
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
  for (facVar in c(pop, geno, plant)) {
    if (!is.factor(inDat[[facVar]])) {
      inDat[[facVar]] <- factor(inDat[[facVar]])
    }
  }
  ## Elements and number of elements by level of the hierarchy.
  l.pop <- levels(inDat[[pop]])
  l.geno <- levels(inDat[[geno]])
  l.plant <- levels(inDat[[plant]])
  n.pop <- nlevels(inDat[[pop]])
  n.plants_p_pop <- apply(X = table(inDat[[pop]], inDat[[plant]]),
                          MARGIN = 1, FUN = function(x) { sum(x!=0) })
  n.geno_p_pop <- apply(X = table(inDat[[pop]], inDat[[geno]]),
                        MARGIN = 1, FUN = function(x) { sum(x!=0) })
  n.geno <- nlevels(inDat[[geno]])
  n.plants_p_geno <- apply(table(inDat[[geno]], inDat[[plant]]),
                           MARGIN = 1, FUN = function(x) { sum(x!=0) })[l.geno]
  n.tot <- nlevels(inDat[[plant]])
  # Time interval
  x <- sort(unique(inDat[[time]]))
  ## Construct design matrices: data in an array
  ## Design matrices for the population curves.
  MM.pop <- MM.basis(x = x, ndx = smooth.pop$nseg, bdeg = smooth.pop$bdeg,
                     pord = smooth.pop$pord)
  X.pop <- MM.pop$X
  Z.pop <- MM.pop$Z
  #$ Design matrices for the genotype curves.
  MM.geno <- MM.basis(x = x, ndx = smooth.geno$nseg, bdeg = smooth.geno$bdeg,
                      pord = smooth.geno$pord)
  X.geno <- MM.geno$X
  Z.geno <- MM.geno$Z
  ## Design matrices for the individual curves.
  MM.ind <- MM.basis(x = x, ndx = smooth.plant$nseg, bdeg = smooth.plant$bdeg,
                     pord = smooth.plant$pord)
  X.ind <- MM.ind$X
  Z.ind <- MM.ind$Z
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
  ind.ind.pop <- rep(1:n.pop, n.plants_p_pop)
  if (n.pop == 1) {
    C.pop <- Matrix::Diagonal(n = n.tot)
  } else {
    xxt <- data.frame(ind.pop = as.factor(ind.ind.pop))
    mft <- model.frame(~ ind.pop - 1, data = xxt, drop.unused.levels = TRUE)
    mtt <- terms(mft)
    f.termst <- attr(mtt, "term.labels")[attr(mtt,"dataClasses") == "factor"]
    C.pop <- Matrix::sparse.model.matrix(mtt, data = mft,
                                         contrasts.arg = lapply(X = mft[, f.termst, drop = FALSE],
                                                                FUN = contrasts,
                                                                contrasts = FALSE))
    attr(mtt, "xlev") <- .getXlevels(mtt, mft)
    attr(mtt, "contrast") <- attr(C.pop, "contrast")
  }
  ## Matrix to assign plants to genotypes.
  ind.ind.geno <- rep(1:n.geno, n.plants_p_geno)
  if (n.geno == 1) {
    C.geno <- Matrix::Diagonal(n = n.tot)
  } else {
    xxg <- data.frame(ind.pop = as.factor(ind.ind.geno))
    mfg <- model.frame(~ ind.pop - 1, xxg, drop.unused.levels = TRUE)
    mtg <- terms(mfg)
    f.termsg <- attr(mtg, "term.labels")[attr(mtg,"dataClasses") == "factor"]
    C.geno <- Matrix::sparse.model.matrix(mtg, data = mfg,
                                          contrasts.arg = lapply(X = mfg[,f.termsg, drop = FALSE],
                                                                 FUN = contrasts,
                                                                 contrasts = FALSE))
    attr(mtg, "xlev") <- .getXlevels(mtg, mfg)
    attr(mtg, "contrast") <- attr(C.geno,"contrast")
  }
  ## GLAM matrices
  GLAM <- TRUE
  X <- list(X1 = C.pop,
            X2 = X.pop) # Parametric population effect
  Z <- list(Z1 = list(Z11 = C.pop,
                      Z12 = Z.pop),
            Z2 = list(Z21 = C.geno,
                      Z22 = X.geno), # Random intercept and slope
            Z3 = list(Z31 = C.geno,
                      Z32 = Z.geno), # Genotype specific curves
            Z4 = list(Z41 = Matrix::Diagonal(n = n.tot),
                      Z42 = X.ind), # Random intercept and slopes (individual)
            Z5 = list(Z51 = Matrix::Diagonal(n = n.tot),
                      Z52 = Z.ind)) # Genotype specific curves (individual)
  ## Number of parameters: fixed and random (for each component)
  np <- c(ncol(X.pop) * n.pop, ncol(Z.pop) * n.pop,
          ncol(X.geno) * n.geno, ncol(Z.geno) * n.geno,
          ncol(X.ind) * n.tot, ncol(Z.ind) * n.tot)
  names(np) <- c("pop.fixed", "pop.random", "geno.x.random", "geno.z.random",
                 "ind.x.random", "ind.z.random")
  np.comp <- c(np[1] + np[2], np[3] + np[4], np[5] + np[6])
  names(np.comp) <- c("pop", "geno", "plant")
  np.e <- cumsum(np.comp)
  np.s <- np.e - np.comp + 1
  ## Construct precision matrix
  ## Smooth main effect: one per population
  g <- rep(x = list(MM.pop$d), times = n.pop)
  if (isTRUE(dif.var$geno)) {
    ## Random intercepts and slopes (genotype).
    for (i in 1:n.pop) {
      g[[n.pop + i]] <- lapply(X = 1:smooth.geno$pord, FUN = function(j) {
        rep(x = sapply(X = 1:smooth.geno$pord, FUN = function(k) {
          ifelse(k == j, 1, 0)
        }), times = n.geno_p_pop[i])
      })
    }
    ## Smooth effects (genotype).
    for (i in 1:n.pop) {
      g[[n.pop * 2 + i]] <- rep(x = MM.geno$d, times = n.geno_p_pop[i])
    }
    if (isTRUE(dif.var$plant)) {
      ## Random intercepts and slopes (individual).
      for (i in 1:n.geno) {
        g[[n.pop * 3 + i]] <- lapply(X = 1:smooth.plant$pord, FUN = function(j) {
          rep(x = sapply(X = 1:smooth.plant$pord, FUN = function(k) {
            ifelse(k == j, 1, 0)
          }), times = n.plants_p_geno[i])
        })
      }
      ## Smooth effects (individual).
      for (i in 1:n.geno) {
        g[[n.pop * 3 + n.geno + i]] <- rep(x = MM.ind$d,
                                           times = n.plants_p_geno[i])
      }
    } else {
      ## Random intercepts and slopes (individual).
      g[[n.pop * 3 + 1]] <- lapply(X = 1:smooth.plant$pord, FUN = function(j) {
        rep(x = sapply(X = 1:smooth.plant$pord, FUN = function(k) {
          ifelse(k == j, 1, 0)
        }), times = n.tot)
      })
      ## Smooth effects (individual).
      g[[n.pop * 3 + 2]] <- list(rep(x = MM.ind$d, times = n.tot))
    }
  } else if (isFALSE(dif.var$geno)) {
    ## Random intercepts and slopes (genotype).
    g[[n.pop + 1]] <- lapply(X = 1:smooth.geno$pord, FUN = function(j) {
      rep(x = sapply(X = 1:smooth.geno$pord, FUN = function(k) {
        ifelse(k == j, 1, 0)
      }), times = n.geno)
    })
    # Smooth effects (genotype)
    g[[n.pop + 2]] <- list(rep(x = MM.geno$d, times = n.geno))
    if (isTRUE(dif.var$plant)) {
      ## Random intercepts and slopes (individual).
      for (i in 1:n.geno) {
        g[[n.pop + 2 + i]] <- lapply(X = 1:smooth.plant$pord, FUN = function(j) {
          rep(x = sapply(X = 1:smooth.plant$pord, FUN = function(k) {
            ifelse(k == j, 1, 0)
          }), times = n.plants_p_geno[i])
        })
      }
      ## Smooth effects (individual).
      for (i in 1:n.geno) {
        g[[n.pop + 2 + n.geno + i]] <- rep(x = MM.ind$d,
                                           times = n.plants_p_geno[i])
      }
    } else {
      ## Random intercepts and slopes (individual).
      g[[n.pop + 3]] <- lapply(X = 1:smooth.plant$pord, FUN = function(j) {
        rep(x = sapply(X = 1:smooth.plant$pord, FUN = function(k) {
          ifelse(k == j, 1, 0)
        }), times = n.tot)
      })
      ## Smooth effects (individual).
      g[[n.pop + 4]] <- list(rep(x = MM.ind$d, times = n.tot))
    }
  }
  ## Construct the components of the precision matrix (as needed by the algorithm).
  g <- construct.capital.lambda(g)
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
    # number of coef for each random component.
    g_ext <- lapply(X = g, FUN = function(x) { c(rep(0, np[1]), x) })
    lG_ext <- lapply(X = g_ext, FUN = spam::diag.spam)
    lP <- append(list(V), lG_ext)
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
        cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
        cat(sprintf("%12.3f", unlist(ed)), sep = "")
        cat('\n')
      }
      if (dla < thr) break
      la <- lanew
      devold <- dev
    }
    eta.old <- eta
    eta <- Xtheta(X, b.fixed) + Ztheta(Z, b.random, np[-1]) + offset.f
    mu <- family$linkinv(eta)
    ## Convergence criterion: linear predictor.
    tol <- sum((eta - eta.old) ^ 2) / sum(eta ^ 2)
    if (tol < 1e-6 || family$family == "gaussian") break
  }
  coeff <- c(b.fixed, b.random)
  ## Plant raw data.
  aux <- matrix(inDat[, trait], ncol = n.tot)
  obs_ind <- lapply(X = split(aux, rep(ind.ind.geno, each = nrow(aux))),
                    FUN = function(x, nobs) {
                      matrix(x, nrow = nobs)
                    }, nobs = length(x))
  names(obs_ind) <- l.geno
  ## Object to be returned
  res <- structure(
    list(y = obs_ind,
         time = raw.time,
         l.geno = l.geno,
         l.pop = l.pop,
         l.plant = l.plant,
         n.plants_p_pop = n.plants_p_pop,
         n.geno_p_pop = n.geno_p_pop,
         n.plants_p_geno = n.plants_p_geno,
         MM = list(MM.pop = MM.pop, MM.geno = MM.geno, MM.plant = MM.ind),
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
         smooth = list(smooth.pop = smooth.pop, smooth.geno = smooth.geno,
                       smooth.plant = smooth.plant)),
    class = c("psHDM", "list")
  )
  ## Make predictions on original data points.
  preds <- predict.psHDM(res,
                         se = list(pop = FALSE, geno = FALSE, plant = FALSE),
                         trace = trace)
  ## Add predictions to results.
  res <- c(res, preds[c("pop.level", "geno.level", "plant.level")])
  class(res) <- c("psHDM", "list")
  return(res)
}

