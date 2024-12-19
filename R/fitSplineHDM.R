#' Fit P-Spline Hierarchical Curve Data Models
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
#' desired. The function outputs are estimated curves (time series of trajectories
#' and deviations) and their first and second derivatives for the three-levels
#' of the hierarchy. The outputs can then be used to estimate relevant parameters
#' from the curves for further analysis (see \code{\link{estimateSplineParameters}}).
#'
#' @param inDat A data.frame with corrected spatial data.
#' @param genotypes A character vector indicating the genotypes for which
#' hierarchical models should be fitted. If \code{NULL}, splines will be fitted
#' for all genotypes.
#' @param plotIds A character vector indicating the plotIds for which
#' hierarchical models should be fitted. If \code{NULL}, splines will be
#' fitted for all plotIds.
#' @param trait A character string indicating the trait for which the spline
#' should be fitted.
#' @param useTimeNumber Should the timeNumber be used instead of the timePoint?.
#' If \code{useTimeNumber = FALSE}, inDat should contain a column called timePoint
#' of class \code{POSIXct}.
#' @param timeNumber If \code{useTimeNumber = TRUE}, a character vector
#' indicating the column containing the numerical time to use.
#' @param pop A character string indicating the the populations to which each
#' genotype/variety belongs. This variable must be a factor in the data frame.
#' @param genotype A character string indicating the populations to which each
#' genotype/variety belongs. This variable must be a factor in the data frame.
#' @param plotId A character string indicating the genotypes/varieties to which
#' each plant/plot/individual belongs. This variable must be a factor in the
#' data frame.
#' @param weights A character string indicating the column in the data containing
#' the weights to be used in the fitting process (for error propagation from
#' first stage to second stage). By default, when \code{weights = NULL}, the
#' weights are considered to be one.
#' @param difVar Should different variances for random effects at genotype
#' (separately for each population) and plant level (separately for each
#' genotype) be considered?.
#' @param smoothPop A list specifying the P-Spline model at the population
#' level (nseg: number of segments; bdeg: degree of the B-spline basis; pord:
#' penalty order).
#' @param smoothGeno A list specifying the P-Spline model at the genotype
#' level.
#' @param smoothPlot A list specifying the P-Spline model at the plant level.
#' @param offset A character string indicating the column in the data with
#' an a priori known component to be included in the linear predictor during
#' fitting. By default, when \code{offset = NULL}, the offset is considered to
#' be zero.
#' @param family An object of class \code{family} specifying the distribution
#' and link function. The default is \code{gaussian()}.
#' @param maxit An optional value that controls the maximum number of iterations
#' of the algorithm. The default is 200.
#' @param trace An optional value that controls the function trace.
#' The default is \code{TRUE}.
#' @param thr An optional value that controls the convergence threshold of the
#' algorithm. The default is 1.e-03.
#' @param minNoTP The minimum number of time points for which data should be
#' available for a plant. Defaults to 60% of all time points present in the
#' TP object. No splines are fitted for plants with less than the minimum number
#' of timepoints.
#'
#' @returns An object of class \code{psHDM}, a list with the following outputs:
#' \code{time}, a numeric vector with the timepoints.
#' \code{popLevs}, a data.frame with the names of the populations
#' \code{genoLevs}, a factor with the names of the genotypes.
#' \code{plotLevs}, a factor with the names of the plants
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
#' \code{cholHn}, the inverse of the variance-covariance matrix for the
#' coefficients.
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
#'   ggplot2::geom_line(na.rm = TRUE) +
#'   ggplot2::facet_grid(~geno.decomp)
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
#'                           useTimeNumber = TRUE,
#'                           timeNumber = "timeNumber",
#'                           genotypes = c("GenoA14_WD", "GenoA51_WD",
#'                                        "GenoB11_WW", "GenoB02_WD",
#'                                        "GenoB02_WW"),
#'                           pop = "geno.decomp",
#'                           genotype = "genoTreat",
#'                           plotId = "plotId",
#'                           weights = "wt",
#'                           difVar = list(geno = FALSE, plot = FALSE),
#'                           smoothPop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothGeno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothPlot = list(nseg = 4, bdeg = 3, pord = 2),
#'                           trace = FALSE)
#'
#' ## Visualize the data.frames with predicted values at the three levels of
#' ## the hierarchy.
#'
#' # Population level
#' head(fit.psHDM$popLevel)
#'
#' # Genotype level
#' head(fit.psHDM$genoLevel)
#'
#' # Plot level
#' head(fit.psHDM$plotLevel)
#'
#' @references Pérez-Valencia, D.M., Rodríguez-Álvarez, M.X., Boer, M.P. et al.
#' A two-stage approach for the spatio-temporal analysis of high-throughput
#' phenotyping data. Sci Rep 12, 3177 (2022). \doi{10.1038/s41598-022-06935-9}
#'
#' @family functions for fitting hierarchical curve data models
#'
#' @export
fitSplineHDM <- function(inDat,
                         genotypes = NULL,
                         plotIds = NULL,
                         trait,
                         useTimeNumber = FALSE,
                         timeNumber = NULL,
                         pop = "pop",
                         genotype = "genotype",
                         plotId = "plotId",
                         weights = NULL,
                         difVar = list(geno = FALSE, plot = FALSE),
                         smoothPop = list(nseg = 10, bdeg = 3, pord = 2),
                         smoothGeno = list(nseg = 10, bdeg = 3, pord = 2),
                         smoothPlot = list(nseg = 10, bdeg = 3, pord = 2),
                         offset = NULL,
                         family = gaussian(),
                         maxit = 200,
                         trace = TRUE,
                         thr = 1e-03,
                         minNoTP = NULL) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(inDat, "data.frame")) {
    stop("inDat should be a data.frame.\n")
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.\n")
  }
  if (!is.numeric(thr) || length(thr) > 1 || thr < 0) {
    stop("thr should be a positive numerical value.\n")
  }
  if (isTRUE(useTimeNumber) &&
      (is.null(timeNumber) || !is.character(timeNumber) ||
       length(timeNumber) > 1)) {
    stop("timeNumber should be a character string of length 1.\n")
  }
  corrCols <- c(trait, if (useTimeNumber) timeNumber else "timePoint",
                genotype, pop, plotId, "colId", "rowId")
  if (!all(hasName(x = inDat, name = corrCols))) {
    stop("inDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!is.null(genotypes) &&
      (!is.character(genotypes) ||
       !all(genotypes %in% inDat[[genotype]]))) {
    stop("genotypes should be a character vector of genotypes in inDat.\n")
  }
  if (!is.null(plotIds) &&
      (!is.character(plotIds) || !all(plotIds %in% inDat[[plotId]]))) {
    stop("plotIds should be a character vector of plotIds in inDat.\n")
  }
  if (!is.null(minNoTP) && (!is.numeric(minNoTP) || length(minNoTP) > 1)) {
    stop("minNoTP should be a numerical value.\n")
  }
  if (!useTimeNumber) {
    if (!inherits(inDat[["timePoint"]], "POSIXct")) {
      stop("Column timePoint should be of class POSIXct.\n")
    }
    ## Convert time point to time number with the first time point as 0.
    minTime <- min(inDat[["timePoint"]], na.rm = TRUE)
    inDat[["timeNumber"]] <- as.numeric(inDat[["timePoint"]] - minTime) / 1000
  } else {
    if (!is.numeric(inDat[[timeNumber]])) {
      stop("timeNumber should be a numerical column.\n")
    }
    inDat[["timeNumber"]] <- inDat[[timeNumber]]
  }
  ## Restrict inDat to selected genotypes and plotIds.
  if (!is.null(genotypes)) {
    inDat <- inDat[inDat[[genotype]] %in% genotypes, ]
  }
  if (!is.null(plotIds)) {
    inDat <- inDat[inDat[["plotId"]] %in% plotIds, ]
  }
  if (nrow(inDat) == 0) {
    stop("At least one valid combination of genotype and plotId should be ",
         "selected.\n")
  }
  if (!is.list(difVar) || length(difVar) != 2 ||
      !(setequal(names(difVar), c("geno", "plot")))) {
    stop("difVar should be a named list of length 2.\n")
  }
  chkSmooth(smoothPop)
  chkSmooth(smoothGeno)
  chkSmooth(smoothPlot)
  ## Unused levels might cause strange behaviour.
  inDat <- droplevels(inDat)
  ## Check that pop - geno - plot structure is unambiguously defined.
  genoPopTab <- table(inDat[[genotype]], inDat[[pop]])
  genoPopCount <- rowSums(genoPopTab > 0)
  dupGeno <- names(genoPopCount[genoPopCount > 1])
  if (length(dupGeno) > 0) {
    stop("The following genotypes are in multiple populations:\n",
         paste(dupGeno, collapse = ", "), "\n")
  }
  plotGenoTab <- table(inDat[[plotId]], inDat[[genotype]])
  plotGenoCount <- rowSums(plotGenoTab > 0)
  dupPlot <- names(plotGenoCount[plotGenoCount > 1])
  if (length(dupPlot) > 0) {
    stop("The following plots are specified for multiple genotypes:\n",
         paste(dupPlot, collapse = ", "), "\n")
  }
  ## Determine minimum number of time points required.
  nTimeNumber <- length(unique(inDat[["timeNumber"]]))
  if (is.null(minNoTP)) {
    minTP <- 0.6 * nTimeNumber
  } else if (minNoTP < 0 || minNoTP > nTimeNumber) {
    stop("minNoTP should be a number bewtween 0 and ", nTimeNumber, ".\n")
  } else {
    minTP <- minNoTP
  }
  ## Check for plotIds that have a limited amount of observations.
  plotTab <- table(inDat[!is.na(inDat[[trait]]), "plotId"])
  plotLimObs <- names(plotTab[plotTab < minTP])
  if (length(plotLimObs) > 5) {
    warning("More than 5 plots have observations for less than the ",
            "minimum number of time points, which is ", round(minTP), ". The  ",
            "first 5 are printed, to see them all run attr(..., 'plotLimObs') ",
            "on the output\n",
            paste(plotLimObs[1:5], collapse = ", "), "\n", call. = FALSE)
  } else if (length(plotLimObs) > 0) {
    warning("The following plots have observations for less than ",
            "the minimum number of time points, which is ", round(minTP), ":\n",
            paste(plotLimObs, collapse = ", "), "\n", call. = FALSE)
  }
  ## Create data.frame with time number and, if present,
  ## time point on prediction scale.
  timeRange <- data.frame(timeNumber = sort(unique(inDat[["timeNumber"]])))
  if (hasName(x = inDat, name = "timePoint")) {
    timeRange[["timePoint"]] <- sort(unique(inDat[["timePoint"]]))
  }
  ## Create a full data set of observations for all combinations of timepoints.
  fullGrid <- merge(unique(inDat[c(pop, genotype, plotId)]),
                    timeRange)
  inDat <- merge(fullGrid, inDat, all.x = TRUE)
  ## Order the data
  inDat <- inDat[order(inDat[, pop], inDat[, genotype], inDat[, plotId],
                       inDat[, "timeNumber"]), ]
  ## Normalize time.
  rawTime <- timeRange[["timeNumber"]]
  inDat[["timeNumber"]] <- inDat[["timeNumber"]] - min(inDat[["timeNumber"]]) + 1
  ## Define offset.
  if (is.null(offset)) {
    offsetf <- rep(0, nrow(inDat))
  } else{
    if (!hasName(x = inDat, name = offset)) {
      stop("offset should be a column in inDat.\n")
    }
    offsetf <- inDat[[offset]]
  }
  ## Specify weights.
  if (is.null(weights)) {
    weightsf <- rep(1, nrow(inDat))
  } else{
    if (!hasName(x = inDat, name = weights)) {
      stop("weights should be a column in inDat.\n")
    }
    weightsf <- inDat[[weights]]
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
  x <- sort(unique(inDat[["timeNumber"]]))
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
  ## Response.
  y <- inDat[[trait]]
  ## Set missing values and corresponding weights to 0.
  nas <- is.na(y)
  y[nas] <- 0
  weightsf[nas] <- 0
  offsetf[nas] <- 0
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
    g <- c(g, lapply(X = nGenoPop, FUN = constructG, ord = smoothGeno$pord))
    ## Smooth effects (genotype).
    for (genoPop in nGenoPop) {
      g <- c(g, list(rep(x = MMGeno$d, times = genoPop)))
    }
    if (isTRUE(difVar$plot)) {
      ## Random intercepts and slopes (individual).
      g <- c(g, lapply(X = nPlotGeno, FUN = constructG, ord = smoothPlot$pord))
      ## Smooth effects (individual).
      for (plotGeno in nPlotGeno) {
        g <- c(g, list(rep(x = MMGeno$d, times = plotGeno)))
      }
    } else {
      ## Random intercepts and slopes (individual).
      g <- c(g, lapply(X = nTot, FUN = constructG, ord = smoothPlot$pord))
      ## Smooth effects (individual).
      g <- c(g, list(rep(x = MMPlot$d, times = nTot)))
    }
  } else if (isFALSE(difVar$geno)) {
    ## Random intercepts and slopes (genotype).
    g <- c(g, lapply(X = nGeno, FUN = constructG, ord = smoothGeno$pord))
    # Smooth effects (genotype)
    g <- c(g, list(rep(x = MMGeno$d, times = nGeno)))
    if (isTRUE(difVar$plot)) {
      ## Random intercepts and slopes (individual).
      g <- c(g, lapply(X = nPlotGeno, FUN = constructG, ord = smoothPlot$pord))
      ## Smooth effects (individual).
      for (plotGeno in nPlotGeno) {
        g <- c(g, list(rep(x = MMPlot$d, times = plotGeno)))
      }
    } else {
      ## Random intercepts and slopes (individual).
      g <- c(g, lapply(X = nTot, FUN = constructG, ord = smoothPlot$pord))
      ## Smooth effects (individual).
      g <- c(g, list(rep(x = MMPlot$d, times = nTot)))
    }
  }
  ## Construct the components of the precision matrix (as needed by the algorithm).
  g <- constructCapitalLambda(g)
  ## Construct names for ed.
  edNames <- c(paste0("p", 1:nPop),
               if (isTRUE(difVar$geno))
                 paste0(paste0("g.", c("int", "slp", "smooth"),
                               rep(1:nPop, each = 3))) else
                                 c("g.int", "g.slp", "g.smooth"),
               if (isTRUE(difVar$plot))
                 paste0(paste0("i.", c("int", "slp", "smooth"),
                               rep(1:nGeno, each = 3))) else
                                 c("i.int", "i.slp", "i.smooth"))
  ## Initialise the parameters.
  la <- rep(x = 1, times = nrow(g) + 1)
  devold <- 1e10
  mustart <- etastart <- NULL
  nobs <- length(y)
  eval(family$initialize)
  mu <- mustart
  eta <- family$linkfun(mustart)
  ## Iteration process to estimate coefficients and variance components.
  for (iter in 1:maxit) {
    deriv <- family$mu.eta(eta)
    z <- (eta - offsetf) + (y - mu) / deriv
    w <- as.vector(deriv ^ 2 / family$variance(mu))
    w <- w * weightsf
    mat <- construct.matrices(X = X, Z = Z, z = z, w = w, GLAM = GLAM)
    V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., mat$ZtZ.) # Does not change from iteration to iteration
    ## V is an object of class Matrix, transform it to spam.
    V <- spam::as.spam.dgCMatrix(V)
    ## number of coef for each random component.
    lGExt <- apply(X = g, MARGIN = 1, FUN = function(x) {
      spam::diag.spam(c(rep(0, np[1]), x))
    }, simplify = FALSE)
    lP <- append(list(V), lGExt)
    ADcholC <- LMMsolver:::ADchol(lP)
    EDmax <- rowSums(g > 1e-10)
    ## First iteration (cholesky decomposition).
    Ginv <- colSums(g)
    G <- 1 / Ginv
    D <- spam::diag.spam(c(rep(x = 0, times = np[1]), Ginv))
    cholHn <- chol(V + D) # The sparsness structure does not change (Ginv is diagonal)
    for (it in 1:maxit) {
      Ginv <- colSums(g / la[-1])
      G <- 1 / Ginv
      D <- spam::diag.spam(c(rep(x = 0, times = np[1]), Ginv))
      Hn <- V / la[1] + D
      cholHn <- update(cholHn, Hn)
      b <- spam::backsolve(cholHn,
                           spam::forwardsolve(cholHn, mat$u / la[1]))
      theta <- 1 / la
      EDc <- theta * LMMsolver:::dlogdet(ADcholC, theta)
      ## Fixed and random coefficients.
      bFixed <- b[1:np[1]]
      bRandom <- b[-(1:np[1])]
      bRandom2 <- bRandom ^ 2
      ## Variance components as in SOP.
      ed <- ifelse(EDmax - EDc[-1] <= 1e-10, 1e-10, EDmax - EDc[-1])
      ssv <- g %*% bRandom2
      tau <- ifelse(ssv / ed <= 1e-10, 1e-10, ssv / ed)
      ssr <- mat$yty. - crossprod(c(bFixed, bRandom), 2 * mat$u - V %*% b)
      dev <- deviance_spam(C = cholHn, G = spam::diag.spam(G), w = w[w != 0],
                           sigma2 = la[1], ssr = ssr,
                           edf = sum(bRandom2 * Ginv))[1]
      if (family$family %in% c("gaussian", "quasipoisson")) {
        phi <- as.numeric((ssr / (length(z[w != 0]) - sum(ed) - np[1])))
      } else {
        phi <- 1
      }
      ## New variance components and convergence check.
      lanew <- c(phi, tau)
      dla <- abs(devold - dev)
      if (trace) {
        if(it == 1){
          ## Print header.
          cat("Effective dimensions\n")
          cat("-------------------------\n")
          cat(sprintf("%1$3s %2$12s","It.","Deviance"), sep = "")
          cat(sprintf("%10s", edNames), sep = "")
          cat("\n")
        }
        cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
        cat(sprintf("%10.3f", ed), sep = "")
        cat('\n')
      }
      if (dla < thr) break
      la <- lanew
      devold <- dev
    }
    etaOld <- eta
    eta <- Xtheta(X, bFixed) + Ztheta(Z, bRandom, np[-1]) + offsetf
    mu <- family$linkinv(eta)
    ## Convergence criterion: linear predictor.
    tol <- sum((eta - etaOld) ^ 2) / sum(eta ^ 2)
    if (tol < 1e-6 || family$family == "gaussian") break
  }
  coeff <- c(bFixed, bRandom)
  ## Plant raw data.
  aux <- matrix(inDat[, trait], ncol = nTot)
  obsPlot <- lapply(X = split(aux, rep(plotPlotGeno, each = nrow(aux))),
                    FUN = function(x, nobs) {
                      matrix(x, nrow = nobs)
                    }, nobs = length(x))
  names(obsPlot) <- genoLevs
  ## Add names to ed and vc.
  vc <- as.vector(tau)
  names(vc) <- names(ed) <- edNames
  ## Check convergence.
  convergence <- it < maxit
  if (!convergence) {
    warning("Model didn't converge in ", it, " iterations.\n")
  }
  ## Object to be returned.
  res <- structure(
    list(y = obsPlot,
         time = timeRange,
         popLevs = popLevs,
         genoLevs = genoLevs,
         plotLevs = plotLevs,
         nPlotPop = nPlotPop,
         nGenoPop = nGenoPop,
         nPlotGeno = nPlotGeno,
         MM = list(MMPop = MMPop, MMGeno = MMGeno, MMPlot = MMPlot),
         ed = ed,
         vc = vc,
         phi = phi,
         coeff = coeff,
         deviance = dev,
         convergence = convergence,
         dim = np,
         family = family,
         cholHn = cholHn,
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
  attr(res, which = "plotLimObs") <- plotLimObs
  return(res)
}

