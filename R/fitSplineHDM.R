##..........................................................
## R-functions required to fit, predict and plot the P-spline Hierarchical Curve Data Model
## used in the second stage of the two-stage approach described in the paper
##
## "A two-stage approach for the spatio-temporal analysis of high-throughput phenotyping data"
##
## Authors: Diana M. Pérez-Valencia, María Xosé Rodríguez-Álvarez, Martin P. Boer, Lukas Kronenberg,
##          Llorenç Cabrera-Bosquet, Andreas Hund, Emilie J. Millet and Fred A. van Eeuwijk
##
## Contact: Diana M. Pérez-Valencia - dperez@bcamath.org
##..........................................................

## To navigate more easily through the code when using RStudio, we recommend to collapse it: Alt+o (Cmd+Alt+o on the Mac)

# Function to FIT the P-spline Hierarchical Curve Data Model --------------
# Input Arguments ---------------------------------------------------------
# response:     a character string with the name of the variable that contains the response variable of interest.
# time:         a character string with the name of the variable that contains the timepoints at which the response is measured for each plant/plot/individual.
# pop:          a character string with the name of the variable that contains the populations to which each genotype/variety belongs. This variable must be a factor in the data frame.
# geno:         a character string with the name of the variable that contains the genotypes/varieties to which each plant/plot/individual belongs. This variable must be a factor in the data frame.
# plant:        a character string with the name of the variable that contains the plants/plots/individuals. This variable must be a factor in the data frame.
# weights:      a character string with the name of the variable that contains the weights to be used in the fitting process. By default, the weights are considered to be one. The default is NULL.
# data:         a data frame containing all needed variables. The data must be complete. That is, all plants must be measured at all timepoints considered. If not, missing data should be included.
# dif.var:      a list that controls the variances for intercepts/slopes and smooth effect for genotypes (separately for each population) and plants (separately for each genotype).  The default is FALSE.
# smooth.pop:   a list specifying the P-Spline at the population level (nseg: number of segments; bdeg: degree of the B-spline basis; pord: penalty order).
# smooth.geno:  a list specifying the P-Spline model at the genotype/variety level.
# smooth.plant: a list specifying the P-Spline model at the plant/plot/individual belongs.
# offset:       an optional numerical vector containing an a priori known component to be included in the linear predictor during fitting. The default is NULL.
# family:       an object of class family specifying the distribution and link function.  The default is gaussian().
# maxit:        an optional value that controls the maximum number of iterations of the algorithm.  The default is 200.
# trace:        an optional value that controls the function trace.  The default is TRUE.
# thr:          an optional value that controls the convergence threshold of the algorithm. The default is 1.e-03.

# Output Arguments --------------------------------------------------------
# y:               a list with the original response for plants in each genotype.
# time:            a numeric vector with the timepoint.
# l.geno:          a vector with the names of the genotypes/varieties.
# l.pop:           a vector with the names of the populations.
# l.plant:         a vector with the names of the plants/plots/individuals.
# n.plants_p_pop:  a numeric vector with the number of plants/plots/individuals per population.
# n.geno_p_pop:    a numeric vector with the number of genotypes/varieties per population.
# n.plants_p_geno: a numeric vector with the number of plants/plots/individuals per genotype/variety.
# MM:              a list with the design matrices at plant, genotype and population levels.
# ed:              a numeric vector with the estimated effective dimension (or effective degrees of freedom) for each random component of the model (intercept, slope and non-linear trend) at each level of the hierarchy (population, genotype and plant)
# tot_ed:          a numeric value with the sum of the effective dimensions for all components of the model.
# vc:              a numeric vector with the (REML) variance component estimates for each random component of the model (intercept, slope and non-linear trend) at each level of the hierarchy (population, genotype and plant)
# phi:             a numeric value with the error variance estimate.
# coeff:           a numeric vector with the estimated fixed and random effect coefficients.
# eta_pop:         a list with the estimated population trajectories (eta_pop), and first (eta_pop_deriv1) and second (eta_pop_deriv2) order derivatives. Each element of the list is a numeric matrix where each column corresponds to one population.
# eta_geno:        a list with the estimated genotype/variety trajectories (eta_geno), and first (eta_geno_deriv1) and second (eta_geno_deriv2) order derivatives. Each element of the list is a list that corresponds to one population and is a matrix with as many columns as genotypes/varieties per population.
# eta_geno_dev:    a list with the estimated genotype/variety deviations (eta_geno_dev), and first (eta_geno_dev_deriv1) and second (eta_geno_dev_deriv2) order derivatives. Each element of the list is a list that corresponds to one population and is a matrix with as many columns as genotypes/varieties per population.
# eta_plant:       a list with the estimated plant/plot/individual trajectories (eta_plant), and first (eta_plant_deriv1) and second (eta_plant_deriv2) order derivatives. Each element of the list is a list that corresponds to one genotype/variety and is a matrix with as many columns as plants per genotype.
# eta_plant_dev:   a list with the estimated plant/plot/individual deviations (eta_plant_dev), and first (eta_plant_dev_deriv1) and second (eta_plant_dev_deriv2) order derivatives. Each element of the list is a list that corresponds to one genotype/variety and is a matrix with as many columns as plants per genotype.
# deviance:        the (REML) deviance at convergence.
# convergence:     a logical value indicating whether the algorithm managed to converge before the given number of iterations.
# dim:             a numeric vector with the (model) dimension of each model component (fixed and/or random) at each level of the hierarchy (population, genotype, and plant). These values correspond to the number of parameters to be estimated.
# family:          an object of class family specifying the distribution and link function.
# time_proc:       the function trace (in seconds)
# Vp:              the variance-covariance matrix for the coefficients.
# smooth:          a list with the information about number of segments (nseg), degree of the B-spline basis (bdeg) and penalty order (pord) for the three levels of the hierarchy.

# Function ----------------------------------------------------------------

#' fitSplineHDM
#'
#' fitSplineHDM
#'
#' @export
fitSplineHDM <- function(response,
                         time,
                         pop,
                         geno,
                         plant,
                         data,
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
  ## Unused levels might cause strange behaviour.
  data <- droplevels(data)
  ## Create a full data set of observations for all combinations of
  ## timepoints, row and column.
  timeDat <- data.frame(time = sort(unique(data[[time]])))
  colnames(timeDat) <- time
  fullGrid <- merge(unique(data[c(pop, geno, plant, "colId", "rowId")]),
                    timeDat)
  data <- merge(fullGrid, data, all.x = TRUE)

  # Order the data - no longer needed, merging in previous step does this.
  # data <- data[order(data[, pop], data[, geno], data[, plant], data[, time]), ]

  ## Normalize time.
  raw.time    <- timeDat[[time]]
  data[[time]] <- data[[time]] - min(data[[time]]) + 1
  ## Define offset.
  if (is.null(offset)) {
    data$offset <- 0
  }
  ## Specify weights.
  if(is.null(weights)) {
    weights <- rep(1, nrow(data))
  } else{
    weights <- data[[weights]] ^ 2
  }
  ## Convert to factors - only if not already factors.
  for (facVar in c(pop, geno, plant)) {
    if (!is.factor(data[[facVar]])) {
      data[[facVar]] <- factor(data[[facVar]])
    }
  }
  ## Elements and number of elements by level of the hierarchy.
  l.pop           <- levels(data[[pop]])
  l.geno          <- levels(data[[geno]])
  l.plant         <- levels(data[[plant]])
  n.pop           <- nlevels(data[[pop]])
  n.plants_p_pop  <- apply(X = table(data[[pop]], data[[plant]]),
                           MARGIN = 1, FUN = function(x) { sum(x!=0) })
  n.geno_p_pop    <- apply(X = table(data[[pop]], data[[geno]]),
                           MARGIN = 1, FUN = function(x) { sum(x!=0) })
  n.geno          <- nlevels(data[[geno]])
  n.plants_p_geno <- apply(table(data[[geno]], data[[plant]]),
                           MARGIN = 1, FUN = function(x) { sum(x!=0) })[l.geno]
  n.tot           <- nlevels(data[[plant]])
  # Time interval
  x <- sort(unique(data[[time]]))

  # Construct matrices: data in an array
  # Design matrices for the population curves
  MM.pop <- MM.basis(x, min(x), max(x), smooth.pop$nseg, smooth.pop$bdeg, smooth.pop$pord)
  X.pop  <- MM.pop$X
  Z.pop  <- MM.pop$Z

  # Design matrices for the genotype curves
  MM.geno <- MM.basis(x, min(x), max(x), smooth.geno$nseg, smooth.geno$bdeg, smooth.geno$pord)
  X.geno  <- MM.geno$X
  Z.geno <- MM.geno$Z

  # Design matrices for the individual curves
  MM.ind <- MM.basis(x, min(x), max(x), smooth.plant$nseg, smooth.plant$bdeg, smooth.plant$pord)
  X.ind  <- MM.ind$X
  Z.ind <- MM.ind$Z

  # Response, offset and weights
  y         <- data[,response]
  offset.f  <- data$offset
  weights.f <- weights

  # In case we have NAs
  nas            <- is.na(y)
  y[nas]         <- 0
  weights.f[nas] <- 0

  # Matrices assigning plants to populations and genotypes
  # Matrix to assign plants to populations
  ind.ind.pop <- rep(1:n.pop, n.plants_p_pop)
  if(n.pop == 1) {
    C.pop <- matrix(rep(1, n.tot), ncol = 1)
  } else {
    xxt <- data.frame(ind.pop = as.factor(ind.ind.pop))
    mft <- model.frame(~ ind.pop - 1, xxt, drop.unused.levels = TRUE)
    mtt <- terms(mft)
    f.termst <- attr(mtt, "term.labels")[attr(mtt,"dataClasses") == "factor"]
    C.pop <- model.matrix(mtt, data = mft, contrasts.arg = lapply(mft[,f.termst, drop = FALSE], contrasts, contrasts=FALSE))
    attr(mtt, "xlev") <- .getXlevels(mtt, mft)
    attr(mtt, "contrast") <- attr(C.pop,"contrast")
  }

  # Matrix to assign plants to genotypes
  ind.ind.geno <- rep(1:n.geno, n.plants_p_geno)
  if(n.geno == 1) {
    C.geno <- matrix(rep(1, n.tot), ncol = 1)
  } else {
    xxg <- data.frame(ind.pop = as.factor(ind.ind.geno))
    mfg <- model.frame(~ ind.pop - 1, xxg, drop.unused.levels = TRUE)
    mtg <- terms(mfg)
    f.termsg <- attr(mtg, "term.labels")[attr(mtg,"dataClasses") == "factor"]
    C.geno <- model.matrix(mtg, data = mfg, contrasts.arg = lapply(mfg[,f.termsg, drop = FALSE], contrasts, contrasts=FALSE))
    attr(mtg, "xlev") <- .getXlevels(mtg, mfg)
    attr(mtg, "contrast") <- attr(C.geno,"contrast")
  }

  # GLAM matrices
  GLAM = TRUE
  X <- list(X1 = Matrix::Matrix(C.pop), X2 = Matrix::Matrix(X.pop)) # Parametric population effect
  Z <- list(Z1 = list(Z11 = Matrix::Matrix(C.pop), Z12 = Matrix::Matrix(Z.pop)),
            Z2 = list(Z21 = Matrix::Matrix(C.geno), Z22 = Matrix::Matrix(X.geno)), 	 	# Random intercept and slope
            Z3 = list(Z31 = Matrix::Matrix(C.geno), Z32 = Matrix::Matrix(Z.geno)),	 	# Genotype specific curves
            Z4 = list(Z41 = Matrix::Matrix(diag(n.tot)), Z42 = Matrix::Matrix(X.ind)), 	# Random intercept and slopes (individual)
            Z5 = list(Z51 = Matrix::Matrix(diag(n.tot)), Z52 = Matrix::Matrix(Z.ind))) 	# Genotype specific curves (individual)

  # Number of parameters: fixed and random (for each component)
  np        <- c(ncol(X.pop)*n.pop, ncol(Z.pop)*n.pop, ncol(X.geno)*n.geno, ncol(Z.geno)*n.geno, ncol(X.ind)*n.tot, ncol(Z.ind)*n.tot)
  names(np) <- c("pop.fixed","pop.random","geno.x.random","geno.z.random","ind.x.random", "ind.z.random")
  np.comp   <- c(np[1]+np[2], np[3]+np[4], np[5]+np[6])
  names(np.comp) <- c("pop", "geno", "plant")
  np.e      <- cumsum(np.comp)
  np.s      <- np.e - np.comp + 1

  # Construct precision matrix
  # Smooth main effect: one per population
  g <- lapply(1:n.pop, function(x) MM.pop$d)

  if(isTRUE(dif.var$geno)){
    # Random intercepts and slopes (genotype)
    for (i in 1:n.pop) {
      g[[n.pop+i]] <- lapply(1:smooth.geno$pord, function(j) rep(sapply(1:smooth.geno$pord, function(k) ifelse(k==j,1,0)), n.geno_p_pop[i]))
    }
    # Smooth effects (genotype)
    for (i in 1:n.pop) {
      g[[n.pop*2+i]] <- rep(MM.geno$d, n.geno_p_pop[i])
    }
    if (isTRUE(dif.var$plant)){
      # Random intercepts and slopes (individual)
      for (i in 1:n.geno) {
        g[[n.pop*3+i]] <- lapply(1:smooth.plant$pord, function(j) rep(sapply(1:smooth.plant$pord, function(k) ifelse(k==j,1,0)), n.plants_p_geno[i]))
      }
      # Smooth effects (individual)
      for (i in 1:n.geno) {
        g[[n.pop*3 + n.geno + i]] <- rep(MM.ind$d, n.plants_p_geno[i])
      }
    } else{
      # Random intercepts and slopes (individual)
      g[[n.pop*3 + 1]] <- lapply(1:smooth.plant$pord, function(j) rep(sapply(1:smooth.plant$pord, function(k) ifelse(k==j,1,0)), n.tot))
      # Smooth effects (individual)
      g[[n.pop*3 + 2]] <- list(rep(MM.ind$d, n.tot))
    }
  } else if (isFALSE(dif.var$geno)){
    # Random intercepts and slopes (genotype)
    g[[n.pop + 1]] <- lapply(1:smooth.geno$pord, function(j) rep(sapply(1:smooth.geno$pord, function(k) ifelse(k==j,1,0)), n.geno))
    # Smooth effects (genotype)
    g[[n.pop + 2]] <- list(rep(MM.geno$d, n.geno))
    if(isTRUE(dif.var$plant)){
      # Random intercepts and slopes (individual)
      for (i in 1:n.geno) {
        g[[n.pop + 2 + i]] <- lapply(1:smooth.plant$pord, function(j) rep(sapply(1:smooth.plant$pord, function(k) ifelse(k==j,1,0)), n.plants_p_geno[i]))
      }
      # Smooth effects (individual)
      for (i in 1:n.geno) {
        g[[n.pop + 2 + n.geno + i]] <- rep(MM.ind$d, n.plants_p_geno[i])
      }
    } else{
      # Random intercepts and slopes (individual)
      g[[n.pop + 3]] <- lapply(1:smooth.plant$pord, function(j) rep(sapply(1:smooth.plant$pord, function(k) ifelse(k==j,1,0)), n.tot))
      # Smooth effects (individual)
      g[[n.pop + 4]] <- list(rep(MM.ind$d, n.tot))
    }
  }

  # Construct the components of the precision matrix (as needed by the algorithm)
  g <- construct.capital.lambda(g)

  # Initialise the parameters
  la = c(1,rep(1, l = length(g)))
  devold = 1e10

  mustart <- etastart <- NULL
  nobs    <- length(y)
  eval(family$initialize)
  mu      <- mustart
  eta     <- family$linkfun(mustart)

  # Iteration process to estimate coefficients and variance components
  for (iter in 1:maxit) {
    deriv <- family$mu.eta(eta)
    z     <- (eta - offset.f) + (y - mu)/deriv
    w     <- as.vector(deriv^2/family$variance(mu))
    w     <- w*weights.f
    mat   <- construct.matrices(X, Z, z, w, GLAM)
    V     <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., mat$ZtZ.) # Does not change from iteration to iteration
    # V is an object of class Matrix, transform it to spam
    V     <- spam::as.spam.dgCMatrix(V)

    # First iteration (cholesky decomposition)
    Ginv  <- vector(length = length(g[[1]]))
    for (i in 1:length(g)) {
      Ginv <- Ginv + g[[i]]
    }
    G      <- 1/Ginv
    D      <- spam::bdiag.spam(spam::diag.spam(rep(0L, np[1]), nrow = np[1], ncol = np[1]), spam::diag.spam(Ginv))
    cholHn <- chol(V + D) # The sparsness structure does not change (Ginv is diagonal)
    for (it in 1:maxit) {
      Ginv <- vector(length = length(g[[1]]))
      for (i in 1:length(g)) {
        Ginv <- Ginv + (1/la[i+1])*g[[i]]
      }
      G  <- 1/Ginv
      D  <- spam::bdiag.spam(spam::diag.spam(rep(0L, np[1]), nrow = np[1], ncol = np[1]), spam::diag.spam(Ginv))
      Hn <- (1/la[1])*V + D
      cholHn     <- update(cholHn, Hn)
      b          <- spam::backsolve(cholHn, spam::forwardsolve(cholHn, (1/la[1])*mat$u))
      cholHn.inv <- spam::forwardsolve.spam(cholHn, spam::diag.spam(rep(1L, nrow(cholHn))))
      Hninv.diag <- colSums(cholHn.inv*cholHn.inv)
      aux        <- G - Hninv.diag[-(1:np[1])]

      # Fixed and random coefficients
      b.fixed  <- b[1:np[1]]
      b.random <- b[-(1:np[1])]

      # Variance components as in SOP
      ssv <- ed <- tau <- vector(mode="list", length=length(g))
      for (i in 1:length(g)) {
        ed[[i]]  <- (1/la[i+1])*sum(aux*g[[i]])
        ed[[i]]  <- ifelse(ed[[i]] <= 1e-10, 1e-10, ed[[i]])
        ssv[[i]] <- sum(b.random^2*g[[i]])
        tau[[i]] <- ssv[[i]]/ed[[i]]
        tau[[i]] <- ifelse(tau[[i]] <= 1e-10, 1e-10, tau[[i]])
      }
      ssr <- mat$yty. - t(c(b.fixed, b.random)) %*% (2*mat$u - V %*% b)
      dev <- deviance_spam(cholHn, spam::diag.spam(G), w[w != 0], la[1], ssr, sum(b.random^2*Ginv))[1]

      if(family$family == "gaussian" | family$family == "quasipoisson") {
        phi <- as.numeric((ssr/(length(z[w != 0]) - sum(unlist(ed)) - np[1])))
      } else {
        phi <- 1
      }
      # New variance components and convergence check
      lanew <- c(phi, unlist(tau))
      dla   <- abs(devold - dev)
      if(trace) {
        cat(sprintf("%1$3d %2$12.6f", it, dev), sep = "")
        cat(sprintf("%12.3f", unlist(ed)), sep = "")
        cat('\n')
      }
      if (dla < thr) break
      la     <- lanew
      devold <- dev
    }
    eta.old <- eta
    eta     <- Xtheta(X, b.fixed) + Ztheta(Z, b.random, np[-1]) + offset.f
    mu      <- family$linkinv(eta)

    # Convergence criterion: linear predictor
    tol <- sum((eta - eta.old)^2)/sum(eta^2)
    if (tol < 1e-6 | family$family == "gaussian") break

  }

  coeff <- c(b.fixed, b.random)

  # # Create the objects to be returned
  #
  # # Growth curves and deviations at population, genotype and plant level
  #
  # # Calculate the estimated population functions
  # # First the model matrices
  # if(n.pop == 1) {
  #   mm.ind.pop <- matrix(1, ncol = 1, nrow = 1)
  # } else {
  #   mfpt       <- model.frame(mtt, data.frame(ind.pop = as.factor(1:n.pop)), xlev = attr(mtt, "xlev"))
  #   mm.ind.pop <- model.matrix(mtt, data = mfpt, contrasts.arg = attr(mtt, "contrast"))
  # }
  # dm.pop <- cbind(mm.ind.pop%x%X.pop, mm.ind.pop%x%Z.pop)
  #
  # # The first coefficients are for the populations. Matrix with as many columns as populations
  # # Population-specific growth curves
  # eta_pop           <- matrix(dm.pop%*%coeff[1:ncol(dm.pop)], ncol = n.pop)
  # colnames(eta_pop) <- l.pop
  #
  # # Calculate the estimated genotype functions
  # # First the model matrices
  # if(n.geno == 1) {
  #   mm.ind.geno <- matrix(1, ncol = 1, nrow = 1)
  # } else {
  #   mfpg        <- model.frame(mtg, data.frame(ind.pop = as.factor(1:n.geno)), xlev = attr(mtg, "xlev"))
  #   mm.ind.geno <- model.matrix(mtg, data = mfpg, contrasts.arg = attr(mtg, "contrast"))
  # }
  #
  # # Genotype-specific deviations
  # ind.geno.pop <- rep(1:n.pop, n.geno_p_pop)
  # dm.geno      <- cbind(kronecker(Matrix::Matrix(mm.ind.geno), X.geno), kronecker(Matrix::Matrix(mm.ind.geno), Z.geno))
  # aux          <- matrix(dm.geno%*%coeff[(sum(np[1:2]) + 1):(sum(np[1:2]) + ncol(dm.geno))], ncol = n.geno)
  # # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per population.
  # eta_geno_dev <- lapply(split(aux, rep(ind.geno.pop, each = nrow(aux))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(eta_geno_dev) <- l.pop
  #
  # # Genotype-specific growth curves
  # if(n.pop == 1) {
  #   mm.ind.pop <- matrix(1, ncol = 1, nrow = n.geno)
  # } else {
  #   mfpt       <- model.frame(mtt, data.frame(ind.pop = rep(as.factor(1:n.pop), n.geno_p_pop)), xlev = attr(mtt, "xlev"))
  #   mm.ind.pop <- model.matrix(mtt, data = mfpt, contrasts.arg = attr(mtt, "contrast"))
  # }
  # dm.pop.geno     <- cbind(mm.ind.pop%x%X.pop, mm.ind.pop%x%Z.pop)
  # dm.geno         <- cbind(dm.pop.geno, kronecker(Matrix::Matrix(mm.ind.geno), Matrix::Matrix(X.geno)), kronecker(Matrix::Matrix(mm.ind.geno), Matrix::Matrix(Z.geno)))
  # aux             <- matrix(dm.geno%*%coeff[1:ncol(dm.geno)], ncol = n.geno)
  # # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per threatment
  # eta_geno        <- lapply(split(aux, rep(ind.geno.pop, each = nrow(aux))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(eta_geno) <- l.pop
  #
  # # Give the name of the genotypes
  # e <- cumsum(n.geno_p_pop)
  # s <- e - n.geno_p_pop + 1
  # for(i in 1:n.pop) {
  #   colnames(eta_geno[[i]])     <- l.geno[s[i]:e[i]]
  #   colnames(eta_geno_dev[[i]]) <- l.geno[s[i]:e[i]]
  # }
  #
  # # Calculate the estimated plant functions
  # # Agregated at the genotype level (list, with as many elements as genotypes. Each element of the list is a matrix, with as many columns as plants per genotype)
  # # Plant-specific deviations
  # aux            <- matrix(eta, ncol = n.tot)
  # eta_ind        <- lapply(split(aux, rep(ind.ind.geno, each = nrow(aux))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(eta_ind) <- l.geno
  #
  # # Plant raw data
  aux            <- matrix(data[,response], ncol = n.tot)
  obs_ind        <- lapply(split(aux, rep(ind.ind.geno, each = nrow(aux))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  names(obs_ind) <- l.geno
  #
  # # Plant-specific growth curves
  # MM           <- cbind(kronecker(Matrix::Matrix(diag(n.tot)), X.ind), kronecker(Matrix::Matrix(diag(n.tot)), Z.ind))
  # aux          <- matrix(MM%*%coeff[(sum(np[1:4]) + 1):(sum(np[1:4]) + ncol(MM))], ncol = n.tot)
  # eta_ind_dev  <- lapply(split(aux, rep(ind.ind.geno, each = nrow(aux))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(eta_ind_dev) <- l.geno
  #
  # # Obtain 1st and 2nd order derivatives: at population, genotype and plant level
  #
  # # Needed:
  # # knots 	  knots used in the fit
  # # x 		    covariate values where to obtain the derivative
  # # bdeg      order of the B-spline basis (used in the fit)
  # # deriv     desired derivative
  # # theta     estimated coefficients for the fit (for the B-spline model, not for the reformulation). f(x) = B*theta -> theta = U.X*b.fixed + U.Z*b.random
  #
  # # Functions at popopulation level
  # # Transformation matrix, theta and B
  # if(n.pop == 1) {
  #   mm.ind.pop <- matrix(1, ncol = 1, nrow = 1)
  # } else {
  #   mfpt       <- model.frame(mtt, data.frame(ind.pop = as.factor(1:n.pop)), xlev = attr(mtt, "xlev"))
  #   mm.ind.pop <- model.matrix(mtt, data = mfpt, contrasts.arg = attr(mtt, "contrast"))
  # }
  # T_pop     <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), MM.pop$U.X), kronecker(Matrix::Matrix(mm.ind.pop), MM.pop$U.Z))
  # theta_pop <- T_pop%*%coeff[np.s[1]:np.e[1]]
  #
  # # First derivative
  # B_pop_deriv1 <- kronecker(Matrix::Matrix(diag(n.pop)), spline.bbase(knots = MM.pop$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 1))
  # f_pop_deriv1 <- matrix(B_pop_deriv1%*%theta_pop, ncol = n.pop)
  # colnames(f_pop_deriv1) <- l.pop
  #
  # # Second derivative
  # B_pop_deriv2 <- kronecker(Matrix::Matrix(diag(n.pop)), spline.bbase(knots = MM.pop$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 2))
  # f_pop_deriv2 <- matrix(B_pop_deriv2%*%theta_pop, ncol = n.pop)
  # colnames(f_pop_deriv2) <- l.pop
  #
  # # Functions at genotype level
  # # Genotype-specific deviations
  # if(n.geno == 1) {
  #   mm.ind.geno <- matrix(1, ncol = 1, nrow = 1)
  # } else {
  #   mfpg        <- model.frame(mtg, data.frame(ind.pop = as.factor(1:n.geno)), xlev = attr(mtg, "xlev"))
  #   mm.ind.geno <- model.matrix(mtg, data = mfpg, contrasts.arg = attr(mtg, "contrast"))
  # }
  # T_geno_dev     <- cbind(kronecker(Matrix::Matrix(diag(n.geno)), MM.geno$U.X), kronecker(Matrix::Matrix(diag(n.geno)), MM.geno$U.Z))
  # theta_geno_dev <- T_geno_dev%*%coeff[np.s[2]:np.e[2]]
  #
  # # First derivative
  # B_geno_dev_deriv1 <- kronecker(Matrix::Matrix(diag(n.geno)), spline.bbase(knots = MM.geno$knots, X. = x, BDEG. = smooth.geno$bdeg, deriv = 1))
  # f_geno_dev_deriv1 <- matrix(B_geno_dev_deriv1%*%theta_geno_dev, ncol = n.geno)
  # f_geno_dev_deriv1 <- lapply(split(f_geno_dev_deriv1, rep(ind.geno.pop, each = nrow(f_geno_dev_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_geno_dev_deriv1) <- l.pop
  #
  # # Second derivative
  # B_geno_dev_deriv2 <- kronecker(Matrix::Matrix(diag(n.geno)), spline.bbase(knots = MM.geno$knots, X. = x, BDEG. = smooth.geno$bdeg, deriv = 2))
  # f_geno_dev_deriv2 <- matrix(B_geno_dev_deriv2%*%theta_geno_dev, ncol = n.geno)
  # f_geno_dev_deriv2 <- lapply(split(f_geno_dev_deriv2, rep(ind.geno.pop, each = nrow(f_geno_dev_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_geno_dev_deriv2) <- l.pop
  #
  # # Genotype-specific growth curves
  # if(n.pop == 1) {
  #   mm.ind.pop <- matrix(1, ncol = 1, nrow = n.geno)
  # } else {
  #   ind.pop    <- rep(as.factor(1:n.pop), n.geno_p_pop)
  #   mm.ind.pop <- model.matrix(~ 0 + ind.pop) # The contrast matrix changes here!!!!!!
  # }
  # T_geno     <- Matrix::bdiag(T_pop, T_geno_dev)
  # theta_geno <- T_geno%*%coeff[np.s[1]:np.e[2]]
  #
  # # First derivative
  # B_geno_deriv1 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = MM.pop$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 1)), B_geno_dev_deriv1)
  # f_geno_deriv1 <- matrix(B_geno_deriv1%*%theta_geno, ncol = n.geno)
  # f_geno_deriv1 <- lapply(split(f_geno_deriv1, rep(ind.geno.pop, each = nrow(f_geno_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_geno_deriv1) <- l.pop
  #
  # # Second derivative
  # B_geno_deriv2 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = MM.pop$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 2)), B_geno_dev_deriv2)
  # f_geno_deriv2 <- matrix(B_geno_deriv2%*%theta_geno, ncol = n.geno)
  # f_geno_deriv2 <- lapply(split(f_geno_deriv2, rep(ind.geno.pop, each = nrow(f_geno_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_geno_deriv2) <- l.pop
  #
  # # Give the name of the genotypes
  # e <- cumsum(n.geno_p_pop)
  # s <- e - n.geno_p_pop + 1
  # for(i in 1:n.pop) {
  #   colnames(f_geno_dev_deriv1[[i]]) <- l.geno[s[i]:e[i]]
  #   colnames(f_geno_dev_deriv2[[i]]) <- l.geno[s[i]:e[i]]
  #   colnames(f_geno_deriv1[[i]])     <- l.geno[s[i]:e[i]]
  #   colnames(f_geno_deriv2[[i]])     <- l.geno[s[i]:e[i]]
  # }
  #
  # # Functions at plant level
  # # Plant-specific deviations
  # # Transformation matrix, theta and B
  # T_plant_dev     <- cbind(kronecker(Matrix::Matrix(diag(n.tot)), Matrix::Matrix(MM.ind$U.X)), kronecker(Matrix::Matrix(diag(n.tot)), Matrix::Matrix(MM.ind$U.Z)))
  # theta_plant_dev <- Matrix::Matrix(T_plant_dev%*%coeff[np.s[3]:np.e[3]])
  #
  # # First derivative
  # B_plant_dev_deriv1 <- kronecker(Matrix::Matrix(diag(n.tot)), spline.bbase(knots = MM.ind$knots, X. = x, BDEG. = smooth.plant$bdeg, deriv = 1))
  # f_plant_dev_deriv1 <- matrix(B_plant_dev_deriv1%*%theta_plant_dev, ncol = n.tot)
  # f_plant_dev_deriv1 <- lapply(split(f_plant_dev_deriv1, rep(ind.ind.geno, each = nrow(f_plant_dev_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_plant_dev_deriv1) <- l.geno
  #
  # # Second derivative
  # B_plant_dev_deriv2 <- kronecker(Matrix::Matrix(diag(n.tot)), spline.bbase(knots = MM.ind$knots, X. = x, BDEG. = smooth.plant$bdeg, deriv = 2))
  # f_plant_dev_deriv2 <- matrix(B_plant_dev_deriv2%*%theta_plant_dev, ncol = n.tot)
  # f_plant_dev_deriv2 <- lapply(split(f_plant_dev_deriv2, rep(ind.ind.geno, each = nrow(f_plant_dev_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_plant_dev_deriv2) <- l.geno
  #
  # # Plant-specific growth curves
  # if(n.pop == 1) {
  #   mm.ind.pop <- matrix(1, ncol = 1, nrow = n.plants_p_pop)
  # } else {
  #   ind.pop    <- rep(as.factor(1:n.pop), n.plants_p_pop)
  #   mm.ind.pop <- model.matrix(~ 0 + ind.pop) # The contrast matrix changes here!!!!!!
  # }
  # if(n.geno == 1) {
  #   mm.ind.geno <- matrix(1, ncol = 1, nrow = n.geno)
  # } else {
  #   ind.geno    <- rep(as.factor(1:n.geno), n.plants_p_geno)
  #   mm.ind.geno <- model.matrix(~ 0 + ind.geno) # The contrast matrix changes here!!!!!!
  # }
  # T_plant     <- Matrix::bdiag(T_pop, T_geno_dev, T_plant_dev)
  # theta_plant <- T_plant%*%coeff[np.s[1]:np.e[3]]
  #
  # # First derivative
  # B_plant_deriv1  <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = MM.pop$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 1)),
  #                          kronecker(Matrix::Matrix(mm.ind.geno), spline.bbase(knots = MM.geno$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 1)),
  #                          B_plant_dev_deriv1)
  # f_plant_deriv1  <- matrix(B_plant_deriv1%*%theta_plant, ncol = n.tot)
  # f_plant_deriv1  <- lapply(split(f_plant_deriv1, rep(ind.ind.geno, each = nrow(f_plant_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_plant_deriv1) <- l.geno
  #
  # # Second derivative
  # B_plant_deriv2 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = MM.pop$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 2)),
  #                         kronecker(Matrix::Matrix(mm.ind.geno), spline.bbase(knots = MM.geno$knots, X. = x, BDEG. = smooth.pop$bdeg, deriv = 2)),
  #                         B_plant_dev_deriv2)
  # f_plant_deriv2 <- matrix(B_plant_deriv2%*%theta_plant, ncol = n.tot)
  # f_plant_deriv2 <- lapply(split(f_plant_deriv2, rep(ind.ind.geno, each = nrow(f_plant_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(x))
  # names(f_plant_deriv2) <- l.geno

  # Object to be returned
  res                 <- list()
  res$y               <- obs_ind
  res$time            <- raw.time
  res$l.geno          <- l.geno
  res$l.pop           <- l.pop
  res$l.plant         <- l.plant
  res$n.plants_p_pop  <- n.plants_p_pop
  res$n.geno_p_pop    <- n.geno_p_pop
  res$n.plants_p_geno <- n.plants_p_geno
  res$MM              <- list(MM.pop = MM.pop, MM.geno = MM.geno, MM.ind = MM.ind)
  res$ed              <- unlist(ed)
  res$tot_ed          <- sum(unlist(ed))
  res$vc              <- unlist(tau)
  res$phi             <- phi
  res$coeff           <- coeff
  # res$eta_pop         <- list(eta_pop = eta_pop, eta_pop_deriv1 = f_pop_deriv1, eta_pop_deriv2 = f_pop_deriv2)
  # res$eta_geno        <- list(eta_geno = eta_geno, eta_geno_deriv1 = f_geno_deriv1, eta_geno_deriv2 = f_geno_deriv2)
  # res$eta_geno_dev    <- list(eta_geno_dev = eta_geno_dev, eta_geno_dev_deriv1 = f_geno_dev_deriv1, eta_geno_dev_deriv2 = f_geno_dev_deriv2)
  # res$eta_plant       <- list(eta_plant = eta_ind, eta_plant_deriv1 = f_plant_deriv1, eta_plant_deriv2 = f_plant_deriv2)
  # res$eta_plant_dev   <- list(eta_plant_dev = eta_ind_dev, eta_plant_dev_deriv1 = f_plant_dev_deriv1, eta_plant_dev_deriv2 = f_plant_dev_deriv2)
  res$deviance        <- dev
  res$convergence     <- ifelse(it < maxit, TRUE, FALSE)
  res$dim             <- np
  res$family          <- family
  res$Vp              <- spam::chol2inv(cholHn) #Hinv
  res$smooth          <- list(smooth.pop = smooth.pop, smooth.geno = smooth.geno, smooth.plant = smooth.plant)
  class(res)          <- "psHDM"
  return(res)
}

spline.bbase <- function(knots, X., BDEG., deriv = 0, eps = 1e-15) {
  B <- splines::spline.des(knots = knots, x = X., derivs = deriv, ord = BDEG. + 1, outer.ok = TRUE)$design
  B
}

bbase <- function(x,
                  xl = min(x),
                  xr = max(x),
                  ndx = 10,
                  bdeg = 3,
                  eps = 1e-15) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- splines::spline.des(knots, x, ord = bdeg + 1, outer.ok = TRUE)$design
  res <- list(B = B, knots = knots)
  res
}

MM.basis <- function (x, xl, xr, ndx, bdeg, pord) {
  Bb <- bbase(x,xl,xr,ndx,bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B)
  n <- nrow(B)
  D <- diff(diag(m), differences=pord)
  P.svd <- svd(crossprod(D))
  d <- (P.svd$d)[1:(m-pord)]
  U.Z <- (P.svd$u)[,1:(m-pord)]
  Z <- B%*%U.Z
  X <- NULL
  for(i in 0:(pord-1)){
    X <- cbind(X,x^i)
  }
  U.X <- NULL
  for(i in 0:(pord-1)){
    U.X <- cbind(U.X, knots[-c((1:(bdeg - 1)),(length(knots)- (bdeg - 1) + 1):length(knots))]^i)
  }
  list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
}

construct.capital.lambda <- function(g) {
  length.eq <- all(sapply(g, function(x) {
    diff(range(unlist(lapply(x, length)))) < .Machine$double.eps ^ 0.5
  }))
  if(length.eq) {
    l <- length(g)
    if(l == 1) {
      if(length(g[[1]]) == 1) {
        res <- g
      } else {
        res <- do.call("c", lapply(g, function(x) x))
      }
    } else {
      dim <- sapply(g, function(x) {
        if(is.list(x))
          unlist(lapply(x, length))[1]
        else
          length(x)
      })
      end <- cumsum(dim)
      init <- end - dim + 1

      res <- do.call("c", lapply(1:length(g), function(x, g, init, end, dim) {
        temp <- g[[x]]
        if(is.list(temp)) {
          lapply(temp, function(y, x, dim) {
            aux <- rep(0L, l = sum(dim))
            aux[init[x]:end[x]] <- y
            aux
          }, x = x, dim = dim)
        } else {
          aux <- rep(0L, l = sum(dim))
          aux[init[x]:end[x]] <- temp
          list(aux)
        }
      }, g = g, init = init, end = end, dim = dim))
    }
  } else {
    stop("Error in construct.capital.lambda")
  }
  res
}

construct.block <- function(A1,A2,A3,A4) {
  block <- rbind(cbind(A1,A2), cbind(A3,A4))
  return(block)
}

Rten <- function(X) {
  one <- Matrix::Matrix(1, 1, ncol(X))
  kronecker(X,one)*kronecker(one,X)
}

Rten2 <- function(X1,X2) {
  one.1 <- Matrix::Matrix(1,1,ncol(X1))
  one.2 <- Matrix::Matrix(1,1,ncol(X2))
  kronecker(X1,one.2)*kronecker(one.1,X2)
}

H <- function(X,A) {
  d <- dim(A)
  M <- Matrix::Matrix(A, nrow = d[1])
  XM <- X%*%M
  slam::as.simple_sparse_array(array(XM, c(nrow(XM),d[-1])))
}

Rotate <- function(A) {
  d <- 1:length(dim(A))
  d1 <- c(d[-1],d[1])
  slam::as.simple_sparse_array(aperm(A, d1))
}

RH <- function(X,A) {
  Rotate(H(X,A))
}

A1.form <- function(l, w = NULL){
  d <- length(l)
  n <- rev(sapply(l,nrow))
  c <- rev(sapply(l,ncol))
  if (is.null(w)) {
    W <- slam::as.simple_sparse_array(array(1, n))
  } else {
    W <- slam::as.simple_sparse_array(array(w, n))
  }
  tmp <- RH(Matrix::t(Rten(l[[d]])), W)
  for (i in (d-1):1) {
    tmp <- RH(Matrix::t(Rten(l[[i]])),tmp)
  }
  dim(tmp)<- rep(c, rep(2,d))
  Fast1 <- slam::as.simple_sparse_array(aperm(tmp, as.vector(matrix(1:(d*2), byrow = TRUE, ncol = 2))))
  Fast <- if(prod(c)) Matrix::Matrix(Fast1, nrow = prod(c)) else Fast1
  return(Fast)
}

A2.form <- function(l1, l2, w = NULL) {
  d1 <- length(l1)
  d2 <- length(l2)
  if(!(d1 == d2)) {
    stop("l1 and l2 should have the same dimension")
  }
  n <- rev(sapply(l1, nrow))
  d <- rev(sapply(l1, ncol))
  c <- rev(sapply(l2, ncol))

  if (is.null(w)) {
    W <- slam::as.simple_sparse_array(array(1, n))
  } else {
    W <- slam::as.simple_sparse_array(array(w, n))
  }
  tmp <- RH(Matrix::t(Rten2(l2[[d1]], l1[[d1]])), W)
  for (i in (d1-1):1) {
    tmp <- RH(Matrix::t(Rten2(l2[[i]], l1[[i]])),tmp)
  }
  dim(tmp)<- as.vector(rbind(d,c))
  Fast1 <- slam::as.simple_sparse_array(aperm(tmp, as.vector(matrix(1:(d1*2), byrow = TRUE, ncol = 2))))
  Fast <- if(prod(d)) Matrix::Matrix(Fast1, nrow = prod(d)) else Fast1
  return(Fast)
}

XtX <- function(X, w = NULL) {
  A1.form(X, w)
}

XtZ <- function(X, Z, w = NULL) {
  d <- length(Z)
  res <- NULL
  for (i in 1:d) {
    res <- cbind(res, A2.form(X, Z[[i]], w))
  }
  res
}

ZtZ <- function(Z, w = NULL) {
  d <- length(Z)
  upper <- list()
  for(i in 1:d) {
    upper[[i]] <- list()
    upper[[i]][[i]] <- A1.form(Z[[i]], w)
  }
  # Obtain the elements of the matrix
  if(d > 1)
    for(i in 1:(d-1)) {
      for(j in (i+1):d) {
        upper[[i]][[j]] <- A2.form(Z[[i]], Z[[j]], w)
      }
    }
  # Create the matrix
  res <- NULL
  for (i in 1:d) {
    if( i == 1) {
      res <- do.call("cbind", upper[[1]])
    } else {
      tmp <- do.call("cbind", upper[[i]])
      for(j in (i-1):1) {
        if(length(upper[[j]][[i]]))
          tmp <- cbind(Matrix::t(upper[[j]][[i]]), tmp)
      }
      if(nrow(tmp))
        res <- rbind(res, tmp)
    }
  }
  res
}

Xty <- function(X, y, w = NULL) {
  d <- length(X)
  n <- rev(sapply(X, nrow))
  if(is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w*y, n)
  }
  tmp <- RH(Matrix::t(X[[d]]),Y)
  for(i in (d-1):1) {
    tmp <- RH(Matrix::t(X[[i]]), tmp)
  }
  as.vector(tmp)
}

Zty <- function(Z, y, w = NULL) {
  d <- length(Z)
  n <- rev(sapply(Z[[1]], nrow))
  if(is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w*y, n)
  }
  res <- NULL
  for(i in 1:d) {
    k <- length(Z[[i]])
    tmp <- RH(Matrix::t(Z[[i]][[k]]),Y)
    for(j in (k-1):1) {
      tmp <- RH(Matrix::t(Z[[i]][[j]]), tmp)
    }
    res <- c(res, as.vector(tmp))
  }
  res
}

Xtheta <- function(X, theta) {
  d <- length(X)
  n <- rev(sapply(X, ncol))
  Theta <- array(theta, n)
  tmp <- RH(X[[d]], Theta)
  for(i in (d-1):1) {
    tmp <- RH(X[[i]], tmp)
  }
  as.vector(tmp)
}

Ztheta <- function(Z, theta, np) {
  d <- length(Z)
  for(i in 1:d) {
    if (i == 1) {
      res <- Xtheta(Z[[i]], theta[1:(np[1])])
    } else {
      init <- sum(np[1:(i-1)])
      fin  <- np[i]
      if(fin) res <- res + Xtheta(Z[[i]], theta[(init+1):(init+fin)])
    }
  }
  res
}

construct.matrices <- function(X, Z, z, w, GLAM) {
  if(GLAM) {
    XtX. <- XtX(X,w)
    XtZ. <- XtZ(X,Z,w)
    ZtX. <- Matrix::t(XtZ.)
    ZtZ. <- ZtZ(Z,w)
    Zty. <- Zty(Z,z,w)
    Xty. <- Xty(X,z,w)
    yty. <- sum((z^2)*w)
    ZtXtZ <- rbind(XtZ., ZtZ.)
    u <- c(Xty.,Zty.)
  } else {
    XtW. = t(X*w)
    XtX. = XtW.%*%X
    XtZ. = XtW.%*%Z
    ZtX. = Matrix::t(XtZ.)
    ZtW. = Matrix::t(Z*w)
    ZtZ. = ZtW.%*%Z
    Xty. = XtW.%*%z
    Zty. = ZtW.%*%z
    yty. <- sum((z^2)*w)
    ZtXtZ = rbind(XtZ., ZtZ.)
    u <- c(Xty.,Zty.)
  }
  res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
}

construct.G.matrix <- function(g, la) {
  la.g <- la[-1]
  g.comp <- length(g)
  g.tot.comp <- unlist(lapply(g, length))
  np.e.g <- cumsum(g.tot.comp)
  np.s.g <- np.e.g - g.tot.comp + 1

  Gl <- list(length = g.comp)
  for(i in 1:g.comp) {
    Gl[[i]] <- Reduce('+', mapply("*", g[[i]], la.g[np.s.g[i]:np.e.g[i]]))
  }
  G <- spam::bdiag.spam(Matrix:::bdiag(Gl))
  G
}

construct.Ginv.matrix <- function(g, la) {
  la.g <- la[-1]
  g.comp <- length(g)
  g.tot.comp <- unlist(lapply(g, length))
  np.e.g <- cumsum(g.tot.comp)
  np.s.g <- np.e.g - g.tot.comp + 1

  Gl <- list(length = g.comp)
  for(i in 1:g.comp) {
    Gl[[i]] <- solve(Reduce('+', mapply("*", g[[i]], la.g[np.s.g[i]:np.e.g[i]])))
  }
  Ginv <- spam::bdiag.spam(Matrix:::bdiag(Gl))
  Ginv
}

deviance_spam <- function(C, G, w, sigma2, ssr, edf) {
  log_det_C <- spam::determinant(C)$modulus*2
  log_det_G <- spam::determinant(G)$modulus
  deviance <- log_det_C + log_det_G + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
  deviance
}









