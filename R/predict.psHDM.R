# Function to PREDICT the P-spline Hierarchical Curve Data Model ----------
# Input Arguments ---------------------------------------------------------
# object:     a model object of class "psHDM" as obtained using function psHDM().
# newtimes:   a numeric vector with timepints at which predictions are required
# pred:       a list that controls the hierarchical levels at which predictions are required (Population, genotypes and plants).  The default is TRUE.
# se:         a list that controls the hierarchical levels at which standard errors are required (Population, genotypes and plants).  The default is TRUE except at the plant level.

# Output Arguments --------------------------------------------------------
# newtimes:        a numeric vector with the timepoint at which predictions and standard errors have been obtained.
# raw.time:        a numeric vector with the original timepoints.
# y:               a list with the original response for plants in each genotype.
# l.geno:          a vector with the names of the genotypes/varieties.
# l.pop:           a vector with the names of the populations.
# l.plant:         a vector with the names of the plants/plots/individuals.
# n.plants_p_pop:  a numeric vector with the number of plants/plots/individuals per population.
# n.geno_p_pop:    a numeric vector with the number of genotypes/varieties per population.
# n.plants_p_geno: a numeric vector with the number of plants/plots/individuals per genotype/variety.
# f_pop:           a list with the estimated population trajectories (fp_pop), and first (fp_pop_deriv1) and second (fp_pop_deriv2) order derivatives. Each element of the list is a numeric matrix where each column corresponds to one population.
# f_geno:          a list with the estimated genotype/variety trajectories (fp_geno), and first (fp_geno_deriv1) and second (fp_geno_deriv2) order derivatives. Each element of the list is a list that corresponds to one population and is a matrix with as many columns as genotypes/varieties per population.
# f_geno_dev:      a list with the estimated genotype/variety deviations (fp_geno_dev), and first (fp_geno_dev_deriv1) and second (fp_geno_dev_deriv2) order derivatives. Each element of the list is a list that corresponds to one population and is a matrix with as many columns as genotypes/varieties per population.
# f_plant:         a list with the estimated plant/plot/individual trajectories (fp_plant), and first (fp_plant_deriv1) and second (fp_plant_deriv2) order derivatives. Each element of the list is a list that corresponds to one genotype/variety and is a matrix with as many columns as plants per genotype.
# f_plant_dev:     a list with the estimated plant/plot/individual deviations (fp_plant_dev), and first (fp_plant_dev_deriv1) and second (fp_plant_dev_deriv2) order derivatives. Each element of the list is a list that corresponds to one genotype/variety and is a matrix with as many columns as plants per genotype.
# se_f_pop:        a list with the estimated population trajectories (se.fp_pop), and first (se.fp_pop_deriv1) and second (se.fp_pop_deriv2) order derivatives. Each element of the list is a numeric matrix where each column corresponds to one population.
# se_f_geno:       a list with the estimated genotype/variety trajectories (se.fp_geno), and first (se.fp_geno_deriv1) and second (se.fp_geno_deriv2) order derivatives. Each element of the list is a list that corresponds to one population and is a matrix with as many columns as genotypes/varieties per population.
# se_f_geno_dev:   a list with the estimated genotype/variety deviations (se.fp_geno_dev), and first (se.fp_geno_dev_deriv1) and second (se.fp_geno_dev_deriv2) order derivatives. Each element of the list is a list that corresponds to one population and is a matrix with as many columns as genotypes/varieties per population.
# se_f_plant:      a list with the estimated plant/plot/individual trajectories (se.fp_plant), and first (se.fp_plant_deriv1) and second (se.fp_plant_deriv2) order derivatives. Each element of the list is a list that corresponds to one genotype/variety and is a matrix with as many columns as plants per genotype.
# se_f_plant_dev:  a list with the estimated plant/plot/individual deviations (se.fp_plant_dev), and first (se.fp_plant_dev_deriv1) and second (se.fp_plant_dev_deriv2) order derivatives. Each element of the list is a list that corresponds to one genotype/variety and is a matrix with as many columns as plants per genotype.

# Function ----------------------------------------------------------------


#' predict.psHDM
#'
#' predict.psHDM
#'
#' @export
predict.psHDM <- function(object,
                          newtimes,
                          pred = list(pop = TRUE, geno = TRUE, plant = TRUE),
                          se = list(pop = TRUE, geno = TRUE, plant = FALSE),
                          ...) {

  if(!inherits(object, "psHDM")) {
    stop("The object class is not correct")
  }

  if(missing(newtimes)) {
    xp <- object$time
  } else {
    if(!is.vector(newtimes) | !is.numeric(newtimes)) {
      stop("newtimes should be a vector")
    }
    xp <- newtimes
  }

  # Output data
  res             <- list()
  res$newtimes     <- xp

  # Normalize time
  xp <- xp - min(object$time) + 1

  # Number of parameters: fixed and random (for each component)
  np        <- object$dim
  np.comp   <- c(np[1]+np[2], np[3]+np[4], np[5]+np[6])
  names(np.comp) <- c("pop", "geno", "plant")
  np.e      <- cumsum(np.comp)
  np.s      <- np.e - np.comp + 1

  # Predictions
  # Functions at population level
  if(isTRUE(pred$pop)) {
    # Transformation matrix, theta and B
    n.pop  <- length(object$l.pop)
    T_pop <- Matrix::Matrix(cbind(kronecker(diag(n.pop), object$MM$MM.pop$U.X),
                          kronecker(diag(n.pop), object$MM$MM.pop$U.Z)))
    theta_pop <- Matrix::crossprod(Matrix::t(T_pop), Matrix::Matrix(object$coeff[np.s[1]:np.e[1]], ncol = 1))
    Bp_pop <- Matrix::Matrix(kronecker(diag(n.pop), spline.bbase(knots = object$MM$MM.pop$knots, X. = xp, BDEG. = object$smooth$smooth.pop$bdeg)))
    fp_pop <- matrix(Matrix::crossprod(Matrix::t(Bp_pop), theta_pop), ncol = n.pop) # Note that f_pop == eta_pop
    colnames(fp_pop) <- object$l.pop

    # First derivative
    Bp_pop_deriv1 <- Matrix::Matrix(kronecker(diag(n.pop), spline.bbase(knots = object$MM$MM.pop$knots, X. = xp, BDEG. = object$smooth$smooth.pop$bdeg, deriv = 1)))
    fp_pop_deriv1 <- matrix(Matrix::crossprod(Matrix::t(Bp_pop_deriv1), theta_pop), ncol = n.pop)
    colnames(fp_pop_deriv1) <- object$l.pop

    # Second derivative
    Bp_pop_deriv2 <- kronecker(diag(n.pop), spline.bbase(knots = object$MM$MM.pop$knots, X. = xp, BDEG. = object$smooth$smooth.pop$bdeg, deriv = 2))
    fp_pop_deriv2 <- matrix(Matrix::crossprod(Matrix::t(Bp_pop_deriv2), theta_pop), ncol = n.pop)
    colnames(fp_pop_deriv2) <- object$l.pop

    print("Population-specific growth curves OK")

    res$f_pop <- list(fp_pop = fp_pop, fp_pop_deriv1 = fp_pop_deriv1, fp_pop_deriv2 = fp_pop_deriv2)
  }

  # Functions at genotype level
  if(isTRUE(pred$geno)) {
    # Genotype-specific deviations
    n.geno     <- length(object$l.geno)
    T_geno_dev <- Matrix::Matrix(cbind(kronecker(Matrix::Matrix(diag(n.geno)), object$MM$MM.geno$U.X),
                               kronecker(diag(n.geno), object$MM$MM.geno$U.Z)))
    theta_geno_dev <- Matrix::crossprod(Matrix::t(T_geno_dev), matrix(object$coeff[np.s[2]:np.e[2]], ncol = 1))
    Bp_geno_dev    <- Matrix::Matrix(kronecker(Matrix::Matrix(diag(n.geno)), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg)))
    fp_geno_dev    <- matrix(Matrix::crossprod(Matrix::t(Bp_geno_dev), theta_geno_dev), ncol = n.geno) # Note that f_geno_dev == eta_geno_dev
    # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per threatment
    ind.geno.pop   <- rep(1:n.pop, object$n.geno_p_pop)
    fp_geno_dev    <- lapply(split(fp_geno_dev, rep(ind.geno.pop, each = nrow(fp_geno_dev))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_geno_dev) <- object$l.pop

    # First derivative
    Bp_geno_dev_deriv1 <- Matrix::Matrix(kronecker(Matrix::Matrix(diag(n.geno)), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg, deriv = 1)))
    fp_geno_dev_deriv1 <- matrix(Matrix::crossprod(Matrix::t(Bp_geno_dev_deriv1), theta_geno_dev), ncol = n.geno)
    fp_geno_dev_deriv1 <- lapply(split(fp_geno_dev_deriv1, rep(ind.geno.pop, each = nrow(fp_geno_dev_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_geno_dev_deriv1) <- object$l.pop

    # Second derivative
    Bp_geno_dev_deriv2 <- Matrix::Matrix(kronecker(Matrix::Matrix(diag(n.geno)), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg, deriv = 2)))
    fp_geno_dev_deriv2 <- matrix(Matrix::crossprod(Matrix::t(Bp_geno_dev_deriv2), theta_geno_dev), ncol = n.geno)
    fp_geno_dev_deriv2 <- lapply(split(fp_geno_dev_deriv2, rep(ind.geno.pop, each = nrow(fp_geno_dev_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_geno_dev_deriv2) <- object$l.pop

    print("Genotype-specific deviations OK")

    # Genotype-specific growth curves
    if(n.pop == 1) {
      mm.ind.pop <- matrix(1, ncol = 1, nrow = n.geno)
    } else {
      ind.pop    <- rep(as.factor(1:n.pop), object$n.geno_p_pop)
      mm.ind.pop <- model.matrix(~ 0 + ind.pop) # The contrast matrix changes here!!!!!!
    }
    T_geno     <- Matrix::bdiag(T_pop, T_geno_dev)
    theta_geno <- Matrix::crossprod(Matrix::t(T_geno), matrix(object$coeff[np.s[1]:np.e[2]], ncol = 1))
    Bp_geno    <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg)), Bp_geno_dev)
    fp_geno    <- matrix(Matrix::crossprod(Matrix::t(Bp_geno), theta_geno), ncol = n.geno) # Note that f_geno == eta_geno
    # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per threatment
    fp_geno    <- lapply(split(fp_geno, rep(ind.geno.pop, each = nrow(fp_geno))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_geno) <- object$l.pop

    # First derivative
    Bp_geno_deriv1 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg, deriv = 1)), Bp_geno_dev_deriv1)
    fp_geno_deriv1 <- matrix(Matrix::crossprod(Matrix::t(Bp_geno_deriv1), theta_geno), ncol = n.geno)
    fp_geno_deriv1 <- lapply(split(fp_geno_deriv1, rep(ind.geno.pop, each = nrow(fp_geno_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_geno_deriv1) <- object$l.pop

    # Second derivative
    Bp_geno_deriv2 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg, deriv = 2)), Bp_geno_dev_deriv2)
    fp_geno_deriv2 <- matrix(Matrix::crossprod(Matrix::t(Bp_geno_deriv2), theta_geno), ncol = n.geno)
    fp_geno_deriv2 <- lapply(split(fp_geno_deriv2, rep(ind.geno.pop, each = nrow(fp_geno_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_geno_deriv2) <- object$l.pop

    print("Genotype-specific growth curves OK")

    # Give the name of the genotypes
    e <- cumsum(object$n.geno_p_pop)
    s <- e - object$n.geno_p_pop + 1
    for(i in 1:n.pop) {
      colnames(fp_geno_dev[[i]])        <- object$l.geno[s[i]:e[i]]
      colnames(fp_geno_dev_deriv1[[i]]) <- object$l.geno[s[i]:e[i]]
      colnames(fp_geno_dev_deriv2[[i]]) <- object$l.geno[s[i]:e[i]]
      colnames(fp_geno[[i]])            <- object$l.geno[s[i]:e[i]]
      colnames(fp_geno_deriv1[[i]])     <- object$l.geno[s[i]:e[i]]
      colnames(fp_geno_deriv2[[i]])     <- object$l.geno[s[i]:e[i]]
    }
    res$f_geno      <- list(fp_geno = fp_geno, fp_geno_deriv1 = fp_geno_deriv1, fp_geno_deriv2 = fp_geno_deriv2)
    res$f_geno_dev  <- list(fp_geno_dev = fp_geno_dev, fp_geno_dev_deriv1 = fp_geno_dev_deriv1, fp_geno_dev_deriv2 = fp_geno_dev_deriv2)
  }

  # Functions at plant level
  if(isTRUE(pred$plant)) {
    # Plant-specific deviations
    # Transformation matrix, theta and B
    n.tot        <- sum(object$n.plants_p_geno)
    T_plant_dev  <- cbind(kronecker(Matrix::Matrix(diag(n.tot)), object$MM$MM.ind$U.X),
                          kronecker(Matrix::Matrix(diag(n.tot)), object$MM$MM.ind$U.Z))
    theta_plant_dev     <- Matrix::Matrix(T_plant_dev%*%object$coeff[np.s[3]:np.e[3]]) #crossprod(t(T_plant_dev), matrix(object$coeff[np.s[3]:np.e[3]], ncol = 1))
    ind.ind.geno        <- rep(1:n.geno, object$n.plants_p_geno)
    Bp_plant_dev        <- Matrix::Matrix(kronecker(Matrix::Matrix(diag(n.tot)), spline.bbase(knots = object$MM$MM.ind$knots, X. = xp, BDEG. = object$smooth$smooth.plant$bdeg)))
    fp_plant_dev        <- matrix(Matrix::crossprod(Matrix::t(Bp_plant_dev), theta_plant_dev), ncol = n.tot) # Note that f_plant_dev == eta_plant_dev
    # List, with as many elements as genotypes. Each element of the list is a matrix, with as many columns as plants per genotype
    fp_plant_dev        <- lapply(split(fp_plant_dev, rep(ind.ind.geno, each = nrow(fp_plant_dev))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_plant_dev) <- object$l.geno

    # First derivative
    Bp_plant_dev_deriv1 <- kronecker(Matrix::Matrix(diag(n.tot)), spline.bbase(knots = object$MM$MM.ind$knots, X. = xp, BDEG. = object$smooth$smooth.plant$bdeg, deriv = 1))
    fp_plant_dev_deriv1 <- matrix(Matrix::crossprod(Matrix::t(Bp_plant_dev_deriv1), theta_plant_dev), ncol = n.tot)
    fp_plant_dev_deriv1 <- lapply(split(fp_plant_dev_deriv1, rep(ind.ind.geno, each = nrow(fp_plant_dev_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_plant_dev_deriv1) <- object$l.geno

    # Second derivative
    Bp_plant_dev_deriv2 <- kronecker(Matrix::Matrix(diag(n.tot)), spline.bbase(knots = object$MM$MM.ind$knots, X. = xp, BDEG. = object$smooth$smooth.plant$bdeg, deriv = 2))
    fp_plant_dev_deriv2 <- matrix(Matrix::crossprod(Matrix::t(Bp_plant_dev_deriv2), theta_plant_dev), ncol = n.tot)
    fp_plant_dev_deriv2 <- lapply(split(fp_plant_dev_deriv2, rep(ind.ind.geno, each = nrow(fp_plant_dev_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_plant_dev_deriv2) <- object$l.geno

    print("Plant-specific deviations OK")

    # Plant-specific growth curves
    if(n.pop == 1) {
      mm.ind.pop <- matrix(1, ncol = 1, nrow = object$n.plants_p_pop)
    } else {
      ind.pop    <- rep(as.factor(1:n.pop), object$n.plants_p_pop)
      mm.ind.pop <- model.matrix(~ 0 + ind.pop) # The contrast matrix changes here!!!!!!
    }
    if(n.geno == 1) {
      mm.ind.geno <- matrix(1, ncol = 1, nrow = n.geno)
    } else {
      ind.geno    <- rep(as.factor(1:n.geno), object$n.plants_p_geno)
      mm.ind.geno <- model.matrix(~ 0 + ind.geno) # The contrast matrix changes here!!!!!!
    }
    T_plant     <- Matrix::bdiag(T_pop, T_geno_dev, T_plant_dev)
    theta_plant <- Matrix::crossprod(Matrix::t(T_plant), matrix(object$coeff[np.s[1]:np.e[3]], ncol = 1))
    Bp_plant    <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = object$MM$MM.pop$knots, X. = xp, BDEG. = object$smooth$smooth.pop$bdeg)),
                         kronecker(Matrix::Matrix(mm.ind.geno), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg)),
                         Bp_plant_dev)
    fp_plant <- matrix(Matrix::crossprod(Matrix::t(Bp_plant), theta_plant), ncol = n.tot) # Note that f_plant == eta_plant
    # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per population
    fp_plant <- lapply(split(fp_plant, rep(ind.ind.geno, each = nrow(fp_plant))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_plant) <- object$l.geno

    # First derivative
    Bp_plant_deriv1 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = object$MM$MM.pop$knots, X. = xp, BDEG. = object$smooth$smooth.pop$bdeg, deriv = 1)),
                             kronecker(Matrix::Matrix(mm.ind.geno), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg, deriv = 1)),
                             Bp_plant_dev_deriv1)
    fp_plant_deriv1 <- matrix(Matrix::crossprod(Matrix::t(Bp_plant_deriv1), theta_plant), ncol = n.tot)
    fp_plant_deriv1 <- lapply(split(fp_plant_deriv1, rep(ind.ind.geno, each = nrow(fp_plant_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_plant_deriv1) <- object$l.geno

    # Second derivative
    Bp_plant_deriv2 <- cbind(kronecker(Matrix::Matrix(mm.ind.pop), spline.bbase(knots = object$MM$MM.pop$knots, X. = xp, BDEG. = object$smooth$smooth.pop$bdeg, deriv = 2)),
                             kronecker(Matrix::Matrix(mm.ind.geno), spline.bbase(knots = object$MM$MM.geno$knots, X. = xp, BDEG. = object$smooth$smooth.geno$bdeg, deriv = 2)),
                             Bp_plant_dev_deriv2)
    fp_plant_deriv2 <- matrix(Matrix::crossprod(Matrix::t(Bp_plant_deriv2), theta_plant), ncol = n.tot)
    fp_plant_deriv2 <- lapply(split(fp_plant_deriv2, rep(ind.ind.geno, each = nrow(fp_plant_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(fp_plant_deriv2) <- object$l.geno

    print("Plant-specific growth curves OK")

    res$f_plant     <- list(fp_plant = fp_plant, fp_plant_deriv1 = fp_plant_deriv1, fp_plant_deriv2 = fp_plant_deriv2)
    res$f_plant_dev <- list(fp_plant_dev = fp_plant_dev, fp_plant_dev_deriv1 = fp_plant_dev_deriv1, fp_plant_dev_deriv2 = fp_plant_dev_deriv2)
  }

  # Stardard errors to construct confidence intervals
  if(isTRUE(se$pop)) {
    # Population-specific growth curves
    se.theta_pop <- Matrix::crossprod(Matrix::t(T_pop), Matrix::tcrossprod(object$Vp[np.s[1]:np.e[1], np.s[1]:np.e[1]], T_pop))
    se.fp_pop    <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_pop) * Matrix::tcrossprod(se.theta_pop, Bp_pop))), ncol = n.pop)
    colnames(se.fp_pop) <- object$l.pop

    # Standard errors first derivative
    se.fp_pop_deriv1 <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_pop_deriv1) * Matrix::tcrossprod(se.theta_pop, Bp_pop_deriv1))), ncol = n.pop)
    colnames(se.fp_pop_deriv1) <- object$l.pop

    # Standard errors second derivative
    se.fp_pop_deriv2 <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_pop_deriv2) * Matrix::tcrossprod(se.theta_pop, Bp_pop_deriv2))), ncol = n.pop)
    colnames(se.fp_pop_deriv2) <- object$l.pop

    print("Standard errors for population-specific growth curves OK")

    res$se.f_pop  <- list(se.fp_pop = se.fp_pop, se.fp_pop_deriv1 = se.fp_pop_deriv1, se.fp_pop_deriv2 = se.fp_pop_deriv2)
  }

  if(isTRUE(se$geno)) {
    # Genotype-specific deviations
    se.theta_geno_dev <- Matrix::crossprod(Matrix::t(T_geno_dev), Matrix::tcrossprod(object$Vp[np.s[2]:np.e[2], np.s[2]:np.e[2]], T_geno_dev))
    se.fp_geno_dev    <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_geno_dev) * Matrix::tcrossprod(se.theta_geno_dev, Bp_geno_dev))), ncol = n.geno)
    se.fp_geno_dev    <- lapply(split(se.fp_geno_dev, rep(ind.geno.pop, each = nrow(se.fp_geno_dev))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_geno_dev) <- object$l.pop

    # Standard errors first derivative
    se.fp_geno_dev_deriv1 <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_geno_dev_deriv1) * Matrix::tcrossprod(se.theta_geno_dev, Bp_geno_dev_deriv1))), ncol = n.geno)
    se.fp_geno_dev_deriv1 <- lapply(split(se.fp_geno_dev_deriv1, rep(ind.geno.pop, each = nrow(se.fp_geno_dev_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_geno_dev_deriv1) <- object$l.pop

    # Standard errors second derivative
    se.fp_geno_dev_deriv2 <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_geno_dev_deriv2) * Matrix::tcrossprod(se.theta_geno_dev, Bp_geno_dev_deriv2))), ncol = n.geno)
    se.fp_geno_dev_deriv2 <- lapply(split(se.fp_geno_dev_deriv2, rep(ind.geno.pop, each = nrow(se.fp_geno_dev_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_geno_dev_deriv2) <- object$l.pop

    print("Standard errors for genotyp-specific deviations OK")

    # Genotype-specific growth curves
    se.theta_geno     <- Matrix::crossprod(Matrix::t(T_geno), Matrix::tcrossprod(object$Vp[np.s[1]:np.e[2], np.s[1]:np.e[2]], T_geno))
    se.fp_geno        <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_geno) * Matrix::tcrossprod(se.theta_geno, Bp_geno))), ncol = n.geno)
    se.fp_geno        <- lapply(split(se.fp_geno, rep(ind.geno.pop, each = nrow(se.fp_geno))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_geno) <- object$l.pop

    # Standard errors first derivative
    se.fp_geno_deriv1  <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_geno_deriv1) * Matrix::tcrossprod(se.theta_geno, Bp_geno_deriv1))), ncol = n.geno)
    se.fp_geno_deriv1  <- lapply(split(se.fp_geno_deriv1, rep(ind.geno.pop, each = nrow(se.fp_geno_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_geno_deriv1) <- object$l.pop

    # Standard errors second derivative
    se.fp_geno_deriv2  <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_geno_deriv2) * Matrix::tcrossprod(se.theta_geno, Bp_geno_deriv2))), ncol = n.geno)
    se.fp_geno_deriv2  <- lapply(split(se.fp_geno_deriv2, rep(ind.geno.pop, each = nrow(se.fp_geno_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_geno_deriv2) <- object$l.pop

    print("Standard errors for genotype-specific growth curves OK")

    # Give the name of the genotypes
    e <- cumsum(object$n.geno_p_pop)
    s <- e - object$n.geno_p_pop + 1
    for(i in 1:n.pop) {
      colnames(se.fp_geno_dev[[i]])        <- object$l.geno[s[i]:e[i]]
      colnames(se.fp_geno_dev_deriv1[[i]]) <- object$l.geno[s[i]:e[i]]
      colnames(se.fp_geno_dev_deriv2[[i]]) <- object$l.geno[s[i]:e[i]]
      colnames(se.fp_geno[[i]])            <- object$l.geno[s[i]:e[i]]
      colnames(se.fp_geno_deriv1[[i]])     <- object$l.geno[s[i]:e[i]]
      colnames(se.fp_geno_deriv2[[i]])     <- object$l.geno[s[i]:e[i]]
    }

    res$se.f_geno      <- list(se.fp_geno = se.fp_geno, se.fp_geno_deriv1 = se.fp_geno_deriv1, se.fp_geno_deriv2 = se.fp_geno_deriv2)
    res$se.f_geno_dev  <- list(se.fp_geno_dev = se.fp_geno_dev, se.fp_geno_dev_deriv1 = se.fp_geno_dev_deriv1, se.fp_geno_dev_deriv2 = se.fp_geno_dev_deriv2)
  }

  if(isTRUE(se$plant)) {
    # Plant-specific deviations
    se.theta_plant_dev     <- Matrix::crossprod(Matrix::t(T_plant_dev), Matrix::tcrossprod(object$Vp[np.s[3]:np.e[3], np.s[3]:np.e[3]], T_plant_dev))
    se.fp_plant_dev        <- matrix(sqrt(Matrix::colSums(t(Bp_plant_dev) * Matrix::tcrossprod(se.theta_plant_dev, Bp_plant_dev))), ncol = n.tot)
    se.fp_plant_dev        <- lapply(split(se.fp_plant_dev, rep(ind.ind.geno, each = nrow(se.fp_plant_dev))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_plant_dev) <- object$l.geno

    # Standard errors first derivative
    se.fp_plant_dev_deriv1  <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_plant_dev_deriv1) * Matrix::tcrossprod(se.theta_plant_dev, Bp_plant_dev_deriv1))), ncol = n.tot)
    se.fp_plant_dev_deriv1  <- lapply(split(se.fp_plant_dev_deriv1, rep(ind.ind.geno, each = nrow(se.fp_plant_dev_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_plant_dev_deriv1) <- object$l.geno

    # Standard errors second derivative
    se.fp_plant_dev_deriv2  <- matrix(sqrt(Matrix::colSums(Matrix::t(Bp_plant_dev_deriv2) * Matrix::tcrossprod(se.theta_plant_dev, Bp_plant_dev_deriv2))), ncol = n.tot)
    se.fp_plant_dev_deriv2  <- lapply(split(se.fp_plant_dev_deriv2, rep(ind.ind.geno, each = nrow(se.fp_plant_dev_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_plant_dev_deriv2) <- object$l.geno

    print("Standard errors for plant-specific deviations OK")

    # Standard errors for plant deviations + geno deviations + population effects
    se.theta_plant     <- Matrix::crossprod(Matrix::t(T_plant), Matrix::tcrossprod(object$Vp[np.s[1]:np.e[3], np.s[1]:np.e[3]], T_plant))
    se.fp_plant        <- matrix(sqrt(colSums(Matrix::t(Bp_plant) * Matrix::tcrossprod(se.theta_plant, Bp_plant))), ncol = n.tot)
    se.fp_plant        <- lapply(split(se.fp_plant, rep(ind.ind.geno, each = nrow(se.fp_plant))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_plant) <- object$l.geno

    # Standard errors first derivative
    se.fp_plant_deriv1  <- matrix(sqrt(colSums(t(Bp_plant_deriv1) * Matrix::tcrossprod(se.theta_plant, Bp_plant_deriv1))), ncol = n.tot)
    se.fp_plant_deriv1  <- lapply(split(se.fp_plant_deriv1, rep(ind.ind.geno, each = nrow(se.fp_plant_deriv1))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_plant_deriv1) <- object$l.geno

    # Standard errors second derivative
    se.fp_plant_deriv2  <- matrix(sqrt(colSums(Matrix::t(Bp_plant_deriv2) * Matrix::tcrossprod(se.theta_plant, Bp_plant_deriv2))), ncol = n.tot)
    se.fp_plant_deriv2  <- lapply(split(se.fp_plant_deriv2, rep(ind.ind.geno, each = nrow(se.fp_plant_deriv2))), function(x, nobs) matrix(x, nrow = nobs), nobs = length(xp))
    names(se.fp_plant_deriv2) <- object$l.geno

    print("Standard errors for plant-specific growth curves OK")

    res$se.f_plant     <- list(se.fp_plant = se.fp_plant, se.fp_plant_deriv1 = se.fp_plant_deriv1, se.fp_plant_deriv2 = se.fp_plant_deriv2)
    res$se.f_plant_dev <- list(se.fp_plant_dev = se.fp_plant_dev, se.fp_plant_dev_deriv1 = se.fp_plant_dev_deriv1, se.fp_plant_dev_deriv2 = se.fp_plant_dev_deriv2)
  }

  res$l.pop            <- object$l.pop
  res$l.geno           <- object$l.geno
  res$l.plant          <- object$l.plant
  res$n.plants_p_pop   <- object$n.plants_p_pop
  res$n.geno_p_pop     <- object$n.geno_p_pop
  res$n.plants_p_geno  <- object$n.plants_p_geno
  res$y                <- object$y
  res$raw.time         <- object$time
  class(res)           <- "pred.psHDM"
  return(res)
}




