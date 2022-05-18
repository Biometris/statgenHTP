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


# Help functions ----------------------------------------------------------
  mixmod_to_bspline_pred <- function(what = c("pop", "geno", "plant"), 
                                     object, np.s, np.e, xp,
                                     Tmat = NULL, mod.mat = NULL, dev = TRUE, 
                                     Bbasis = NULL){
    # Predictions
    # Transformation matrix Tm, Mixed model coefficients MM.coeff, theta and B, and fitted/predicted values f

    what.s <- what.e <- what
    if(isFALSE(dev) & what != "pop"){
      if(what == "geno"){
        what.s = "pop"
        what.e = "geno"
      }else {
        what.s = "pop"
        what.e = "plant"
      }
    }
    if(is.null(Tmat)) {
      Tm    <- cbind(Matrix::kronecker(Matrix::Diagonal(length(object[[paste0("l.",what)]])), object$MM[[paste0("MM.",what)]]$U.X),
                     Matrix::kronecker(Matrix::Diagonal(length(object[[paste0("l.",what)]])), object$MM[[paste0("MM.",what)]]$U.Z))
    } else{
      Tm <- Tmat
    }
    if(is.null(mod.mat)) {
      mmat <- Matrix::Diagonal(length(object[[paste0("l.",what)]]))
    }else if(is.list(mod.mat)){
      mmat <- list(mmat1 = mod.mat[[1]], mmat2 = mod.mat[[2]])
    }else{
      mmat <- mod.mat
    }
    MM.coeff <- matrix(object$coeff[np.s[what.s]:np.e[what.e]], ncol = 1)
    theta <- Matrix::Matrix(Tm %*% MM.coeff)
    B.full <- function(mmat, what, deriv){
      Matrix::kronecker(mmat, spline.bbase(knots = object$MM[[paste0("MM.",what)]]$knots, X. = xp, BDEG. = object$smooth[[paste0("smooth.",what)]]$bdeg, deriv = deriv))
    }
    if(is.null(Bbasis)) {
      B     <- B.full(mmat, what, deriv = 0)
      B.d1  <- B.full(mmat, what, deriv = 1)
      B.d2  <- B.full(mmat, what, deriv = 2)
    } else{
      if(what == "geno"){
        B    <- cbind(B.full(mmat, what = "pop", deriv = 0), Bbasis$B)
        B.d1 <- cbind(B.full(mmat, what = "pop", deriv = 1), Bbasis$B.d1)
        B.d2 <- cbind(B.full(mmat, what = "pop", deriv = 2), Bbasis$B.d2)
      }else if(what == "plant"){
        B    <- cbind(B.full(mmat$mmat1, what = "pop", deriv = 0), 
                      B.full(mmat$mmat2, what = "geno", deriv = 0),
                      Bbasis$B)
        B.d1 <- cbind(B.full(mmat$mmat1, what = "pop", deriv = 1), 
                      B.full(mmat$mmat2, what = "geno", deriv = 1),
                      Bbasis$B.d1)
        B.d2 <- cbind(B.full(mmat$mmat1, what = "pop", deriv = 2), 
                      B.full(mmat$mmat2, what = "geno", deriv = 2),
                      Bbasis$B.d2)
      }
    }
    f     <- matrix(Matrix::crossprod(Matrix::t(B), theta), ncol = length(object[[paste0("l.",what)]])) # Note that f == eta_what
    f.d1  <- matrix(Matrix::crossprod(Matrix::t(B.d1), theta), ncol = length(object[[paste0("l.",what)]])) 
    f.d2  <- matrix(Matrix::crossprod(Matrix::t(B.d2), theta), ncol = length(object[[paste0("l.",what)]])) 
    
    aux <- lapply(list(f = f, f.d1 = f.d1, f.d2 = f.d2), function(x){
      colnames(x) <- object[[paste0("l.",what)]]
      x
    })
    
    # Object to be returned
    res <- list(level = what,
                Tm = Tm,
                MM.coeff = MM.coeff,
                B = B,
                B.d1 = B.d1,
                B.d2 = B.d2,
                pred = aux)
    return(res)
  }
  standard_errors        <- function(Tm, B, B.d1, B.d2, 
                                     what = c("pop","geno","plant"), dev = TRUE, 
                                     object, np.s, np.e){
    what.s <- what.e <- what
    if(isFALSE(dev) & what != "pop"){
      if(what == "geno"){
        what.s = "pop"
        what.e = "geno"
      }else {
        what.s = "pop"
        what.e = "plant"
      }
    }
    se.theta <- Matrix::Matrix(Tm) %*% Matrix::Matrix(object$Vp[np.s[what.s]:np.e[what.e], np.s[what.s]:np.e[what.e]]) %*% Matrix::t(Tm)
    se.f     <- matrix(sqrt(Matrix::colSums(Matrix::t(B) * (se.theta %*% Matrix::t(B)))), ncol = length(object[[paste0("l.",what)]]))
    se.f.d1  <- matrix(sqrt(Matrix::colSums(Matrix::t(B.d1) * (se.theta %*% Matrix::t(B.d1)))), ncol = length(object[[paste0("l.",what)]]))
    se.f.d2  <- matrix(sqrt(Matrix::colSums(Matrix::t(B.d2) * (se.theta %*% Matrix::t(B.d2)))), ncol = length(object[[paste0("l.",what)]]))
    res <- list(se.f = se.f,
                se.f.d1 = se.f.d1,
                se.f.d2 = se.f.d2) 
    res <- lapply(res, function(x){
      colnames(x) <- object[[paste0("l.",what)]]
      x
    })
    return(res)
  }
  format_pred            <- function(what = c("geno", "plant"), 
                                     object, curves, xp){
    if (what == "geno") {
      my.each <- object$n.geno_p_pop
      my.x    <- object$l.pop
    } else {
      my.each <- object$n.plants_p_geno
      my.x    <- object$l.geno
    }
    res      <- tidyr::gather(as.data.frame(curves), what, f, factor_key = TRUE)
    res$time <- rep(xp, length(object[[paste0("l.", what)]]))
    res      <- lapply(split(res, f = rep(my.x, my.each*length(xp))), 
                       function(x) tidyr::spread(x, what, f)[,-1])
    return(res)
  }
  
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
  res$newtimes    <- xp
  
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
    pop_level <- mixmod_to_bspline_pred(what = "pop", object = object, np.s, np.e, xp, dev = FALSE)
    print("Population-specific growth curves OK")
    res$f_pop <- pop_level$pred
  }
  
  # Functions at genotype level
  if(isTRUE(pred$geno)) {
    # Genotype-specific deviations and first- and second-order derivatives
      geno_level <- mixmod_to_bspline_pred(what = "geno", object = object, np.s, np.e, xp, dev = TRUE)
      # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per population
      res$f_geno_dev <- lapply(geno_level$pred, format_pred, what = "geno", object = object, xp = xp)
      print("Genotype-specific deviations OK")
      
    # Genotype-specific growth curves and first- and second-order derivatives
      # Contrast matrix: Assign genotypes to populations
      if(length(object$l.pop) == 1) {
        mm.geno.pop <- matrix(1, ncol = 1, nrow = length(object$l.geno))
      } else {
        geno.pop    <- rep(as.factor(1:length(object$l.pop)), object$n.geno_p_pop)
        mm.geno.pop <- Matrix::sparse.model.matrix(~ 0 + geno.pop) # The contrast matrix changes here!!!!!!
      }
      T_geno     <- Matrix::bdiag(pop_level$Tm, geno_level$Tm)
      geno_tra_level <- mixmod_to_bspline_pred(what = "geno", object = object, np.s, np.e, xp, Tmat = T_geno, mod.mat = mm.geno.pop, dev = FALSE, 
                                               Bbasis = list(B = geno_level$B, B.d1 = geno_level$B.d1, B.d2 = geno_level$B.d2))
      # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per population
      res$f_geno <- lapply(geno_tra_level$pred, format_pred, what = "geno", object = object, xp = xp)
      print("Genotype-specific growth curves OK")
  }
  
  # Functions at plant level
  if(isTRUE(pred$plant)) {
    # Plant-specific deviations and first- and second-order derivatives
      plant_level <- mixmod_to_bspline_pred(what = "plant", object = object, np.s, np.e, xp, dev = TRUE)
      # List, with as many elements as genotypes. Each element of the list is a matrix, with as many columns as plants per genotype
      res$f_plant_dev <- lapply(plant_level$pred, format_pred, what = "plant", object = object, xp = xp)
      print("Plant-specific deviations OK")
      
    # Plant-specific growth curves
      # Contrast matrix: Assign plants to populations
      if(length(object$l.pop) == 1) {
        mm.ind.pop <- matrix(1, ncol = 1, nrow = object$n.plants_p_pop)
      } else {
        ind.pop    <- rep(as.factor(1:length(object$l.pop)), object$n.plants_p_pop)
        mm.ind.pop <- Matrix::sparse.model.matrix(~ 0 + ind.pop) # The contrast matrix changes here!!!!!!
      }
      # Contrast matrix: Assign plants to genotypes
      if(length(object$l.geno) == 1) {
        mm.ind.geno <- matrix(1, ncol = 1, nrow = length(object$l.geno))
      } else {
        ind.geno    <- rep(as.factor(1:length(object$l.geno)), object$n.plants_p_geno)
        mm.ind.geno <- Matrix::sparse.model.matrix(~ 0 + ind.geno) # The contrast matrix changes here!!!!!!
      }
      T_plant     <- Matrix::bdiag(pop_level$Tm, geno_level$Tm, plant_level$Tm)
      plant_tra_level <- mixmod_to_bspline_pred(what = "plant", object = object, np.s, np.e, xp, Tmat = T_plant, mod.mat = list(mm.ind.pop, mm.ind.geno), dev = FALSE, 
                                               Bbasis = list(B = plant_level$B, B.d1 = plant_level$B.d1, B.d2 = plant_level$B.d2))
      # List, with as many elements as populations. Each element of the list is a matrix, with as many columns as genotypes per population
      res$f_plant <- lapply(plant_tra_level$pred, format_pred, what = "plant", object = object, xp = xp)
      print("Plant-specific growth curves OK")
  }
  
  # Stardard errors to construct confidence intervals
  if(isTRUE(se$pop)) {
    # Population-specific growth curves
    res$se.f_pop <- standard_errors(Tm = pop_level$Tm, B = pop_level$B, B.d1 = pop_level$B.d1, B.d2 = pop_level$B.d2, what = "pop", dev = FALSE, object = object, np.s, np.e)
    print("Standard errors for population-specific growth curves OK")
  }
  
  if(isTRUE(se$geno)) {
    # Genotype-specific deviations
      res$se.f_geno_dev <- standard_errors(Tm = geno_level$Tm, B = geno_level$B, B.d1 = geno_level$B.d1, B.d2 = geno_level$B.d2, what = "geno", dev = TRUE, object = object, np.s, np.e)
      res$se.f_geno_dev <- lapply(res$se.f_geno_dev, format_pred, what = "geno", object = object, xp = xp)
      print("Standard errors for genotype-specific deviations OK")
    
    # Genotype-specific growth curves
      res$se.f_geno <- standard_errors(Tm = T_geno, B = geno_tra_level$B, B.d1 = geno_tra_level$B.d1, B.d2 = geno_tra_level$B.d2, what = "geno", dev = FALSE, object = object, np.s, np.e)
      res$se.f_geno <- lapply(res$se.f_geno, format_pred, what = "geno", object = object, xp = xp)
      print("Standard errors for genotype-specific growth curves OK")
  }
  
  if(isTRUE(se$plant)) {
    # Plant-specific deviations
      res$se.f_plant_dev <- standard_errors(Tm = plant_level$Tm, B = plant_level$B, B.d1 = plant_level$B.d1, B.d2 = plant_level$B.d2, what = "plant", dev = TRUE, object = object, np.s, np.e)
      res$se.f_plant_dev <- lapply(res$se.f_plant_dev, format_pred, what = "plant", object = object, xp = xp)
      print("Standard errors for plant-specific deviations OK")
    
    # Standard errors for plant deviations + geno deviations + population effects
      res$se.f_plant <- standard_errors(Tm = T_plant, B = plant_tra_level$B, B.d1 = plant_tra_level$B.d1, B.d2 = plant_tra_level$B.d2, what = "plant", dev = FALSE, object = object, np.s, np.e)
      res$se.f_plant <- lapply(res$se.f_plant, format_pred, what = "plant", object = object, xp = xp)
      print("Standard errors for plant-specific growth curves OK")
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
