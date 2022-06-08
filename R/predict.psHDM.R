# Help functions ----------------------------------------------------------
mixmod_to_bspline_pred <- function(what = c("pop", "geno", "plant"),
                                   object,
                                   np.s,
                                   np.e,
                                   xp,
                                   Tmat = NULL,
                                   mod.mat = NULL,
                                   dev = TRUE,
                                   Bbasis = NULL){
  ## Predictions.
  ## Transformation matrix Tm, Mixed model coefficients MM.coeff,
  ## theta and B, and fitted/predicted values f.
  what.s <- what.e <- what
  if (isFALSE(dev) && what != "pop") {
    if (what == "geno") {
      what.s <- "pop"
      what.e <- "geno"
    } else {
      what.s <- "pop"
      what.e <- "plant"
    }
  }
  if (is.null(Tmat)) {
    Tm <- cbind(Matrix::kronecker(Matrix::Diagonal(length(object[[paste0("l.", what)]])),
                                  object$MM[[paste0("MM.", what)]]$U.X),
                Matrix::kronecker(Matrix::Diagonal(length(object[[paste0("l.", what)]])),
                                  object$MM[[paste0("MM.",what)]]$U.Z))
  } else{
    Tm <- Tmat
  }
  if (is.null(mod.mat)) {
    mmat <- Matrix::Diagonal(length(object[[paste0("l.",what)]]))
  } else if (is.list(mod.mat)) {
    mmat <- list(mmat1 = mod.mat[[1]], mmat2 = mod.mat[[2]])
  } else{
    mmat <- mod.mat
  }
  MM.coeff <- matrix(object$coeff[np.s[what.s]:np.e[what.e]], ncol = 1)
  theta <- Matrix::Matrix(Tm %*% MM.coeff)
  B.full <- function(mmat, what, deriv) {
    Matrix::kronecker(mmat,
                      spline.bbase(knots = object$MM[[paste0("MM.",what)]]$knots,
                                   X. = xp,
                                   BDEG. = object$smooth[[paste0("smooth.", what)]]$bdeg,
                                   deriv = deriv))
  }
  if (is.null(Bbasis)) {
    B     <- B.full(mmat, what, deriv = 0)
    B.d1  <- B.full(mmat, what, deriv = 1)
    B.d2  <- B.full(mmat, what, deriv = 2)
  } else{
    if (what == "geno"){
      B    <- cbind(B.full(mmat, what = "pop", deriv = 0), Bbasis$B)
      B.d1 <- cbind(B.full(mmat, what = "pop", deriv = 1), Bbasis$B.d1)
      B.d2 <- cbind(B.full(mmat, what = "pop", deriv = 2), Bbasis$B.d2)
    } else if(what == "plant") {
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
  f     <- matrix(Matrix::crossprod(Matrix::t(B), theta),
                  ncol = length(object[[paste0("l.",what)]])) # Note that f == eta_what
  f.d1  <- matrix(Matrix::crossprod(Matrix::t(B.d1), theta),
                  ncol = length(object[[paste0("l.",what)]]))
  f.d2  <- matrix(Matrix::crossprod(Matrix::t(B.d2), theta),
                  ncol = length(object[[paste0("l.",what)]]))

  aux <- lapply(list(f = f, f.d1 = f.d1, f.d2 = f.d2), function(x) {
    colnames(x) <- object[[paste0("l.", what)]]
    return(x)
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

standard_errors <- function(Tm,
                            B,
                            B.d1,
                            B.d2,
                            what = c("pop","geno","plant"),
                            dev = TRUE,
                            object,
                            np.s,
                            np.e) {
  what.s <- what.e <- what
  if (isFALSE(dev) && what != "pop"){
    if (what == "geno") {
      what.s <- "pop"
      what.e <- "geno"
    } else {
      what.s <- "pop"
      what.e <- "plant"
    }
  }
  se.theta <- Matrix::Matrix(Tm) %*%
    Matrix::Matrix(object$Vp[np.s[what.s]:np.e[what.e],
                             np.s[what.s]:np.e[what.e]]) %*% Matrix::t(Tm)
  se.f     <- matrix(sqrt(Matrix::colSums(Matrix::t(B) *
                                            (se.theta %*% Matrix::t(B)))),
                     ncol = length(object[[paste0("l.",what)]]))
  se.f.d1  <- matrix(sqrt(Matrix::colSums(Matrix::t(B.d1) *
                                            (se.theta %*% Matrix::t(B.d1)))),
                     ncol = length(object[[paste0("l.",what)]]))
  se.f.d2  <- matrix(sqrt(Matrix::colSums(Matrix::t(B.d2) *
                                            (se.theta %*% Matrix::t(B.d2)))),
                     ncol = length(object[[paste0("l.",what)]]))
  res <- list(se.f = se.f,
              se.f.d1 = se.f.d1,
              se.f.d2 = se.f.d2)
  res <- lapply(res, function(x){
    colnames(x) <- object[[paste0("l.", what)]]
    return(x)
  })
  return(res)
}

list.to.df <- function(object1,
                       object2,
                       what,
                       xp) {
  if(!inherits(object2, "psHDM")) {
    stop("The object class is not correct")
  }
  # list.to.vec <- function(x) {
  #   c(do.call("cbind", lapply(X = x, FUN = function(i) as.matrix(i))))
  # }
  res <- list()
  ## Population-specific growth curves.
  if (what == "pop") {
    df.pop.tra <- data.frame(timeNumber = rep(xp, length(object2$l.pop)),
                             pop = rep(object2$l.pop, each = length(xp)),
                             f_pop = c(object1$f),
                             f_pop_deriv1 = c(object1$f.d1),
                             f_pop_deriv2 = c(object1$f.d2))
    if (!is.null(object1$se.f)) {
      df.pop.tra$se_pop <- c(object1$se.f)
      df.pop.tra$se_pop_deriv1 <- c(object1$se.f.d1)
      df.pop.tra$se_pop_deriv2 <- c(object1$se.f.d2)
    }
    res$pop.level <- df.pop.tra
  }
  ## Genotypic-specific growth curves and deviations
  if (what == "geno") {
    df.geno.tra <- data.frame(timeNumber = rep(xp, length(object2$l.geno)),
                              pop = rep(object2$l.pop,
                                        object2$n.geno_p_pop * length(xp)),
                              geno = rep(object2$l.geno, each = length(xp)),
                              f_geno = as.vector(object1$geno_tra$f),
                              f_geno_deriv1 = as.vector(object1$geno_tra$f.d1),
                              f_geno_deriv2 = as.vector(object1$geno_tra$f.d2),
                              f_geno_dev = as.vector(object1$geno_dev$f),
                              f_geno_dev_deriv1 = as.vector(object1$geno_dev$f.d1),
                              f_geno_dev_deriv2 = as.vector(object1$geno_dev$f.d2))
    if (!is.null(object1$geno_dev$se.f)) {
      df.geno.tra$se_geno <- as.vector(object1$geno_tra$se.f)
      df.geno.tra$se_geno_deriv1 <- as.vector(object1$geno_tra$se.f.d1)
      df.geno.tra$se_geno_deriv2 <- as.vector(object1$geno_tra$se.f.d2)
      df.geno.tra$se_geno_dev <- as.vector(object1$geno_dev$se.f)
      df.geno.tra$se_geno_dev_deriv1 <- as.vector(object1$geno_dev$se.f.d1)
      df.geno.tra$se_geno_dev_deriv2 <- as.vector(object1$geno_dev$se.f.d2)
    }
    res$geno.level  <- df.geno.tra
  }
  ## Plant-specific growth curves and deviations
  if (what == "plant") {
    df.plant.tra <- data.frame(timeNumber = rep(xp, sum(object2$n.plants_p_geno)),
                               pop = rep(object2$l.pop,
                                         object2$n.plants_p_pop * length(xp)),
                               geno = rep(object2$l.geno,
                                          object2$n.plants_p_geno * length(xp)),
                               plant = rep(object2$l.plant, each = length(xp)),
                               f_plant = as.vector(object1$plant_tra$f),
                               f_plant_deriv1 = as.vector(object1$plant_tra$f.d1),
                               f_plant_deriv2 = as.vector(object1$plant_tra$f.d2),
                               f_plant_dev = as.vector(object1$plant_dev$f),
                               f_plant_dev_deriv1 = as.vector(object1$plant_dev$f.d1),
                               f_plant_dev_deriv2 = as.vector(object1$plant_dev$f.d2))
    if(!is.null(object1$se.f_plant)){
      df.plant.tra$se_plant <- as.vector(object1$plant_tra$se.f$se.f)
      df.plant.tra$se_plant_deriv1 <- as.vector(object1$plant_tra$se.f.d1)
      df.plant.tra$se_plant_deriv2 <- as.vector(object1$plant_tra$se.f.d2)
      df.plant.tra$se_plant_dev <- as.vector(object1$plant_dev$se.f)
      df.plant.tra$se_plant_dev_deriv1 <- as.vector(object1$plant_dev$se.f.d1)
      df.plant.tra$se_plant_dev_deriv2 <- as.vector(object1$plant_dev$se.f.d2)
    }
    res$plant.level <- df.plant.tra

    if (isTRUE(setequal(xp, object2$time))) {
      res$plant.level$obs_plant <- c(do.call("cbind", object2$y))
    } else {
      # Raw plant growth curves (Observed data: in the raw time not in newtimes (xp))
      df.plant.obs <-
        data.frame(timeNumber = rep(object2$time, sum(object2$n.plants_p_geno)),
                   pop = rep(object2$l.pop,
                             object2$n.plants_p_pop * length(object2$time)),
                   geno = rep(object2$l.geno,
                              object2$n.plants_p_geno * length(object2$time)),
                   plant = rep(object2$l.plant, each = length(object2$time)),
                   obs_plant = c(do.call("cbind", object2$y)))
      res$plant.obs <- df.plant.obs
    }
  }
  res
}

# Function ----------------------------------------------------------------

#' predict.psHDM
#'
#' Function that predicts the P-spline Hierarchical Curve Data Model (see
#' \code{\link{fitSplineHDM}}) on a dense grid. It provides standard errors
#' for curves at each level of the hierarchy. User has to be aware that
#' standard errors at the plant level demand large memory. We suggest set
#' that option at the \code{FALSE} level
#'
#' @param object An object of class "psHDM" as obtained after fitting
#' (\code{\link{fitSplineHDM}}) the P-spline Hierarchical Curve Data Model
#' @param newtimes A numeric vector with timepoints at which predictions are
#' desired
#' @param pred A list that controls the hierarchical levels at which
#' predictions are desired (population/genotypes/plants).  The default is
#' \code{TRUE}.
#' @param se A list that controls the hierarchical levels at which standard
#' errors are desired (population/genotypes/plants).  The default is
#' \code{TRUE} except at the plant level.
#' @param ... Not used.
#'
#' @return An object of class \code{psHDM}, a list with the following outputs:
#' predict.psHDM
#' \code{newtimes} A numeric vector with the timepoints at which predictions
#' and/or standard errors have been obtained.
#' \code{pop.level} A data.frame with the estimated population trajectories
#' and first and second order derivatives, and if required their respective
#' standard errors, at the \code{newtimes}.
#' \code{geno.level} A data.frame with the estimated genotype-specific
#' deviations and trajectories and their respective first and second order
#' derivatives, and if required their respective standard errors,
#' at the \code{newtimes}.
#' \code{plant.level} A data.frame with the estimated plant-specific
#' deviations and trajectories and their respective first and second order
#' derivatives, and if required their respective standard errors,
#' at the \code{newtimes}.
#' \code{plant.obs} A data.frame with the raw data at the original timepoints.
#'
#' #' @examples
#' ## Predict the P-Splines Hierarchical Curve Data Model in a dense grid
#' ## with standard errors at the population and genotype levels
#' pred.psHDM <- predict.psHDM(object = fit.psHDM,
#'                             newtimes = seq(min(fit.psHDM$time),
#'                                            max(fit.psHDM$time), length = 100),
#'                             pred = list(pop = TRUE, geno = TRUE, plant = TRUE),
#'                             se = list(pop = TRUE, geno = TRUE, plant = FALSE))
#'
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#' ## Plots at plant level for some genotypes (as illustration)
#' plot(object = pred.psHDM,
#'     geno.sub = c("GenoA14_WD","GenoA51_WD","GenoB11_WW","GenoB02_WD","GenoB02_WW"),
#'     my.theme = my.theme())
#'
#' @export
predict.psHDM <- function(object,
                          newtimes,
                          pred = list(pop = TRUE, geno = TRUE, plant = TRUE),
                          se = list(pop = TRUE, geno = TRUE, plant = FALSE),
                          trace = TRUE,
                          ...) {

  if (!inherits(object, "psHDM")) {
    stop("The object class is not correct")
  }
  if (missing(newtimes)) {
    xp <- object$time
  } else {
    if (!is.vector(newtimes) || !is.numeric(newtimes)) {
      stop("newtimes should be a vector")
    }
    xp <- newtimes
  }
  ## Output data.
  res <- list(newtimes = xp)
  ## Normalize time
  xp <- xp - min(xp) + 1
  ## Number of parameters: fixed and random (for each component)
  np        <- object$dim
  np.comp   <- c(np[1] + np[2], np[3] + np[4], np[5] + np[6])
  names(np.comp) <- c("pop", "geno", "plant")
  np.e      <- cumsum(np.comp)
  np.s      <- np.e - np.comp + 1
  ## Predictions
  ## Functions at population level
  if (isTRUE(pred$pop)) {
    pop_level <- mixmod_to_bspline_pred(what = "pop", object = object,
                                        np.s, np.e, xp, dev = FALSE)
    if (trace) {
      print("Population-specific growth curves OK")
    }
    if (isTRUE(se$pop)) {
      pop_level$pred <-
        append(pop_level$pred,
               standard_errors(Tm = pop_level$Tm, B = pop_level$B,
                               B.d1 = pop_level$B.d1, B.d2 = pop_level$B.d2,
                               what = "pop", dev = FALSE, object = object,
                               np.s, np.e))
      if (trace) {
        print("Standard errors for population-specific growth curves OK")
      }
    }
    ## Data frame with all the information at population level
    res <- append(res, list.to.df(object1 = pop_level$pred, object2 = object,
                                  what = "pop", xp))
  }
  ## Functions at genotype level.
  if (isTRUE(pred$geno)) {
    ## Genotype-specific deviations and first- and second-order derivatives
    geno_dev <- mixmod_to_bspline_pred(what = "geno", object = object,
                                       np.s, np.e, xp, dev = TRUE)
    geno_level <- list(geno_dev = geno_dev$pred)
    if (trace) {
      print("Genotype-specific deviations OK")
    }
    ## Genotype-specific growth curves and first- and second-order derivatives
    ## Contrast matrix: Assign genotypes to populations
    if(length(object$l.pop) == 1) {
      mm.geno.pop <- matrix(1, ncol = 1, nrow = length(object$l.geno))
    } else {
      geno.pop    <- rep(as.factor(1:length(object$l.pop)), object$n.geno_p_pop)
      mm.geno.pop <- Matrix::sparse.model.matrix(~ 0 + geno.pop) # The contrast matrix changes here!!!!!!
    }
    T_geno   <- Matrix::bdiag(pop_level$Tm, geno_dev$Tm)
    geno_tra <- mixmod_to_bspline_pred(what = "geno", object = object, np.s,
                                       np.e, xp, Tmat = T_geno,
                                       mod.mat = mm.geno.pop, dev = FALSE,
                                       Bbasis = list(B = geno_dev$B,
                                                     B.d1 = geno_dev$B.d1,
                                                     B.d2 = geno_dev$B.d2))
    geno_level$geno_tra <- geno_tra$pred
    if (trace) {
      print ("Genotype-specific growth curves OK")
    }
    if (isTRUE(se$geno)) {
      ## Genotype-specific deviations
      geno_level$geno_dev <-
        append(geno_level$geno_dev,
               standard_errors(Tm = geno_dev$Tm, B = geno_dev$B,
                               B.d1 = geno_dev$B.d1, B.d2 = geno_dev$B.d2,
                               what = "geno", dev = TRUE, object = object,
                               np.s, np.e))
      if (trace) {
        print("Standard errors for genotype-specific deviations OK")
      }
      ## Genotype-specific growth curves
      geno_level$geno_tra <-
        append(geno_level$geno_tra,
               standard_errors(Tm = T_geno, B = geno_tra$B,
                               B.d1 = geno_tra$B.d1, B.d2 = geno_tra$B.d2,
                               what = "geno", dev = FALSE, object = object,
                               np.s, np.e))
      if (trace) {
        print("Standard errors for genotype-specific growth curves OK")
      }
    }
    ## Data frame with all the information at genotype level.
    res <- append(res, list.to.df(object1 = geno_level, object2 = object,
                                  what = "geno", xp))
  }
  ## Functions at plant level.
  if (isTRUE(pred$plant)) {
    # Plant-specific deviations and first- and second-order derivatives
    plant_dev <- mixmod_to_bspline_pred(what = "plant", object = object, np.s,
                                        np.e, xp, dev = TRUE)
    plant_level <- list(plant_dev = plant_dev$pred)
    if (trace) {
      print("Plant-specific deviations OK")
    }
    ## Plant-specific growth curves.
    ## Contrast matrix: Assign plants to populations.
    if (length(object$l.pop) == 1) {
      mm.ind.pop <- matrix(1, ncol = 1, nrow = object$n.plants_p_pop)
    } else {
      ind.pop    <- rep(as.factor(1:length(object$l.pop)), object$n.plants_p_pop)
      mm.ind.pop <- Matrix::sparse.model.matrix(~ 0 + ind.pop) # The contrast matrix changes here!!!!!!
    }
    ## Contrast matrix: Assign plants to genotypes.
    if (length(object$l.geno) == 1) {
      mm.ind.geno <- matrix(1, ncol = 1, nrow = length(object$l.geno))
    } else {
      ind.geno    <- rep(as.factor(1:length(object$l.geno)), object$n.plants_p_geno)
      mm.ind.geno <- Matrix::sparse.model.matrix(~ 0 + ind.geno) # The contrast matrix changes here!!!!!!
    }
    T_plant   <- Matrix::bdiag(pop_level$Tm, geno_dev$Tm, plant_dev$Tm)
    plant_tra <- mixmod_to_bspline_pred(what = "plant", object = object, np.s,
                                        np.e, xp, Tmat = T_plant,
                                        mod.mat = list(mm.ind.pop, mm.ind.geno),
                                        dev = FALSE,
                                        Bbasis = list(B = plant_dev$B,
                                                      B.d1 = plant_dev$B.d1,
                                                      B.d2 = plant_dev$B.d2))
    plant_level$plant_tra <- plant_tra$pred
    if (trace) {
      print("Plant-specific growth curves OK")
    }
    if (isTRUE(se$plant)) {
      ## Plant-specific deviations.
      plant_level$plant_dev <-
        append(plant_level$plant_dev,
               standard_errors(Tm = plant_dev$Tm, B = plant_dev$B,
                               B.d1 = plant_dev$B.d1, B.d2 = plant_dev$B.d2,
                               what = "plant", dev = TRUE, object = object,
                               np.s, np.e))
      if (trace) {
        print("Standard errors for plant-specific deviations OK")
      }
      ## Standard errors for plant deviations + geno deviations + population effects.
      plant_level$plant_tra <-
        append(plant_level$plant_tra,
               standard_errors(Tm = T_plant, B = plant_tra$B,
                               B.d1 = plant_tra$B.d1, B.d2 = plant_tra$B.d2,
                               what = "plant", dev = FALSE, object = object,
                               np.s, np.e))
      if (trace) {
        print("Standard errors for plant-specific growth curves OK")
      }
    }
    ## Data frame with all the information at genotype level.
    res <- append(res,
                  list.to.df(object1 = plant_level, object2 = object,
                             what = "plant", xp))
  }
  class(res) <- "psHDM"
  return(res)
}
