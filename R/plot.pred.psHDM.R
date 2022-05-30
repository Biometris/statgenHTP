# Function to PLOT the P-spline Hierarchical Curve Data Model -------------
# Aim ---------------------------------------------------------------------
# This plot function provides six plots:
# Population-specific growth curves
# First-order derivative of the population-specific growth curves
# Genotype-specific growth curves (for all genotypes)
# First-order derivative of the genotype-specific growth curves (for all genotypes)
# Genotype-specific deviations (for all genotypes)
# Plant- and genotype-specific growth curves (for a selection of genotypes)

# Input arguments ---------------------------------------------------------
# object:         an object of class "pred.psHDM" as obtained using function predict.HDMss()
# geno.sub:       a character vector with the genotypes for which plots are desired.
# geno.sub.names: a character vector with the names for the genotypes selected (geno.sub). The default is NULL.
# geno.sub.order: a vector with the order of the genotypes selected (geno.sub). The default is NULL.
# xlab:           the x label of the plot. The default is NULL.
# ylab:           the y label of the plot. The default is NULL.
# my.theme:       theme to be used for the ggplots. The default is my.theme().
# global.main:    overall title for plots at each level of the hierarchy. The default is NULL for all of them.
# ask:            a logical value. If TRUE, the default, the user is asked for confirmation, before a new figure is drawn.

# Function ----------------------------------------------------------------

#' plot.pred.psHDM
#'
#' plot.pred.psHDM
#'
#' @export
plot.pred.psHDM <- function(object,
                            geno.sub,
                            geno.sub.names = NULL,
                            geno.sub.order = NULL,
                            xlab = NULL,
                            ylab = NULL,
                            my.theme = my.theme(),
                            ask = TRUE,
                            global.main = list(pop.tra = NULL,
                                               pop.tra.deriv1 = NULL,
                                               geno.tra = NULL,
                                               geno.tra.deriv1 = NULL,
                                               geno.dev = NULL,
                                               plant.tra = NULL)){

  if(!inherits(object, "pred.psHDM")) {
    stop("The object class is not correct")
  }

  if(is.null(xlab)){
    xlab = "Time"
  }

  if(is.null(ylab)){
    ylab = expression(tilde(y)[pgi](t))
  }

  # object as data frames
  dataDF <- list.to.df(object = object)

  ntime <- length(unique(dataDF$plant.obs$timepoint))
  min.t <- min(dataDF$plant.obs$timepoint)
  max.t <- max(dataDF$plant.obs$timepoint)

  if(!is.null(object$f_pop)){
    # Population-specific growth curves
    my.cols <- c("1" = "gray", "2" = "blue", "3" = "blue")
    aa1 <- ggplot2::ggplot(data = dataDF$plant.obs) +
      ggplot2::geom_line(data = dataDF$plant.obs[!is.na(dataDF$plant.obs$obs_plant),],
                         ggplot2::aes(timepoint, obs_plant, group = plant, colour = "1")) +
      ggplot2::geom_point(ggplot2::aes(timepoint, obs_plant, group = plant, colour = "1"), shape = 20) +
      ggplot2::geom_line(data = dataDF$pop.tra, ggplot2::aes(timepoint, eta_pop, group = pop, colour="2")) +
      ggplot2::geom_ribbon(data = dataDF$pop.tra, ggplot2::aes(x = timepoint, ymin = eta_pop - 1.96 * se_pop, ymax = eta_pop + 1.96 * se_pop, group = pop), fill = ggplot2::alpha("blue", 0.3)) +
      ggplot2::geom_rug(ggplot2::aes(x = timepoint, y = NULL), color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols, labels = c(expression(tilde(y)[pgi](t)),expression(hat(f)[p](t)),"95% CI"))+
      ggplot2::labs(x = xlab, y = ylab, title = global.main$pop.tra, color = "") +
      my.theme + ggplot2::theme(legend.position = c(0.08,1.05)) +
      ggplot2::facet_grid( ~ pop)
    print(aa1)

    # First derivative of the population-specific growth curves
    my.cols <- c("1" = "blue", "2" = "blue")
    aa2 <- ggplot2::ggplot(data = dataDF$pop.tra) +
      ggplot2::geom_line(ggplot2::aes(timepoint, eta_pop_deriv1, group = pop, colour = "1")) +
      ggplot2::geom_ribbon(ggplot2::aes(x = timepoint, ymin = eta_pop_deriv1 - 1.96 * se_pop_deriv1, ymax = eta_pop_deriv1 + 1.96 * se_pop_deriv1, group = pop), fill = ggplot2::alpha("blue", 0.3)) +
      ggplot2::geom_rug(data = dataDF$plant.obs, ggplot2::aes(x = timepoint, y = NULL), color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols, labels = c(expression(hat(f)*minute[p](t)),"95% CI"))+
      ggplot2::labs(x = xlab, y = "First derivative", color = "", title = global.main$pop.tra.deriv1) +
      my.theme + ggplot2::theme(legend.position = c(0.06,1.05)) +
      ggplot2::facet_grid( ~ pop)

    if(ask) readline("Press return for next page....")
    print(aa2)
  }

  if(!is.null(object$f_geno)){
    # Genotype-specific growth curves
    my.cols <- c("1" = "blue", "2" = "orange")
    bb1 <- ggplot2::ggplot(data = dataDF$geno.tra) +
      ggplot2::geom_line(ggplot2::aes(timepoint, eta_geno, group = geno, colour="1")) +
      ggplot2::geom_line(data = dataDF$pop.tra, ggplot2::aes(timepoint, eta_pop, group = pop, colour="2"), size = 0.8) +
      ggplot2::geom_ribbon(data = dataDF$pop.tra, ggplot2::aes(x = timepoint, ymin = eta_pop - 1.96 * se_pop, ymax = eta_pop + 1.96 * se_pop, group = pop), fill = ggplot2::alpha("orange", 0.5)) +
      ggplot2::geom_rug(data = dataDF$plant.obs, ggplot2::aes(x = timepoint, y = NULL), color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols, labels = c(expression(hat(f)[p](t)+hat(f)[pg](t)),expression(hat(f)[p](t))))+
      ggplot2::labs(x = xlab, y = ylab, color = "", title = global.main$geno.tra) +
      my.theme + ggplot2::theme(legend.position = c(0.12,1.05)) +
      ggplot2::facet_grid( ~ pop)
    if(ask) readline("Press return for next page....")
    print(bb1)

    # First derivative of the genotype-specific growth curves
    my.cols <- c("1" = "blue", "2" = "orange")
    bb2 <- ggplot2::ggplot(data = dataDF$geno.tra) +
      ggplot2::geom_line(ggplot2::aes(timepoint, eta_geno_deriv1, group = geno, colour = "1")) +
      ggplot2::geom_line(data = dataDF$pop.tra, ggplot2::aes(timepoint, eta_pop_deriv1, group = pop, colour = "2")) +
      ggplot2::geom_ribbon(data = dataDF$pop.tra, ggplot2::aes(x = timepoint, ymin = eta_pop_deriv1 - 1.96 * se_pop_deriv1, ymax = eta_pop_deriv1 + 1.96 * se_pop_deriv1, group = pop), fill = ggplot2::alpha("orange", 0.5)) +
      ggplot2::geom_rug(data = dataDF$plant.obs, ggplot2::aes(x = timepoint, y = NULL), color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::ylim(min(dataDF$geno.tra$eta_geno_deriv1), max(dataDF$geno.tra$eta_geno_deriv1)*1.2) +
      ggplot2::scale_color_manual(values = my.cols, labels = c(expression((hat(f)[p](t)+hat(f)[pg](t))*minute),expression((hat(f)[p](t))*minute)))+
      ggplot2::labs(x = xlab, y = "First derivative", color = "", title = global.main$geno.tra.deriv1) +
      my.theme + ggplot2::theme(legend.position = c(0.13,1.05)) +
      ggplot2::facet_grid( ~ pop)
    if(ask) readline("Press return for next page....")
    print(bb2)

    # Estimated genotypic deviations
    my.cols <- c("1" = "blue")
    cc1 <- ggplot2::ggplot(data = dataDF$geno.tra) +
      ggplot2::geom_line(ggplot2::aes(timepoint, eta_geno_dev, group = geno, colour="1")) +
      ggplot2::geom_rug(data = dataDF$plant.obs, ggplot2::aes(x = timepoint, y = NULL), color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols, labels = c(expression(hat(f)[pg](t))))+
      ggplot2::labs(x = xlab, y = ylab, title = global.main$geno.dev ,color = "") +
      my.theme + ggplot2::theme(legend.position = c(0.07,1.05)) +
      ggplot2::facet_grid( ~ pop)
    if(ask) readline("Press return for next page....")
    print(cc1)
  }

  # Estimated plant trajectories
  if(!is.null(object$f_plant)) {
    df.tra.HDM.sub           <- droplevels(dataDF$plant.tra[dataDF$plant.tra$geno %in% geno.sub,])
    df.obs.HDM.sub           <- droplevels(dataDF$plant.obs[dataDF$plant.obs$geno %in% geno.sub,])
    df.tra.HDM.geno.sub      <- droplevels(dataDF$geno.tra[dataDF$geno.tra$geno %in% geno.sub,])

    if(!is.null(geno.sub.names)){
      df.tra.HDM.sub$geno      <- factor(df.tra.HDM.sub$geno[drop=TRUE], labels = geno.sub.names)
      df.obs.HDM.sub$geno      <- factor(df.obs.HDM.sub$geno[drop=TRUE], labels = geno.sub.names)
      df.tra.HDM.geno.sub$geno <- factor(df.tra.HDM.geno.sub$geno[drop=TRUE], labels = geno.sub.names)
    }
    if(!is.null(geno.sub.order)){
      df.tra.HDM.sub$geno      <- factor(df.tra.HDM.sub$geno, levels = levels(df.tra.HDM.sub$geno)[geno.sub.order])
      df.obs.HDM.sub$geno      <- factor(df.obs.HDM.sub$geno, levels = levels(df.obs.HDM.sub$geno)[geno.sub.order])
      df.tra.HDM.geno.sub$geno <- factor(df.tra.HDM.geno.sub$geno, levels = levels(df.tra.HDM.geno.sub$geno)[geno.sub.order])
    }

    my.cols <- c("1" = "gray", "2" = "blue", "3" = "orange")
    dd1 <- ggplot2::ggplot(data = df.obs.HDM.sub) +
      ggplot2::geom_ribbon(data = df.tra.HDM.geno.sub, ggplot2::aes(x = timepoint, ymin = eta_geno - 1.96 * se_geno, ymax = eta_geno + 1.96 * se_geno, group = pop), fill = ggplot2::alpha("orange", 0.5)) +
      ggplot2::geom_line(data = df.obs.HDM.sub[!is.na(df.obs.HDM.sub$obs_plant),], ggplot2::aes(timepoint, obs_plant, group = plant, colour = "1")) +
      ggplot2::geom_point(ggplot2::aes(timepoint, obs_plant, group = plant, colour = "1"), shape = 20) +
      ggplot2::geom_line(data = df.tra.HDM.sub, ggplot2::aes(timepoint, eta_plant, group = plant, colour = "2"), linetype = 2) +
      ggplot2::geom_line(data = df.tra.HDM.geno.sub, ggplot2::aes(timepoint, eta_geno, group = geno, colour = "3")) +
      ggplot2::geom_rug(data = df.obs.HDM.sub, ggplot2::aes(x = timepoint, y = NULL), color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      # scale_y_continuous(breaks = round(seq(min(df.obs.HDM.sub$obs_plant, na.rm = TRUE), max(df.obs.HDM.sub$obs_plant, na.rm = TRUE)*1.5, length.out = 5), 2)) +
      ggplot2::ylim(min(df.obs.HDM.sub$obs_plant, na.rm = TRUE), max(df.obs.HDM.sub$obs_plant, na.rm = TRUE)*1.2) +
      ggplot2::scale_color_manual(values = my.cols, labels = c(expression(tilde(y)[pgi](t)),expression(hat(f)[p](t)+hat(f)[pg](t)+hat(f)[pgi](t)),expression(hat(f)[p](t)+hat(f)[pg](t))))+
      ggplot2::labs(x = xlab, y = ylab, title = global.main$plant.tra, color = "") +
      my.theme + ggplot2::theme(legend.position = c(0.17,1.05)) +
      ggplot2::facet_grid( ~ geno)
    if(ask) readline("Press return for next page....")
    print(dd1)
  }
}

# Needed functions for plotting --------------------------------------------------------
# Lists to data frames  ---------------------------------------------------
list.to.df <- function(object) {
  if(!inherits(object, "pred.psHDM")) {
    stop("The object class is not correct")
  }
  
  list.to.vec <- function(x) {c(do.call("cbind", lapply(x, function(i) as.matrix(i))))}
  
  res <- list()
  # Population-specific growth curves
  if(!is.null(object$f_pop)) {
    df.pop.tra <- data.frame(eta_pop = c(object$f_pop$f),
                             eta_pop_deriv1 = c(object$f_pop$f.d1),
                             eta_pop_deriv2 = c(object$f_pop$f.d2),
                             se_pop = c(object$se.f_pop$se.f),
                             se_pop_deriv1 = c(object$se.f_pop$se.f.d1),
                             se_pop_deriv2 = c(object$se.f_pop$se.f.d2),
                             pop = rep(object$l.pop, each = length(object$newtimes)),
                             timepoint = rep(object$newtimes, length(object$l.pop)))
    res$pop.tra   <- df.pop.tra
  }

  # Genotypic-specific growth curves and deviations
  if(!is.null(object$f_geno)) {
    df.geno.tra <- data.frame(eta_geno = list.to.vec(object$f_geno$f),
                              eta_geno_deriv1 = list.to.vec(object$f_geno$f.d1),
                              eta_geno_deriv2 = list.to.vec(object$f_geno$f.d2),
                              se_geno = list.to.vec(object$se.f_geno$se.f),
                              se_geno_deriv1 = list.to.vec(object$se.f_geno$se.f.d1),
                              se_geno_deriv2 = list.to.vec(object$se.f_geno$se.f.d2),
                              eta_geno_dev = list.to.vec(object$f_geno_dev$f),
                              eta_geno_dev_deriv1 = list.to.vec(object$f_geno_dev$f.d1),
                              eta_geno_dev_deriv2 = list.to.vec(object$f_geno_dev$f.d2),
                              se_geno_dev = list.to.vec(object$se.f_geno_dev$se.f),
                              se_geno_dev_deriv1 = list.to.vec(object$se.f_geno_dev$se.f.d1),
                              se_geno_dev_deriv2 = list.to.vec(object$se.f_geno_dev$se.f.d2),
                              pop = rep(object$l.pop, object$n.geno_p_pop*length(object$newtimes)),
                              geno = rep(object$l.geno, each = length(object$newtimes)),
                              timepoint = rep(object$newtimes, length(object$l.geno)))
    res$geno.tra  <- df.geno.tra
  }

  # Plant-specific growth curves and deviations
  if(!is.null(object$f_plant)) {
    df.plant.tra <- data.frame(eta_plant = list.to.vec(object$f_plant$f),
                               eta_plant_deriv1 = list.to.vec(object$f_plant$f.d1),
                               eta_plant_deriv2 = list.to.vec(object$f_plant$f.d2),
                               eta_plant_dev = list.to.vec(object$f_plant_dev$f),
                               eta_plant_dev_deriv1 = list.to.vec(object$f_plant_dev$f.d1),
                               eta_plant_dev_deriv2 = list.to.vec(object$f_plant_dev$f.d2),
                               pop = rep(object$l.pop, object$n.plants_p_pop*length(object$newtimes)),
                               geno = rep(object$l.geno, object$n.plants_p_geno*length(object$newtimes)),
                               plant = rep(object$l.plant, each = length(object$newtimes)),
                               timepoint = rep(object$newtimes, sum(object$n.plants_p_geno)))
    res$plant.tra <- df.plant.tra
  }

  # Raw Plant growth curves (Observed data: in the raw.time not in newtimes)
  df.plant.obs <- data.frame(obs_plant = c(do.call("cbind", object$y)),
                             pop = rep(object$l.pop, object$n.plants_p_pop*length(object$raw.time)),
                             geno = rep(object$l.geno, object$n.plants_p_geno*length(object$raw.time)),
                             plant = rep(object$l.plant, each = length(object$raw.time)),
                             timepoint = rep(object$raw.time, sum(object$n.plants_p_geno)))
  res$plant.obs <- df.plant.obs

  class(res)    <- "HDMss"
  res
}

# Theme for the ggplot plots ----------------------------------------------

#' my.theme
#'
#' my.theme
#'
#' @export
my.theme <- function(my.size = 20) {
  ggplot2::theme(strip.text.x = ggplot2::element_text(size = my.size + 5),
                 strip.text.y = ggplot2::element_text(size = my.size + 5),
                 plot.title = ggplot2::element_text(hjust = 0.5, size = my.size + 5, face = "bold"),
                 axis.text = ggplot2::element_text(size = my.size),
                 axis.title = ggplot2::element_text(size = my.size + 5),
                 legend.title = ggplot2::element_text(size = my.size),
                 legend.text = ggplot2::element_text(size = my.size + 5, hjust = 0),
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                 legend.key = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                 legend.justification = c(1,1),
                 panel.background = ggplot2::element_rect(fill = "gray95"))
}
