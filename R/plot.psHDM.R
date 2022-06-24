# Function ----------------------------------------------------------------

#' plot.pred.psHDM
#'
#' This plot function provides five plots for objects of the class \code{psHDM}
#' after fitting (\code{\link{fitSplineHDM}}) or predicting
#' (\code{\link{predict.psHDM}}): (1) Population-specific growth curves,
#' (2) Genotype-specific growth curves (for all genotypes), (3) First-order
#' derivative of the genotype-specific growth curves (for all genotypes),
#' (4) Genotype-specific deviations (for all genotypes), and (5) Plant- and
#' genotype-specific growth curves (for a selection of genotypes). If standard
#' error are available, 95% pointwise confidence intervals are depicted.
#'
#' @param x An object of class "psHDM" as obtained after fitting
#' (\code{\link{fitSplineHDM}}) or predicting (\code{\link{predict.psHDM}}),
#' @param ... Not used.
#' @param geno.sub A character vector with the genotypes for which plots are
#' desired.
#' @param geno.sub.names A character vector with the names for the genotypes
#' selected (geno.sub). The default is \code{NULL}
#' @param geno.sub.order A vector with the order of the genotypes selected
#' (geno.sub). The default is \code{NULL}
#' @param xlab The x label of the plot. The default is \code{NULL}
#' @param ylab The y label of the plot. The default is \code{NULL}
#' @param my.theme Theme to be used for the ggplots. The default is \code{my.theme()}.
#' @param global.main A list with overall titles for plots at each level of the
#' hierarchy. The default is \code{NULL} for all of them.
#' @param ask Should the next figure be drawn?
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
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#' ## Plots at plant level for some genotypes (as illustration)
#' plot(fit.psHDM,
#'     geno.sub = c("GenoA14_WD","GenoA51_WD","GenoB11_WW","GenoB02_WD","GenoB02_WW"),
#'     my.theme = my.theme())
#'
#' @references Pérez-Valencia, D.M., Rodríguez-Álvarez, M.X., Boer, M.P. et al.
#' A two-stage approach for the spatio-temporal analysis of high-throughput
#' phenotyping data. Sci Rep 12, 3177 (2022). \doi{10.1038/s41598-022-06935-9}
#'
#' @export
plot.psHDM <- function(x,
                       ...,
                       geno.sub,
                       geno.sub.names = NULL,
                       geno.sub.order = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       my.theme = my.theme(),
                       ask = TRUE,
                       global.main = list(pop.tra = NULL,
                                          geno.tra = NULL,
                                          geno.tra.deriv1 = NULL,
                                          geno.dev = NULL,
                                          plant.tra = NULL)){

  if (!inherits(x, "psHDM")) {
    stop("The object class is not correct")
  }
  if (is.null(xlab)) {
    xlab <- "Time"
  }
  if (is.null(ylab)) {
    ylab <- expression(tilde(y)[pgi](t))
  }
  if(is.null(x$plant.obs)){
    x$plant.obs <- x$plant.level
  }
  ntime <- length(unique(x$plant.obs$timeNumber))
  min.t <- min(x$plant.obs$timeNumber)
  max.t <- max(x$plant.obs$timeNumber)
  if (!is.null(x$pop.level)) {
    ## Population-specific growth curves.
    my.cols <- c("1" = "gray", "2" = "blue")
    aa1 <- ggplot2::ggplot(data = x$plant.obs) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, obs_plant, group = plant,
                                      colour = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$pop.level,
                         ggplot2::aes(timeNumber, f_pop, group = pop,
                                      colour = "2"), na.rm = TRUE) +
      ggplot2::geom_rug(ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols,
                                  labels = c(expression(tilde(y)[pgi](t)),
                                             expression(hat(f)[p](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = global.main$pop.tra, color = "") +
      my.theme() +
      ggplot2::facet_grid( ~ pop)
    if (!is.null(x$pop.level$se_pop)) {
      aa1 <- aa1 +
        ggplot2::geom_ribbon(data = x$pop.level,
                             ggplot2::aes(x = timeNumber,
                                          ymin = f_pop - 1.96 * se_pop,
                                          ymax = f_pop + 1.96 * se_pop,
                                          group = pop),
                             fill = ggplot2::alpha("blue", 0.3))
    }
    print(aa1)
  }
  if (!is.null(x$geno.level)) {
    ## Genotype-specific growth curves.
    my.cols <- c("1" = "blue", "2" = "orange")
    bb1 <- ggplot2::ggplot(data = x$geno.level) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, f_geno, group = geno,
                                      colour = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$pop.level,
                         ggplot2::aes(timeNumber, f_pop, group = pop,
                                      colour="2"), size = 0.8, na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plant.obs,
                        ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols,
                                  labels = c(expression(hat(f)[p](t)+hat(f)[pg](t)),
                                             expression(hat(f)[p](t)))) +
      ggplot2::labs(x = xlab, y = ylab, color = "",
                    title = global.main$geno.tra) +
      my.theme() +
      ggplot2::facet_grid( ~ pop)
    if (!is.null(x$pop.level$se_pop)) {
      bb1 <- bb1 +
        ggplot2::geom_ribbon(data = x$pop.level,
                             ggplot2::aes(x = timeNumber,
                                          ymin = f_pop - 1.96 * se_pop,
                                          ymax = f_pop + 1.96 * se_pop,
                                          group = pop),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    if (ask) readline("Press return for next page....")
    print(bb1)
    ## First derivative of the genotype-specific growth curves.
    bb2 <- ggplot2::ggplot(data = x$geno.level) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, f_geno_deriv1, group = geno,
                                      colour = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$pop.level,
                         ggplot2::aes(timeNumber, f_pop_deriv1, group = pop,
                                      colour = "2"), na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plant.obs,
                        ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::ylim(min(x$geno.level$f_geno_deriv1),
                    max(x$geno.level$f_geno_deriv1) * 1.2) +
      ggplot2::scale_color_manual(values = my.cols,
                                  labels = c(expression((hat(f)[p](t)+hat(f)[pg](t))*minute),
                                             expression((hat(f)[p](t))*minute))) +
      ggplot2::labs(x = xlab, y = "First derivative", color = "",
                    title = global.main$geno.tra.deriv1) +
      my.theme() +
      ggplot2::facet_grid( ~ pop)
    if (!is.null(x$pop.level$se_pop_deriv1)) {
      bb2 <- bb2 +
        ggplot2::geom_ribbon(data = x$pop.level,
                             ggplot2::aes(x = timeNumber,
                                          ymin = f_pop_deriv1 - 1.96 * se_pop_deriv1,
                                          ymax = f_pop_deriv1 + 1.96 * se_pop_deriv1,
                                          group = pop),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    if(ask) readline("Press return for next page....")
    print(bb2)
    ## Estimated genotypic deviations.
    my.cols <- c("1" = "blue")
    cc1 <- ggplot2::ggplot(data = x$geno.level) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, f_geno_dev, group = geno,
                                      colour = "1"), na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plant.obs,
                        ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = my.cols,
                                  labels = c(expression(hat(f)[pg](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = global.main$geno.dev,
                    color = "") +
      my.theme() +
      ggplot2::facet_grid( ~ pop)
    if(ask) readline("Press return for next page....")
    print(cc1)
  }
  ## Estimated plant trajectories.
  if (!is.null(x$plant.level)) {
    df.tra.HDM.sub <- droplevels(x$plant.level[x$plant.level$geno %in% geno.sub, ])
    df.obs.HDM.sub <- droplevels(x$plant.obs[x$plant.obs$geno %in% geno.sub, ])
    df.tra.HDM.geno.sub <- droplevels(x$geno.level[x$geno.level$geno %in% geno.sub, ])
    if (!is.null(geno.sub.names)) {
      df.tra.HDM.sub$geno      <- factor(df.tra.HDM.sub$geno[drop = TRUE],
                                         labels = geno.sub.names)
      df.obs.HDM.sub$geno      <- factor(df.obs.HDM.sub$geno[drop = TRUE],
                                         labels = geno.sub.names)
      df.tra.HDM.geno.sub$geno <- factor(df.tra.HDM.geno.sub$geno[drop = TRUE],
                                         labels = geno.sub.names)
    }
    if(!is.null(geno.sub.order)){
      df.tra.HDM.sub$geno      <- factor(df.tra.HDM.sub$geno,
                                         levels = levels(df.tra.HDM.sub$geno)[geno.sub.order])
      df.obs.HDM.sub$geno      <- factor(df.obs.HDM.sub$geno,
                                         levels = levels(df.obs.HDM.sub$geno)[geno.sub.order])
      df.tra.HDM.geno.sub$geno <- factor(df.tra.HDM.geno.sub$geno,
                                         levels = levels(df.tra.HDM.geno.sub$geno)[geno.sub.order])
    }
    my.cols <- c("1" = "gray", "2" = "blue", "3" = "red")
    dd1 <- ggplot2::ggplot(data = df.obs.HDM.sub) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, obs_plant, group = plant,
                                      colour = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = df.tra.HDM.sub,
                         ggplot2::aes(timeNumber, f_plant, group = plant,
                                      colour = "2"), linetype = 2, na.rm = TRUE) +
      ggplot2::geom_line(data = df.tra.HDM.geno.sub,
                         ggplot2::aes(timeNumber, f_geno, group = geno,
                                      colour = "3"), na.rm = TRUE) +
      ggplot2::geom_rug(data = df.obs.HDM.sub,
                        ggplot2::aes(x = timeNumber, y = NULL),
                        color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(min.t, max.t, length.out = 5), 0)) +
      ggplot2::ylim(min(df.obs.HDM.sub$obs_plant, na.rm = TRUE),
                    max(df.obs.HDM.sub$obs_plant, na.rm = TRUE) * 1.2) +
      ggplot2::scale_color_manual(values = my.cols,
                                  labels = c(expression(tilde(y)[pgi](t)),
                                             expression(hat(f)[p](t)+hat(f)[pg](t)+hat(f)[pgi](t)),
                                             expression(hat(f)[p](t)+hat(f)[pg](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = global.main$plant.tra,
                    color = "") +
      my.theme() +
      ggplot2::facet_grid( ~ geno)
    if (!is.null(x$geno.level$se_geno)) {
      dd1 <- dd1 +
        ggplot2::geom_ribbon(data = df.tra.HDM.geno.sub,
                             ggplot2::aes(x = timeNumber,
                                          ymin = f_geno - 1.96 * se_geno,
                                          ymax = f_geno + 1.96 * se_geno,
                                          group = pop),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    if(ask) readline("Press return for next page....")
    print(dd1)
  }
}

# Needed functions for plotting --------------------------------------------------------
# Theme for the ggplot plots ----------------------------------------------

#' my.theme
#'
#' my.theme
#'
#' @export
my.theme <- function(my.size = 20) {
  ggplot2::theme(strip.text.x = ggplot2::element_text(size = my.size + 5),
                 strip.text.y = ggplot2::element_text(size = my.size + 5),
                 plot.title = ggplot2::element_text(hjust = 0.5,
                                                    size = my.size + 5,
                                                    face = "bold"),
                 axis.text = ggplot2::element_text(size = my.size),
                 axis.title = ggplot2::element_text(size = my.size + 5),
                 legend.title = ggplot2::element_text(size = my.size),
                 legend.text = ggplot2::element_text(size = my.size + 5,
                                                     hjust = 0),
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                 legend.key = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                 legend.justification = c(1, 1),
                 legend.position = "top",
                 panel.background = ggplot2::element_rect(fill = "gray95"))
}
