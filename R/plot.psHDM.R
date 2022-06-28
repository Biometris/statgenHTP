# Function ----------------------------------------------------------------

#' plot.pred.psHDM
#'
#' This plot function provides five plots for objects of the class \code{psHDM}
#' after fitting (\code{\link{fitSplineHDM}}) or predicting
#' (\code{\link{predict.psHDM}}): (1) Population-specific growth curves,
#' (2) Genotype-specific growth curves (for all genotypes), (3) First-order
#' derivative of the genotype-specific growth curves (for all genotypes),
#' (4) Genotype-specific deviations (for all genotypes), and (5) Plot- and
#' genotype-specific growth curves (for a selection of genotypes). If standard
#' error are available, 95% pointwise confidence intervals are depicted.
#'
#' @param x An object of class "psHDM" as obtained after fitting
#' (\code{\link{fitSplineHDM}}) or predicting (\code{\link{predict.psHDM}}),
#' @param ... Not used.
#' @param genotypes A character vector with the genotypes for which plots are
#' desired.
#' @param genotypeNames A character vector with the names for the genotypes
#' selected (genotypes). The default is \code{NULL}
#' @param genotypeOrder A vector with the order of the genotypes selected
#' (genotypes). The default is \code{NULL}
#' @param xlab The x label of the plot. The default is \code{NULL}
#' @param ylab The y label of the plot. The default is \code{NULL}
#' @param themeHDM Theme to be used for the ggplots. The default is \code{themeHDM()}.
#' @param title A list with overall titles for plots at each level of the
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
#'                           genotype = "genobytreat",
#'                           plotId = "plotId",
#'                           difVar = list(geno = FALSE, plant = FALSE),
#'                           smoothPop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothGeno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothPlot = list(nseg = 4, bdeg = 3, pord = 2),
#'                           weights = "wt")
#'
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#' ## Plots at plant level for some genotypes (as illustration)
#' plot(fit.psHDM,
#'     genotypes = c("GenoA14_WD","GenoA51_WD","GenoB11_WW","GenoB02_WD","GenoB02_WW"),
#'     themeHDM = themeHDM())
#'
#' @references Pérez-Valencia, D.M., Rodríguez-Álvarez, M.X., Boer, M.P. et al.
#' A two-stage approach for the spatio-temporal analysis of high-throughput
#' phenotyping data. Sci Rep 12, 3177 (2022). \doi{10.1038/s41598-022-06935-9}
#'
#' @export
plot.psHDM <- function(x,
                       ...,
                       genotypes,
                       genotypeNames = NULL,
                       genotypeOrder = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       themeHDM = themeHDM(),
                       ask = TRUE,
                       title = list(popTra = NULL,
                                    genoTra = NULL,
                                    genoTraDeriv1 = NULL,
                                    genoDev = NULL,
                                    plotTra = NULL)){

  if (!inherits(x, "psHDM")) {
    stop("The object class is not correct")
  }
  if (is.null(xlab)) {
    xlab <- "Time"
  }
  if (is.null(ylab)) {
    ylab <- expression(tilde(y)[pgi](t))
  }
  if(is.null(x$plotObs)){
    x$plotObs <- x$plotLevel
  }
  ntime <- length(unique(x$plotObs$timeNumber))
  minT <- min(x$plotObs$timeNumber)
  maxT <- max(x$plotObs$timeNumber)
  if (!is.null(x$popLevel)) {
    ## Population-specific growth curves.
    plotCols <- c("1" = "gray", "2" = "blue")
    aa1 <- ggplot2::ggplot(data = x$plotObs) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, obsPlot, group = plotId,
                                      color = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$popLevel,
                         ggplot2::aes(timeNumber, fPop, group = pop,
                                      color = "2"), na.rm = TRUE) +
      ggplot2::geom_rug(ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(tilde(y)[pgi](t)),
                                             expression(hat(f)[p](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = title$popTra, color = "") +
      themeHDM() +
      ggplot2::facet_grid(~ pop)
    if (!is.null(x$popLevel$sePop)) {
      aa1 <- aa1 +
        ggplot2::geom_ribbon(data = x$popLevel,
                             ggplot2::aes(x = timeNumber,
                                          ymin = fPop - 1.96 * sePop,
                                          ymax = fPop + 1.96 * sePop,
                                          group = pop),
                             fill = ggplot2::alpha("blue", 0.3))
    }
    print(aa1)
  }
  if (!is.null(x$genoLevel)) {
    ## Genotype-specific growth curves.
    plotCols <- c("1" = "blue", "2" = "orange")
    bb1 <- ggplot2::ggplot(data = x$genoLevel) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, fGeno, group = genotype,
                                      color = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$popLevel,
                         ggplot2::aes(timeNumber, fPop, group = pop,
                                      color="2"), size = 0.8, na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plotObs,
                        ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(hat(f)[p](t)+hat(f)[pg](t)),
                                             expression(hat(f)[p](t)))) +
      ggplot2::labs(x = xlab, y = ylab, color = "",
                    title = title$genoTra) +
      themeHDM() +
      ggplot2::facet_grid(~ pop)
    if (!is.null(x$popLevel$sePop)) {
      bb1 <- bb1 +
        ggplot2::geom_ribbon(data = x$popLevel,
                             ggplot2::aes(x = timeNumber,
                                          ymin = fPop - 1.96 * sePop,
                                          ymax = fPop + 1.96 * sePop,
                                          group = pop),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    if (ask) readline("Press return for next page....")
    print(bb1)
    ## First derivative of the genotype-specific growth curves.
    bb2 <- ggplot2::ggplot(data = x$genoLevel) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, fGenoDeriv1, group = genotype,
                                      color = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$popLevel,
                         ggplot2::aes(timeNumber, fPopDeriv1, group = pop,
                                      color = "2"), na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plotObs,
                        ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT, length.out = 5), 0)) +
      ggplot2::ylim(min(x$genoLevel$fGenoDeriv1),
                    max(x$genoLevel$fGenoDeriv1) * 1.2) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression((hat(f)[p](t)+hat(f)[pg](t))*minute),
                                             expression((hat(f)[p](t))*minute))) +
      ggplot2::labs(x = xlab, y = "First derivative", color = "",
                    title = title$genoTraDeriv1) +
      themeHDM() +
      ggplot2::facet_grid(~ pop)
    if (!is.null(x$popLevel$sePopDeriv1)) {
      bb2 <- bb2 +
        ggplot2::geom_ribbon(data = x$popLevel,
                             ggplot2::aes(x = timeNumber,
                                          ymin = fPopDeriv1 - 1.96 * sePopDeriv1,
                                          ymax = fPopDeriv1 + 1.96 * sePopDeriv1,
                                          group = pop),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    if(ask) readline("Press return for next page....")
    print(bb2)
    ## Estimated genotypic deviations.
    plotCols <- c("1" = "blue")
    cc1 <- ggplot2::ggplot(data = x$genoLevel) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, fGenoDev, group = genotype,
                                      color = "1"), na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plotObs,
                        ggplot2::aes(x = timeNumber, y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT, length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(hat(f)[pg](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = title$genoDev,
                    color = "") +
      themeHDM() +
      ggplot2::facet_grid(~ pop)
    if(ask) readline("Press return for next page....")
    print(cc1)
  }
  ## Estimated plot trajectories.
  if (!is.null(x$plotLevel)) {
    dfTraSub <- droplevels(x$plotLevel[x$plotLevel$genotype %in% genotypes, ])
    dfObsSub <- droplevels(x$plotObs[x$plotObs$genotype %in% genotypes, ])
    dfTraGenoSub <- droplevels(x$genoLevel[x$genoLevel$genotype %in% genotypes, ])
    if (!is.null(genotypeNames)) {
      dfTraSub$genotype <- factor(dfTraSub$genotype[drop = TRUE],
                                  labels = genotypeNames)
      dfObsSub$genotype <- factor(dfObsSub$genotype[drop = TRUE],
                                  labels = genotypeNames)
      dfTraGenoSub$genotype <- factor(dfTraGenoSub$geno[drop = TRUE],
                                      labels = genotypeNames)
    }
    if (!is.null(genotypeOrder)) {
      dfTraSub$genotype <- factor(dfTraSub$genotype,
                                  levels = levels(dfTraSub$genotype)[genotypeOrder])
      dfObsSub$genotype <- factor(dfObsSub$genotype,
                                  levels = levels(dfObsSub$genotype)[genotypeOrder])
      dfTraGenoSub$genotype <- factor(dfTraGenoSub$genotype,
                                      levels = levels(dfTraGenoSub$genotype)[genotypeOrder])
    }
    plotCols <- c("1" = "gray", "2" = "blue", "3" = "red")
    dd1 <- ggplot2::ggplot(data = dfObsSub) +
      ggplot2::geom_line(ggplot2::aes(timeNumber, obsPlot, group = plotId,
                                      color = "1"), na.rm = TRUE) +
      ggplot2::geom_line(data = dfTraSub,
                         ggplot2::aes(timeNumber, fPlot, group = plotId,
                                      color = "2"), linetype = 2, na.rm = TRUE) +
      ggplot2::geom_line(data = dfTraGenoSub,
                         ggplot2::aes(timeNumber, fGeno, group = genotype,
                                      color = "3"), na.rm = TRUE) +
      ggplot2::geom_rug(data = dfObsSub,
                        ggplot2::aes(x = timeNumber, y = NULL),
                        color = "gray", length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT, length.out = 5), 0)) +
      ggplot2::ylim(min(dfObsSub$obsPlot, na.rm = TRUE),
                    max(dfObsSub$obsPlot, na.rm = TRUE) * 1.2) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(tilde(y)[pgi](t)),
                                             expression(hat(f)[p](t)+hat(f)[pg](t)+hat(f)[pgi](t)),
                                             expression(hat(f)[p](t)+hat(f)[pg](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = title$plotTra,
                    color = "") +
      themeHDM() +
      ggplot2::facet_grid(~ genotype)
    if (!is.null(x$genoLevel$seGeno)) {
      dd1 <- dd1 +
        ggplot2::geom_ribbon(data = dfTraGenoSub,
                             ggplot2::aes(x = timeNumber,
                                          ymin = fGeno - 1.96 * seGeno,
                                          ymax = fGeno + 1.96 * seGeno,
                                          group = pop),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    if(ask) readline("Press return for next page....")
    print(dd1)
  }
}

# Needed functions for plotting --------------------------------------------------------
# Theme for the ggplot plots ----------------------------------------------

#' themeHDM
#'
#' themeHDM
#'
#' @export
themeHDM <- function(textSize = 20) {
  ggplot2::theme(strip.text.x = ggplot2::element_text(size = textSize + 5),
                 strip.text.y = ggplot2::element_text(size = textSize + 5),
                 plot.title = ggplot2::element_text(hjust = 0.5,
                                                    size = textSize + 5,
                                                    face = "bold"),
                 axis.text = ggplot2::element_text(size = textSize),
                 axis.title = ggplot2::element_text(size = textSize + 5),
                 legend.title = ggplot2::element_text(size = textSize),
                 legend.text = ggplot2::element_text(size = textSize + 5,
                                                     hjust = 0),
                 legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                 legend.key = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                 legend.justification = c(1, 1),
                 legend.position = "top",
                 panel.background = ggplot2::element_rect(fill = "gray95"))
}
