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
#' errors are available, 95% point wise confidence intervals are depicted.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class "psHDM" as obtained after fitting
#' (\code{\link{fitSplineHDM}}) or predicting (\code{\link{predict.psHDM}}),
#' @param ... Not used.
#' @param plotType A character string indicating which plot should be made.
#' @param genotypes A character vector with the genotypes for which plots at
#' plot level are desired. Only used when \code{plotType == "genoPlotTra"}.
#' @param genotypeNames A character vector with alternative names for the
#' plotted genotypes (genotypes). If \code{NULL} the names of the genotypes
#' are used. Only used when \code{plotType == "genoPlotTra"}.
#' @param genotypeOrder A vector with the order of the selected genotypes
#' (genotypes). If \code{NULL} then the order in the data is preserved.
#' Only used when \code{plotType == "genoPlotTra"}.
#' @param xlab The x-axis label of the plot.
#' @param ylab The y-axis label of the plot.
#'
#' @examples
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#' head(spatCorrectedArch)
#' ggplot2::ggplot(data = spatCorrectedArch,
#'                ggplot2::aes(x= timeNumber, y = LeafArea_corr, group = plotId)) +
#'   ggplot2::geom_line(na.rm = TRUE) +
#'   ggplot2::facet_grid(~geno.decomp)
#'
#' ## We need to specify the genotype-by-treatment interaction
#' ## Treatment: water regime (WW, WD)
#' spatCorrectedArch$treat <- factor(spatCorrectedArch$geno.decomp,
#'                                   labels = substr(levels(spatCorrectedArch$geno.decomp), 1, 2))
#' spatCorrectedArch$genobytreat <- paste0(spatCorrectedArch$genotype, "_", spatCorrectedArch$treat)
#'
#' ## Fit P-Splines Hierarchical Curve Data Model for selection of genotypes.
#' fit.psHDM  <- fitSplineHDM(inDat = spatCorrectedArch,
#'                           trait = "LeafArea_corr",
#'                           genotypes = c("GenoA14", "GenoA51", "GenoB11",
#'                                        "GenoB02"),
#'                           time = "timeNumber",
#'                           pop = "geno.decomp",
#'                           genotype = "genobytreat",
#'                           plotId = "plotId",
#'                           difVar = list(geno = FALSE, plant = FALSE),
#'                           smoothPop = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothGeno = list(nseg = 4, bdeg = 3, pord = 2),
#'                           smoothPlot = list(nseg = 4, bdeg = 3, pord = 2),
#'                           weights = "wt",
#'                           trace = FALSE)
#'
#' ## Plot the P-Spline predictions at the three levels of the hierarchy
#'
#' ## Plots at population level.
#' plot(fit.psHDM,
#'     plotType = "popTra")
#'
#' ## Plots at genotype level.
#' plot(fit.psHDM,
#'     plotType = "popGenoTra")
#'
#' ## Plots of derivatives at genotype level.
#' plot(fit.psHDM,
#'     plotType = "popGenoDeriv")
#'
#' ## Plots of deviations at genotype level.
#' plot(fit.psHDM,
#'     plotType = "genoDev")
#'
#' ## Plots at plot level for selection of genotypes.
#' plot(fit.psHDM,
#'     plotType = "genoPlotTra",
#'     genotypes = c("GenoA14_WD", "GenoA51_WD", "GenoB11_WW",
#'                  "GenoB02_WD","GenoB02_WW"))
#'
#' @references Pérez-Valencia, D.M., Rodríguez-Álvarez, M.X., Boer, M.P. et al.
#' A two-stage approach for the spatio-temporal analysis of high-throughput
#' phenotyping data. Sci Rep 12, 3177 (2022). \doi{10.1038/s41598-022-06935-9}
#'
#' @family functions for fitting hierarchical curve data models
#'
#' @export
plot.psHDM <- function(x,
                       ...,
                       plotType = c("popTra", "popGenoTra", "popGenoDeriv",
                                    "genoDev", "genoPlotTra"),
                       genotypes = NULL,
                       genotypeNames = NULL,
                       genotypeOrder = NULL,
                       xlab = "Time",
                       ylab = expression(tilde(y)[pgi](t)),
                       title = NULL,
                       output = TRUE,
                       outFile = NULL,
                       outFileOpts = NULL) {
  if (!inherits(x, "psHDM")) {
    stop("The object class is not correct")
  }
  if (!is.null(xlab) && length(xlab) > 1) {
    stop("xlab should have length 1.\n")
  }
  if (!is.null(ylab) && length(ylab) > 1) {
    stop("ylab should have length 1.\n")
  }
  if (!is.null(title) && length(title) > 1) {
    stop("title should have length 1.\n")
  }
  plotType <- match.arg(plotType)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    on.exit(dev.off(), add = TRUE)
    do.call(pdf, args = outFileOpts)
  }
  ## Define plotting theme.
  themeHDM <-
    ggplot2::theme(strip.text.x = ggplot2::element_text(size = 25),
                   strip.text.y = ggplot2::element_text(size = 25),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 25,
                                                      face = "bold"),
                   axis.text = ggplot2::element_text(size = 20),
                   axis.title = ggplot2::element_text(size = 25),
                   legend.title = ggplot2::element_text(size = 20),
                   legend.text = ggplot2::element_text(size = 25, hjust = 0),
                   legend.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                   legend.key = ggplot2::element_rect(fill = ggplot2::alpha("white", 0)),
                   legend.justification = c(1, 1),
                   legend.position = "top",
                   panel.background = ggplot2::element_rect(fill = "gray95"))
  if (is.null(x$plotObs)) {
    x$plotObs <- x$plotLevel
  }
  minT <- min(x$plotObs[["timeNumber"]])
  maxT <- max(x$plotObs[["timeNumber"]])
  if (plotType == "popTra") {
    if (is.null(x$popLevel)) {
      stop("Population-specific growth curves can only be plotted if ",
           "predictions were made at population level.\n")
    }
    ## Population-specific growth curves.
    plotCols <- c("1" = "gray", "2" = "blue")
    p <- ggplot2::ggplot(data = x$plotObs,
                         ggplot2::aes_string(x = "timeNumber")) +
      ggplot2::geom_line(ggplot2::aes_string(y = "obsPlot", group = "plotId",
                                             color = "'1'"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$popLevel,
                         ggplot2::aes_string(y = "fPop", group = "pop",
                                             color = "'2'"), na.rm = TRUE) +
      ggplot2::geom_rug(ggplot2::aes(y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT,
                                                     length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(tilde(y)[pgi](t)),
                                             expression(hat(f)[p](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = title, color = "") +
      themeHDM +
      ggplot2::facet_grid(~ pop)
    if (!is.null(x$popLevel$sePop)) {
      p <- p +
        ggplot2::geom_ribbon(data = x$popLevel,
                             ggplot2::aes_string(ymin = "fPop - 1.96 * sePop",
                                                 ymax = "fPop + 1.96 * sePop",
                                                 group = "pop"),
                             fill = ggplot2::alpha("blue", 0.3))
    }
  } else if (plotType == "popGenoTra") {
    if (is.null(x$genoLevel) || is.null(x$popLevel)) {
      stop("Genotype-specific growth curves can only be plotted if ",
           "predictions were made at genotype and population level.\n")
    }
    ## Genotype-specific growth curves.
    plotCols <- c("1" = "blue", "2" = "orange")
    p <- ggplot2::ggplot(data = x$genoLevel,
                         ggplot2::aes_string(x = "timeNumber")) +
      ggplot2::geom_line(ggplot2::aes_string(y = "fGeno", group = "genotype",
                                             color = "'1'"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$popLevel,
                         ggplot2::aes_string(y = "fPop", group = "pop",
                                             color = "'2'"),
                         size = 0.8, na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plotObs,
                        ggplot2::aes(y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT,
                                                     length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(hat(f)[p](t)+hat(f)[pg](t)),
                                             expression(hat(f)[p](t)))) +
      ggplot2::labs(x = xlab, y = ylab, color = "", title = title) +
      themeHDM +
      ggplot2::facet_grid(~ pop)
    if (!is.null(x$popLevel$sePop)) {
      p <- p +
        ggplot2::geom_ribbon(data = x$popLevel,
                             ggplot2::aes_string(ymin = "fPop - 1.96 * sePop",
                                                 ymax = "fPop + 1.96 * sePop",
                                                 group = "pop"),
                             fill = ggplot2::alpha("orange", 0.5))
    }
  } else if (plotType == "popGenoDeriv") {
    if (is.null(x$genoLevel)) {
      stop("First order derivatives of genotype-specific growth curves can ",
           "only be plotted if predictions were made at genotype level.\n")
    }
    ## First derivative of the genotype-specific growth curves.
    plotCols <- c("1" = "blue", "2" = "orange")
    p <- ggplot2::ggplot(data = x$genoLevel,
                         ggplot2::aes_string(x = "timeNumber")) +
      ggplot2::geom_line(ggplot2::aes_string(y = "fGenoDeriv1", group = "genotype",
                                             color = "'1'"), na.rm = TRUE) +
      ggplot2::geom_line(data = x$popLevel,
                         ggplot2::aes_string(y = "fPopDeriv1", group = "pop",
                                             color = "'2'"), na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plotObs,
                        ggplot2::aes(y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT,
                                                     length.out = 5), 0)) +
      ggplot2::ylim(min(x$genoLevel$fGenoDeriv1),
                    max(x$genoLevel$fGenoDeriv1) * 1.2) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression((hat(f)[p](t)+hat(f)[pg](t))*minute),
                                             expression((hat(f)[p](t))*minute))) +
      ggplot2::labs(x = xlab, y = "First derivative", color = "",
                    title = title) +
      themeHDM +
      ggplot2::facet_grid(~ pop)
    if (!is.null(x$popLevel$sePopDeriv1)) {
      p <- p +
        ggplot2::geom_ribbon(data = x$popLevel,
                             ggplot2::aes_string(
                               ymin = "fPopDeriv1 - 1.96 * sePopDeriv1",
                               ymax = "fPopDeriv1 + 1.96 * sePopDeriv1",
                               group = "pop"),
                             fill = ggplot2::alpha("orange", 0.5))
    }
  } else if (plotType == "genoDev") {
    if (is.null(x$genoLevel)) {
      stop("Genotype-specific deviations can only be plotted if ",
           "predictions were made at genotype level.\n")
    }
    ## Estimated genotypic deviations.
    plotCols <- c("1" = "blue")
    p <- ggplot2::ggplot(data = x$genoLevel,
                         ggplot2::aes_string(x = "timeNumber")) +
      ggplot2::geom_line(ggplot2::aes_string(y = "fGenoDev", group = "genotype",
                                             color = "'1'"), na.rm = TRUE) +
      ggplot2::geom_rug(data = x$plotObs,
                        ggplot2::aes(y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT,
                                                     length.out = 5), 0)) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(hat(f)[pg](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = title,
                    color = "") +
      themeHDM +
      ggplot2::facet_grid(~ pop)
  } else if (plotType == "genoPlotTra") {
    if (is.null(x$genoLevel) || is.null(x$plotLevel)) {
      stop("Plot and Genotype-specific growth curves can only be plotted if ",
           "predictions were made at genotype and plot level.\n")
    }
    if (!is.null(genotypes)) {
      if (!all(genotypes %in% x$plotLevel[["genotype"]])) {
        stop("All genotypes should be in fitted model.\n")
      }
    } else {
      genotypes <- levels(x$plotLevel[["genotype"]])
    }
    ## Estimated plot trajectories.
    dfTraSub <- droplevels(x$plotLevel[x$plotLevel[["genotype"]] %in% genotypes, ])
    dfObsSub <- droplevels(x$plotObs[x$plotObs[["genotype"]] %in% genotypes, ])
    dfTraGenoSub <- droplevels(x$genoLevel[x$genoLevel[["genotype"]] %in% genotypes, ])
    if (!is.null(genotypeNames)) {
      dfTraSub[["genotype"]] <- factor(dfTraSub[["genotype"]],
                                       labels = genotypeNames)
      dfObsSub[["genotype"]] <- factor(dfObsSub[["genotype"]],
                                       labels = genotypeNames)
      dfTraGenoSub[["genotype"]] <- factor(dfTraGenoSub[["genotype"]],
                                           labels = genotypeNames)
    }
    if (!is.null(genotypeOrder)) {
      dfTraSub[["genotype"]] <- factor(dfTraSub[["genotype"]],
                                       levels = levels(dfTraSub[["genotype"]])[genotypeOrder])
      dfObsSub[["genotype"]] <- factor(dfObsSub[["genotype"]],
                                       levels = levels(dfObsSub[["genotype"]])[genotypeOrder])
      dfTraGenoSub[["genotype"]] <- factor(dfTraGenoSub[["genotype"]],
                                           levels = levels(dfTraGenoSub[["genotype"]])[genotypeOrder])
    }
    plotCols <- c("1" = "gray", "2" = "blue", "3" = "red")
    p <- ggplot2::ggplot(data = dfObsSub,
                         ggplot2::aes_string(x = "timeNumber")) +
      ggplot2::geom_line(ggplot2::aes_string(y = "obsPlot", group = "plotId",
                                             color = "'1'"), na.rm = TRUE) +
      ggplot2::geom_line(data = dfTraSub,
                         ggplot2::aes_string(y = "fPlot", group = "plotId",
                                             color = "'2'"),
                         linetype = 2, na.rm = TRUE) +
      ggplot2::geom_line(data = dfTraGenoSub,
                         ggplot2::aes_string(y = "fGeno", group = "genotype",
                                             color = "'3'"), na.rm = TRUE) +
      ggplot2::geom_rug(ggplot2::aes(y = NULL), color = "gray",
                        length = ggplot2::unit(0.01, "npc")) +
      ggplot2::scale_x_continuous(breaks = round(seq(minT, maxT,
                                                     length.out = 5), 0)) +
      ggplot2::ylim(min(dfObsSub$obsPlot, na.rm = TRUE),
                    max(dfObsSub$obsPlot, na.rm = TRUE) * 1.2) +
      ggplot2::scale_color_manual(values = plotCols,
                                  labels = c(expression(tilde(y)[pgi](t)),
                                             expression(hat(f)[p](t)+hat(f)[pg](t)+hat(f)[pgi](t)),
                                             expression(hat(f)[p](t)+hat(f)[pg](t)))) +
      ggplot2::labs(x = xlab, y = ylab, title = title,
                    color = "") +
      themeHDM
    if (!is.null(x$genoLevel$seGeno)) {
      p <- p +
        ggplot2::geom_ribbon(data = dfTraGenoSub,
                             ggplot2::aes_string(
                               ymin = "fGeno - 1.96 * seGeno",
                               ymax = "fGeno + 1.96 * seGeno",
                               group = "pop"),
                             fill = ggplot2::alpha("orange", 0.5))
    }
    ## Calculate the total number of plots.
    nPlots <- nlevels(dfTraSub[["genotype"]])
    ## 25 Plots per page.
    nPag <- ceiling(nPlots / 25)
    if (nPlots >= 25) {
      ## More than 25 plots.
      ## For identical layout on all pages use 5 x 5 plots throughout.
      rowPag <- colPag <- rep(x = 5, times = nPag)
      # 28-7-2020. ggforce has a bug that prevents this identical layout
      # https://github.com/thomasp85/ggforce/issues/201
      # When fixed the code above can be reactivated and the three lines below
      # removed.
      # plotsLastPag <- nPlots %% 25
      # rowPag <- c(rep(x = 5, times = nPag - 1), min(plotsLastPag %/% 5 + 1, 5))
      # colPag <- c(rep(x = 5, times = nPag - 1),
      #             ifelse(plotsLastPag >= 5, 5, plotsLastPag))
    } else {
      ## Less than 25 plots.
      ## Fill page by row of 5 plots.
      plotsPag <- nPlots %% 25
      rowPag <- min(ceiling(plotsPag / 5), 5)
      colPag <- ifelse(plotsPag >= 5, 5, plotsPag)
    }
    ## Build pages of plots.
    pPag <- vector(mode = "list", length = nPag)
    for (i in 1:nPag) {
      pPag[[i]] <- p +
        ggforce::facet_wrap_paginate(facets = "genotype", nrow = rowPag[i],
                                     ncol = colPag[i],
                                     labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
                                     page = i)
      if (output) {
        suppressMessages(plot(pPag[[i]]))
      }
    }
    return(invisible(pPag))
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}

# Needed functions for plotting --------------------------------------------------------
# Theme for the ggplot plots ----------------------------------------------

#' themeHDM
#'
#' Custom theme for creating HDM plots.
#'
#' @param textSize Numeric value, font size of the texts in the plots.
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
