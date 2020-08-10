#' Outlier detection
#'
#' Function for detecting outliers.
#'
#' @param corrDat A data.frame with corrected spatial data.
#' @param predDat A data.frame with predicted data from a fitted spline.
#' @param coefDat A data.frame with coefficients from a fitted spline.
#' @param trait A character string indicating the trait for which to detect
#' outliers.
#' @param genotypes A character vector indicating the genotypes for which to
#' detect outliers. If \code{NULL}, outlier detection will be done for all
#' genotypes.
#' @param geno.decomp A character vector indicating the variables to use to
#' group the genotypic variance in the model.
#' @param thrCor A numerical value used as threshold for determining outliers
#' based on correlation between plots.
#' @param thrPca A numerical value used as threshold for determining outliers
#' based on PCA scores.
#'
#' @importFrom stats dist
#'
#' @export
detectTimeCourseOutliers <- function(corrDat,
                                     predDat,
                                     coefDat,
                                     trait,
                                     genotypes = NULL,
                                     geno.decomp = NULL,
                                     thrCor = 0.9,
                                     thrPca = 1) {
  ## Checks.
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length 1.\n")
  }
  if (!inherits(corrDat, "data.frame")) {
    stop("corrDat should be a data.frame.\n")
  }
  corrCols <- c("plotId", "genotype", trait, geno.decomp)
  if (!all(hasName(x = corrDat, name = corrCols))) {
    stop("corrDat should at least contain the following columns: ",
         paste(corrCols, collapse = ", "))
  }
  if (!inherits(predDat, "data.frame")) {
    stop("predDat should be a data.frame.\n")
  }
  predCols <- c("plotId", "genotype", "pred.value")
  if (!all(hasName(x = predDat, name = predCols))) {
    stop("predDat should at least contain the following columns: ",
         paste(predCols, collapse = ", "))
  }
  if (!inherits(coefDat, "data.frame")) {
    stop("coefDat should be a data.frame.\n")
  }
  coefCols <- c("plotId", "genotype", "type", "obj.coefficients")
  if (!all(hasName(x = coefDat, name = coefCols))) {
    stop("coefDat should at least contain the following columns: ",
         paste(coefCols, collapse = ", "))
  }
  if (is.null(genotypes)) {
    genotypes = unique(as.character(predDat[["genotype"]]))
  } else {
    if (!is.character(genotypes)) {
      stop("genotypes should be a character vector.\n")
    }
    if (!all(genotypes %in% predDat[["genotype"]])) {
      stop("all genotypes should be in predDat")
    }
  }
  if (!all(genotypes %in% coefDat[["genotype"]])) {
    stop("all genotypes should be in coefDat")
  }
  if (!all(genotypes %in% corrDat[["genotype"]])) {
    stop("all genotypes should be in corrDat")
  }
  if (!is.numeric(thrCor) || length(thrCor) > 1 || thrCor < 0 || thrCor > 1) {
    stop("thrCor should be a numerical value between 0 and 1.\n")
  }
  if (!is.numeric(thrPca) || length(thrPca) > 1 || thrPca < 0) {
    stop("thrPca should be a positive numerical value.\n")
  }
  ## Restrict corrDat, predDat and coefDat to genotypes.
  corrDat <- corrDat[corrDat[["genotype"]] %in% genotypes, ]
  ## Get corrected and predicted data per genotype.
  ## Merge geno.decomp to coefDat.
  if (!is.null(geno.decomp)) {
    predDat <- merge(predDat, unique(corrDat[c("plotId", geno.decomp)]))
  }
  genoPreds <- split(x = predDat,
                     f = predDat[c("genotype", geno.decomp)], drop = TRUE)
  ## Restrict to corrected plots that are also in predictions.
  ## Some plots are removed while predicting the splines.
  corrDatPred <- corrDat[corrDat[["plotId"]] %in% predDat[["plotId"]], ]
  genoDats <- split(x = corrDatPred,
                    f = corrDatPred[c("genotype", geno.decomp)], drop = TRUE)
  ## Reshape the spline coeficients per plant x geno.
  ## First restrict to selected genotypes.
  coefDat <- coefDat[coefDat[["genotype"]] %in% genotypes, ]
  ## Merge geno.decomp to coefDat.
  if (!is.null(geno.decomp)) {
    coefDat <- merge(coefDat, unique(corrDatPred[c("plotId", geno.decomp)]))
  }
  plantDats <- lapply(X = split(x = coefDat,
                                f = coefDat[c("genotype", geno.decomp)],
                                drop = TRUE),
                      FUN = function(dat) {
                        plantDat <- reshape2::acast(dat, formula = type ~ plotId,
                                                    value.var = "obj.coefficients")
                        ## Remove intercept.
                        return(plantDat[-1, ])
                      })
  ## Compute correlation matrix.
  cormats <- lapply(X = plantDats, FUN = function(plantDat) {
    if (!is.null(dim(plantDat))) {
      ## if there are plants, estimate the correlation...
      cormat <- cor(plantDat)
      diag(cormat) <- NA
      return(cormat)
    } else {
      ## ... if not, return a null matrix
      return(NULL)
    }
  })
  plantPcas <- lapply(X = plantDats, FUN = function(plantDat) {
    ## Perform a PCA on the spline coefficients per genotype.
    ## Run the PCA.
    if (!is.null(dim(plantDat))) {
      ## if there are plants, perform the PCA...
      return(prcomp(plantDat, center = TRUE, scale. = TRUE))
    } else {
      ## ... if not, return a null object
      return(NULL)
    }
  })
  annotatePlantsCor <- lapply(X = names(cormats), FUN = function(geno) {
    if (!is.null(cormats[[geno]])) {
      meanCor <- rowMeans(cormats[[geno]], na.rm = TRUE)
    } else {
      meanCor <- NULL
    }
    if (any(meanCor < thrCor)) {
      ## Create data.frame with info on plants with average correlation
      ## below threshold.
      annPlotsCorr <- meanCor[meanCor < thrCor]
      return(data.frame(plotId = names(annPlotsCorr), reason = "mean corr",
                        value = annPlotsCorr))
    } else {
      return(NULL)
    }
  })
  cormats <- lapply(X = cormats, FUN = function(cormat) {
    ## Set lower part of cormat to NA for plotting.
    cormat[lower.tri(cormat)] <- NA
    ## Melt to format used by ggplot.
    meltedCormat <- reshape2::melt(cormat, na.rm = TRUE)
  })
  annotatePlantsPca <- lapply(X = names(plantPcas), FUN = function(geno) {
    ## Calculate the pairwise difference of coordinates on the 2nd axis and
    ## annotate plant with average diff larger than threshold.
    diffPc2 <- as.matrix(dist(plantPcas[[geno]]$rotation[, "PC2"]))
    diag(diffPc2) <- NA
    meanPc2 <- rowMeans(diffPc2, na.rm = TRUE)
    ## Create shape for plotting.
    if (any(meanPc2 >= thrPca)) {
      ## Create data.frame with info on plants with avarage difference.
      ## above threshold.
      annPlotsPc2 <- meanPc2[meanPc2 >= thrPca]
      return(data.frame(plotId = names(annPlotsPc2), reason = "pc2",
                        value = annPlotsPc2))
    } else {
      return(NULL)
    }
  })
  ## Create full data.frame with annotated plants.
  annotatePlants <- do.call(rbind, c(annotatePlantsCor, annotatePlantsPca))
  if (!is.null(annotatePlants)) {
    ## Merge genotype and geno.decomp to annotated plants.
    annotatePlants <- merge(unique(corrDatPred[c("genotype", geno.decomp, "plotId")]),
                            annotatePlants)
    ## Order by genotype, geno.decomp and plotId.
    if (!is.null(geno.decomp)) {
      annOrd <- order(annotatePlants[["genotype"]],
                      annotatePlants[[geno.decomp]], annotatePlants[["plotId"]])
    } else {
      annOrd <- order(annotatePlants[["genotype"]], annotatePlants[["plotId"]])
    }
    annotatePlants <- annotatePlants[annOrd, ]
  } else {
    annotatePlants <- data.frame()
  }
  plotInfo <- unique(corrDatPred[c("genotype", geno.decomp)])
  class(annotatePlants) <- c("timeCourseOutliers", class(annotatePlants))
  attr(x = annotatePlants, which = "thrCor") <- thrCor
  attr(x = annotatePlants, which = "thrPca") <- thrPca
  attr(x = annotatePlants, which = "trait") <- trait
  attr(x = annotatePlants, which = "geno.decomp") <- geno.decomp
  attr(x = annotatePlants, which = "plotInfo") <- plotInfo
  attr(x = annotatePlants, which = "cormats") <- cormats
  attr(x = annotatePlants, which = "plantPcas") <- plantPcas
  attr(x = annotatePlants, which = "genoPreds") <- genoPreds
  attr(x = annotatePlants, which = "genoDats") <- genoDats
  return(annotatePlants)
}

#' plot.timeCourseOutliers
#'
#' @inheritParams detectTimeCourseOutliers
#' @inheritParams plot.TP
#'
#' @param x An object of class timeCourseOutliers.
#'
#' @export
plot.timeCourseOutliers <- function(x,
                                    ...,
                                    genotypes = NULL,
                                    title = NULL,
                                    output = TRUE) {
  thrCor <- attr(x = x, which = "thrCor")
  thrPca <- attr(x = x, which = "thrPca")
  trait <- attr(x = x, which = "trait")
  geno.decomp <- attr(x = x, which = "geno.decomp")
  plotInfo <- attr(x = x, which = "plotInfo")
  if (!is.null(genotypes) && (!is.character(genotypes) ||
                              !all(genotypes %in% plotInfo[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes used for ",
         "outlier detection.\n")
  }
  if (is.null(genotypes)) {
    genotypes <- interaction(plotInfo[c("genotype", geno.decomp)])
  } else {
    genotypes <- interaction(plotInfo[plotInfo[["genotype"]] %in% genotypes,
                                      c("genotype", geno.decomp)])
  }
  genotypes <- as.character(genotypes)
  cormats <- cormats <- attr(x = x, which = "cormats")[genotypes]
  plantPcas <- attr(x = x, which = "plantPcas")[genotypes]
  genoPreds <- attr(x = x, which = "genoPreds")[genotypes]
  genoDats <- attr(x = x, which = "genoDats")[genotypes]
  ## Get minimum correlation. Cannot be higher than 0.8.
  minCor <- min(c(unlist(cormats, use.names = FALSE), 0.8), na.rm = TRUE)
  ## Compute the number of breaks for the time scale based on all plants.
  ## If there are less than 3 time points use the number of time points.
  ## Otherwise use 3.
  useTimePoint <- hasName(x = genoPreds[[1]], name = "timePoint")
  timeVar <- if (useTimePoint) "timePoint" else "timeNumber"
  ## Create plots.
  if (!output) {
    ## When calling arrangeGrob a blank page is opened if no plotting
    ## device is currently open.
    ## Adding a call to an empty pdf file prevents this empty plot from
    ## being opened.
    pdf(file = NULL)
    on.exit(dev.off(), add = TRUE)
  }
  p <- lapply(X = genotypes, FUN = function(genotype) {
    ## Create shape for plotting.
    ## Defaults to open circle.
    plotIds <- unique(genoPreds[[genotype]][["plotId"]])
    plotShapes <- setNames(rep(1, times = length(plotIds)), plotIds)
    ## Annotated plants get a closed circle.
    plotShapes[names(plotShapes) %in% x[["plotId"]]] <- 19
    ## Plot of time course per genotype: corrected data + spline per plant.
    kinetic <- ggplot(genoDats[[genotype]], aes_string(x = timeVar, y = trait,
                                                       color = "plotId")) +
      geom_point(aes_string(shape = "plotId"), size = 2, na.rm = TRUE) +
      geom_line(data = genoPreds[[genotype]],
                aes_string(y = "pred.value"), size = 0.5) +
      scale_shape_manual(values = plotShapes) +
      theme_light() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 13))
    if (useTimePoint) {
      ## Format the time scale to Month + day.
      kinetic <- kinetic +
        ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                  labels = scales::date_format("%B %d"))
    }
    ## Correlation plot.
    correl <- ggplot(data = cormats[[genotype]],
                     aes_string("Var2", "Var1", fill = "value")) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colors = c("red", "white", "blue"),
                           values = scales::rescale(c(minCor, thrCor, 1)),
                           limits = c(minCor, 1),
                           name = "Pearson\nCorrelation") +
      ## Move y-axis to the right.
      scale_y_discrete(position = "right") +
      ## Use coord fixed to create a square shaped output.
      coord_fixed() +
      theme(panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),
            panel.grid = element_line(color = "grey92"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1),
            legend.position = "left") +
      labs(title = "Correl of coef", x = NULL, y = NULL)
    ## PCA biplot.
    pcaplot <- factoextra::fviz_pca_var(plantPcas[[genotype]])
    ## Arrange plots.
    lay <- rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 3), c(2, 3))
    ## grid arrange always plots results.
    titleGeno <- paste("Geno", genotype,
                       if (!is.null(title)) paste("-", title), "\n")
    if (genotype == genotypes[1] && output) {
      ## Arrange grob always needs an open device.
      ## This creates a blank first page when first plotting.
      ## By opening a new page for the first plot and then using
      ## newpage = FALSE in the actual plot this blank page is overwritten.
      grid::grid.newpage()
    }
    pGeno <- gridExtra::arrangeGrob(kinetic, correl, pcaplot,
                                    layout_matrix = lay,
                                    top = grid::textGrob(label = titleGeno,
                                                         gp = grid::gpar(fontsize = 15,
                                                                         fontface = 2)))
    if (output) {
      gridExtra::grid.arrange(pGeno, newpage = (genotype != genotypes[1]))
    }
    return(pGeno)
  })
  invisible(p)
}



