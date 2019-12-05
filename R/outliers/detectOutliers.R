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
#' detect outliers.
#' @param thrCor A numerical value used as threshold for determining outliers
#' based on correlation between plots.
#' @param thrPca A numerical value used as threshold for determining outliers
#' based on PCA scores.
#' @param title A character string, the main title added to all plots in the
#' output.
#' @param outFile A character string indicating the .csv file to which the
#' results should be written. If \code{NULL} no file is written and all plots
#' are made within R. Warning: this is potentially very slow for large numbers
#' of genotypes.
#' @param outFileOpts A named list of extra options for the pdf output file,
#' e.g. width and height. See \code{\link[grDevices]{pdf}} for all possible
#' options.
#'
#' @importFrom stats dist
#'
#' @export
detectOutliers <- function(corrDat,
                           predDat,
                           coefDat,
                           trait,
                           genotypes = unique(corrDat[["genotype"]]),
                           thrCor = 0.9,
                           thrPca = 1,
                           title = NULL,
                           outFile = NULL,
                           outFileOpts = NULL) {
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "pdf")
    ## Add outFile to output file options.
    outFileOpts <- c(list(file = outFile), outFileOpts)
    tryCatch({
      ## Open pdf.
      do.call(pdf, args = outFileOpts)
      ## Turn off device when exiting function.
      on.exit(dev.off(), add = TRUE)
    }, error = function(e) stop("Cannot open file", outFile, call. = FALSE))
  }
  ## Get corrected and predicted data per genotype.
  predDat <- predDat[predDat[["genotype"]] %in% genotypes, ]
  genoPreds <- split(x = predDat, f = predDat[["genotype"]], drop = TRUE)
  ## Restrict to corrected plots that are also in predictions.
  ## Some plots are removed while predicting the splines.
  corrDatPred <- corrDat[corrDat[["plotId"]] %in% predDat[["plotId"]], ]
  genoDats <- split(x = corrDatPred, f = corrDatPred[["genotype"]], drop = TRUE)
  ## Reshape the spline coeficients per plant x geno.
  ## First restrict to selected genotypes.
  coefDat <- coefDat[coefDat[["genotype"]] %in% genotypes, ]
  plantDats <- lapply(X = split(x = coefDat, f = coefDat[["genotype"]],
                                drop = TRUE),
                      FUN = function(dat) {
                        plantDat <- reshape2::acast(dat, formula = type ~ plotId,
                                                    value.var = "obj.coefficients")
                        ## Remove intercept.
                        return(plantDat[-1, ])
                      })
  ## Compute correlation matrix.
  cormats <- lapply(X = plantDats, FUN = function(plantDat) {
    cormat <- cor(plantDat)
    diag(cormat) <- NA
    return(cormat)
  })
  plantPcas <- lapply(X = plantDats, FUN = function(plantDat) {
    ## Perform a PCA on the spline coefficients per genotype.
    ## Run the PCA.
    plantPca <- prcomp(plantDat, center = TRUE, scale. = TRUE)
  })
  annotatePlantsCor <- lapply(X = genotypes, FUN = function(geno) {
    meanCor <- rowMeans(cormats[[geno]], na.rm = TRUE)
    if (any(meanCor < thrCor)) {
      ## Create data.frame with info on plants with average correlation
      ## below threshold.
      annPlotsCorr <- meanCor[meanCor < thrCor]
      return(data.frame(genotype = geno, plotId = names(annPlotsCorr),
                        reason = "mean corr", value = annPlotsCorr))
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
  annotatePlantsPca <- lapply(X = genotypes, FUN = function(geno) {
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
      return(data.frame(genotype = geno, plotId = names(annPlotsPc2),
                        reason = "pc2", value = annPlotsPc2))
    } else {
      return(NULL)
    }
  })
  ## Create full data.frame with annotated plants.
  annotatePlants <- do.call(rbind, c(annotatePlantsCor, annotatePlantsPca))
  ## Order by genotype and plotId.
  annotatePlants <- with(annotatePlants,
                         annotatePlants[order(genotype, plotId), ])
  ## Get minimum correlation. Cannot be higher than 0.8.
  minCor <- min(c(unlist(cormats, use.names = FALSE), 0.8), na.rm = TRUE)
  ## Create plots.
  for (geno in genotypes) {
    ## Create shape for plotting.
    ## Defaults to open circle.
    plotShapes <- setNames(rep(1, times = ncol(plantDats[[geno]])),
                           colnames(plantDats[[geno]]))
    ## Annotated plants get a closed circle.
    plotShapes[names(plotShapes) %in% annotatePlants[["plotId"]]] <- 19
    ## Plot of time course per genotype: corrected data + spline per plant.
    kinetic <- ggplot(genoDats[[geno]], aes_string(x = "timeNumber", y = trait,
                                                   color = "plotId")) +
      geom_point(aes_string(shape = "plotId"), size = 2) +
      geom_line(data = genoPreds[[geno]],
                aes_string(y = "pred.value"), size = 0.5) +
      scale_shape_manual(values = plotShapes) +
      plotTheme() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 13))
    ## Correlation plot.
    ## maybe we need to adjust the scale limits per dataset)
    correl <- ggplot(data = cormats[[geno]],
                     aes_string("Var2", "Var1", fill = "value")) +
      geom_tile(color = "white") +
      scale_fill_gradientn(colors = c("red", "white", "blue"),
                           values = scales::rescale(c(minCor, thrCor, 1)),
                           limits = c(minCor, 1),
                           name = "Pearson\nCorrelation") +
      ## Use coord fixed to create a square shaped output.
      coord_fixed() +
      theme(panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),
            panel.grid = element_line(color = "grey92"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1)) +
      labs(title = "Correl of coef", x = NULL, y = NULL)
    ## PCA biplot.
    pcaplot <- factoextra::fviz_pca_var(plantPcas[[geno]])
    ## Arrange plots.
    lay <- rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 3), c(2, 3))
    ## grid arrange always plots results.
    titleGeno <- paste("Geno", geno,
                       if (!is.null(title)) paste("-", title), "\n")
    gridExtra::grid.arrange(kinetic, correl, pcaplot,
                            layout_matrix = lay,
                            top = grid::textGrob(label = titleGeno,
                                                 gp = grid::gpar(fontsize = 15,
                                                                 fontface = 2)))
  }
  return(annotatePlants)
}



