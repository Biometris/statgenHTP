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
                           outFile = NULL,
                           outFileOpts = NULL) {
  if (!is.null(outFile)) {
    ## Check if file exists and is writable.
    chkFile(outFile, fileType = "pdf")
    ## Add outFile to output file options.
    outFileOpts <- c(list(file = outFile), outFileOpts)
    ## Turn off device when exiting function.
    on.exit(dev.off(), add = TRUE)
    ## Open pdf.
    do.call(pdf, args = outFileOpts)
  }
  ## Prepare table for annotated plants.
  annotatePlant <- data.frame()
  for (geno in genotypes) {
    ## Get corrected and predicted data for current genotype.
    genoPred <- predDat[predDat[["genotype"]] == geno, ]
    genoDat <- corrDat[corrDat[["genotype"]] == geno &
                         corrDat[["plotId"]] %in% genoPred[["plotId"]], ]
    ## Reshape the spline coeficients per plant x geno.
    plantDat <- reshape2::acast(coefDat[coefDat[["genotype"]] == geno, ],
                                type ~ plotId, value.var = "obj.coefficients")
    ## Compute correlation matrix, exclude intercept.
    cormat <- round(cor(plantDat[-1, ]), 2)
    diag(cormat) <- NA
    ## Annotate plants with average correl < 0.9.
    ## (needs manual adjustment per dataset)
    thrCor <- 0.9
    ## Get avarage correlation per plant.
    meanCor <- rowMeans(cormat, na.rm = TRUE)
    if (any(meanCor < thrCor)) {
      ## Create data.frame with info on plants with average correlation
      ## below threshold.
      annPlotsCorr <- meanCor[meanCor < thrCor]
      annotateCorr <- data.frame(genotype = geno, plotId = names(annPlotsCorr),
                                 reason = "mean corr", value = annPlotsCorr)
      annotatePlant <- rbind(annotatePlant, annotateCorr)
    }
    ## Perform a PCA on the spline coefficients per genotype.
    ## Run the PCA.
    plantPca <- prcomp(plantDat[-1, ], center = TRUE, scale. = TRUE)
    ## Calculate the pairwise difference of coordinates on the 2nd axis and
    ## annotate plant with average diff > 1.
    ## (manual adjustment of the threshold)
    thrPca <- 1
    diffPc2 <- as.matrix(dist(plantPca$rotation[, "PC2"]))
    diag(diffPc2) <- NA
    meanPc2 <- rowMeans(diffPc2, na.rm = TRUE)
    ## Create shape for plotting.
    plotShapes <- setNames(rep(1, times = nrow(cormat)), rownames(cormat))
    if (any(meanPc2 >= thrPca)) {
      ## Create data.frame with info on plants with avarage difference.
      ## above threshold.
      annPlotsPc2 <- meanPc2[meanPc2 >= thrPca]
      annotatePc2 <- data.frame(genotype = geno, plotId = names(annPlotsPc2),
                                reason = "pc2", value = annPlotsPc2)
      ## Annotated plants get a filled circle.
      plotShapes[names(plotShapes) %in% names(annPlotsPc2)] <- 19
      annotatePlant <- rbind(annotatePlant, annotatePc2)
    }
    ## Plot of time course per genotype: corrected data + spline per plant.
    kinetic <- ggplot(genoDat, aes_string(x = "timeNumber", y = trait,
                                          color = "plotId")) +
      geom_point(aes_string(shape = "plotId"), size = 2) +
      geom_line(data = genoPred, aes_string(y = "pred.value"), size = 0.5) +
      scale_shape_manual(values = plotShapes) +
      plotTheme() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 13))
    ## Correlation plot.
    ## maybe we need to adjust the scale limits per dataset)
    ## Set lower part of cormat to NA for plotting.
    cormat[lower.tri(cormat)] <- NA
    ## Melt to format used by ggplot.
    meltedCormat <- reshape2::melt(cormat, na.rm = TRUE)
    correl <- ggplot(data = meltedCormat,
                     aes_string("Var2", "Var1", fill = "value")) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                           midpoint = 0.9, limit = c(0.8, 1), space = "Lab",
                           name = "Pearson\nCorrelation") +
      coord_fixed() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1)) +
      labs(title = "Correl of coef", x = NULL, y = NULL)
    ## PCA biplot.
    pcaplot <- ggbiplot::ggbiplot(plantPca, obs.scale = 1) +
      plotTheme() +
      theme(aspect.ratio = 1) +
      ggtitle("PCA of coef")
    ## Arrange plots.
    lay <- rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 3), c(2, 3))
    ## grid arrange always plots results.
    gridExtra::grid.arrange(kinetic, correl, pcaplot,
                            layout_matrix = lay,
                            top = grid::textGrob(paste0("Geno ", geno, " - Phenvator Rene - PAM\n"),
                                                 gp = grid::gpar(fontsize = 15, fontface = 2)))
  }
  return(annotatePlant)
}



