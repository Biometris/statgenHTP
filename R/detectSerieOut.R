#' Detect outliers for series of observations
#'
#' Function for detecting strange time courses. The function uses the estimates
#' for the spline coefficients per time course (typically per plant).
#' Correlations between those coefficient vectors are calculated to identify
#' outlying time courses, i.e., plants. An outlying time course will have low
#' correlation to the majority of time courses. To support the analysis by
#' correlations, a principal component analysis is done on the plant
#' (time course) by spline coefficient matrix. A PCA plot of the plant scores
#' will show the outlying plants. Finally the pairwise-ratios of the slopes of
#' a linear model fitted through the spline coefficients are computed. Plants
#' are tagged when the average pairwise-ratio is lower the a given threshold
#' (`thrSlope`).
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
#' based on angles (in degrees) between PCA scores.
#' @param thrSlope A numerical value used as threshold for determining outliers
#' based on slopes.
#'
#' @returns An object of class \code{serieOut}, a \code{data.frame} with outlying
#' series of observations.
#'
#' @examples
#' \donttest{
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## Fit P-splines on a subset of genotypes
#' subGenoVator <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGenoVator,
#'                         knots = 50)
#'
#' ## Extract the data.frames with predicted values and P-Spline coefficients.
#' predDat <- fit.spline$predDat
#' coefDat <- fit.spline$coefDat
#'
#' ## The coefficients are then used to tag suspect time courses.
#' outVator <- detectSerieOut(corrDat = spatCorrectedVator,
#'                            predDat = predDat,
#'                            coefDat = coefDat,
#'                            trait = "EffpsII_corr",
#'                            genotypes = subGenoVator,
#'                            thrCor = 0.9,
#'                            thrPca = 30,
#'                            thrSlope = 0.7)
#'
#' ## The `outVator` can be visualized for selected genotypes.
#' plot(outVator, genotypes = "G151")
#' }
#'
#' @family functions for detecting outliers for series of observations
#'
#' @importFrom utils combn
#' @export
detectSerieOut <- function(corrDat,
                           predDat,
                           coefDat,
                           trait,
                           genotypes = NULL,
                           geno.decomp = NULL,
                           thrCor = 0.9,
                           thrPca = 30,
                           thrSlope = 0.7) {
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
    genotypes <- unique(as.character(predDat[["genotype"]]))
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
  if (!is.numeric(thrCor) || any(thrCor < -1) || any(thrCor > 1)) {
    stop("thrCor should be a numerical vector with values between -1 and 1.\n")
  }
  if (!is.numeric(thrPca) || any(thrPca < 0)) {
    stop("thrPca should be a numerical vector with positive values.\n")
  }
  if (!is.numeric(thrSlope) || any(thrSlope < 0) || any(thrSlope > 1)) {
    stop("thrSlope should be a numerical vector with values between 0 and 1.\n")
  }
  ## Restrict corrDat to genotypes.
  corrDat <- corrDat[corrDat[["genotype"]] %in% genotypes, ]
  ## Get number of values for geno.decomp.
  if (!is.null(geno.decomp)) {
    genoDecompLevs <- unique(corrDat[[geno.decomp]])
    nGenoDecomp <- length(genoDecompLevs)
    if (length(thrCor) == 1) {
      thrCor <- setNames(rep(thrCor, times = nGenoDecomp), genoDecompLevs)
    } else if (is.null(names(thrCor)) ||
               !all(genoDecompLevs %in% names(thrCor))) {
      stop("thrCor should be a named vector, with names matching the levels ",
           "in geno.decomp.\n")
    }
    if (length(thrPca) == 1) {
      thrPca <- setNames(rep(thrPca, times = nGenoDecomp), genoDecompLevs)
    } else if (is.null(names(thrPca)) ||
               !all(genoDecompLevs %in% names(thrPca))) {
      stop("thrPca should be a named vector, with names matching the levels ",
           "in geno.decomp.\n")
    }
    if (length(thrSlope) == 1) {
      thrSlope <- setNames(rep(thrSlope, times = nGenoDecomp), genoDecompLevs)
    } else if (is.null(names(thrSlope)) ||
               !all(genoDecompLevs %in% names(thrSlope))) {
      stop("thrSlope should be a named vector, with names matching the levels ",
           "in geno.decomp.\n")
    }
  } else {
    if (length(thrCor) != 1) {
      stop("thrCor should be a vector of length 1.\n")
    }
    if (length(thrPca) != 1) {
      stop("thrPca should be a vector of length 1.\n")
    }
    if (length(thrSlope) != 1) {
      stop("thrSlope should be a vector of length 1.\n")
    }
  }
  ## Restrict predDat and coefDat to genotypes.
  ## Get corrected and predicted data per genotype.
  ## Merge geno.decomp to coefDat.
  if (!is.null(geno.decomp)) {
    predDat <- merge(predDat, unique(corrDat[c("plotId", geno.decomp)]))
  }
  ## Restrict to corrected plots that are also in predictions.
  ## Some plots are removed while predicting the splines.
  corrDatPred <- corrDat[corrDat[["plotId"]] %in% predDat[["plotId"]], ]
  NAplants <- tapply(corrDatPred[[trait]], droplevels(corrDatPred[["plotId"]]),
                     FUN = function(x) {
                       all(is.na(x))
                     })
  corrDatPred <- corrDatPred[corrDatPred[["plotId"]] %in%
                               names(NAplants)[!NAplants], ]
  genoDats <- split(x = corrDatPred,
                    f = corrDatPred[c("genotype", geno.decomp)], drop = TRUE)
  predDat <- predDat[predDat[["plotId"]] %in%
                       names(NAplants)[!NAplants], ]
  genoPreds <- split(x = predDat,
                     f = predDat[c("genotype", geno.decomp)], drop = TRUE)
  ## Reshape the spline coefficients per plant x geno.
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
                        plantDat <- reshape2::acast(dat,
                                                    formula = type ~ plotId,
                                                    value.var = "obj.coefficients")
                        ## Remove intercept.
                        #plantDat <- plantDat[-1, , drop = FALSE]
                        ## Remove plants with only NA.
                        NAplants <- apply(X = plantDat, MARGIN = 2,
                                          FUN = function(plant) {
                                            all(is.na(plant))
                                          })
                        plantDat <- plantDat[, !NAplants, drop = FALSE]
                        ## Add geno.decomp as attribute for later use.
                        if (!is.null(geno.decomp)) {
                          attr(plantDat, which = "genoDecomp") <-
                            as.character(unique(dat[[geno.decomp]]))
                        }
                        return(plantDat)
                      })
  genoPlotId <- sapply(X = plantDats, FUN = ncol)
  genoPlotIdLim <- names(genoPlotId[genoPlotId < 3])
  if (length(genoPlotIdLim) > 0) {
    warning("The following genotypes have less than 3 plotIds and are skipped ",
            "in the outlier detection:\n",
            paste(genoPlotIdLim, collapse = ", "), "\n", call. = FALSE)
    plantDats[genoPlotIdLim] <- NULL
    if (length(plantDats) == 0) {
      stop("No genotypes left for performing outlier detection.\n")
    }
  }
  ## Compute correlation matrices.
  cormats <- lapply(X = plantDats, FUN = function(plantDat) {
    if (!is.null(dim(plantDat))) {
      ## if there are plants, estimate the correlation...
      cormat <- cor(plantDat)
      diag(cormat) <- NA
      attr(cormat, which = "genoDecomp") <- attr(plantDat, which = "genoDecomp")
      if (!is.null(geno.decomp)) {
        attr(cormat, which = "thrCor") <-
          thrCor[attr(plantDat, which = "genoDecomp")]
      } else {
        attr(cormat, which = "thrCor") <- thrCor
      }
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
      plantPca <- prcomp(plantDat, center = TRUE, scale. = TRUE)
      attr(plantPca, which = "genoDecomp") <- attr(plantDat, which = "genoDecomp")
      ## if there are plants, perform the PCA...
      return(plantPca)
    } else {
      ## ... if not, return a null object
      return(NULL)
    }
  })
  ## Compute slope matrices.
  slopemats <- lapply(X = plantDats, FUN = function(plantDat) {
    if (!is.null(dim(plantDat))) {
      ## if there are plants, estimate the slopes.
      ## Convert matrix to data.frame for lm.
      plantDatLm <- as.data.frame(plantDat)
      ## Construct empty matrix for storing results.
      slopemat <- matrix(nrow = ncol(plantDatLm), ncol = ncol(plantDatLm),
                         dimnames = list(colnames(plantDatLm),
                                         colnames(plantDatLm)))
      ## Compute slope per pair of plots.
      slopemat[lower.tri(slopemat)] <-
        combn(x = colnames(plantDatLm), m = 2, FUN = function(plants) {
          ## Wrap plants in ` to allow for irregular names in formula.
          plants <- paste0("`", plants, "`")
          ## Fit linear model and extract slope.
          modForm <- formula(paste(plants, collapse = "~"))
          slope <- abs(coef(lm(modForm, data = plantDatLm))[2])
          if (slope > 1) slope <- 1 / slope
          return(slope)
        }, simplify = TRUE)
      ## Fill upper triangle of slopemat matrix with copy of lower triangle.
      slopemat[upper.tri(slopemat)] <- t(slopemat)[upper.tri(slopemat)]
      attr(slopemat, which = "genoDecomp") <-
        attr(plantDat, which = "genoDecomp")
      if (!is.null(geno.decomp)) {
        attr(slopemat, which = "thrSlope") <-
          thrSlope[attr(plantDat, which = "genoDecomp")]
      } else {
        attr(slopemat, which = "thrSlope") <- thrSlope
      }
      return(slopemat)
    } else {
      ## ... if not, return a null matrix
      return(NULL)
    }
  })
  ## Find outlying plots based on correlation.
  annotatePlantsCor <- lapply(X = names(cormats), FUN = function(geno) {
    if (!is.null(cormats[[geno]])) {
      ## Compute mean based on all but the worst correlation per plant.
      ## This prevents one bad correlation causing all plants to be outliers.
      meanCor <- apply(cormats[[geno]], MARGIN = 1, FUN = function(plant) {
        mean(plant[-which(rank(plant) == 1)], na.rm = TRUE)
      })
      thrCorPlant <- attr(cormats[[geno]], which = "thrCor")
    } else {
      meanCor <- NULL
    }
    if (any(meanCor < thrCorPlant)) {
      ## Create data.frame with info on plants with average correlation
      ## below threshold.
      annPlotsCorr <- meanCor[meanCor < thrCorPlant]
      return(data.frame(plotId = names(annPlotsCorr), reason = "mean corr",
                        value = annPlotsCorr))
    } else {
      return(NULL)
    }
  })
  ## Convert cormats to format used by ggplot.
  cormats <- lapply(X = cormats, FUN = function(cormat) {
    ## Set lower part of cormat to NA for plotting.
    cormat[lower.tri(cormat)] <- NA
    ## Melt to format used by ggplot.
    meltedCormat <- reshape2::melt(cormat, na.rm = TRUE)
    ## Convert Var1 and Var2 to factor needed for plotting.
    if (!is.factor(meltedCormat[["Var1"]])) {
      meltedCormat[["Var1"]] <- as.factor(meltedCormat[["Var1"]])
      meltedCormat[["Var2"]] <- as.factor(meltedCormat[["Var2"]])
    }
    attr(meltedCormat, which = "thrCor") <- attr(cormat, which = "thrCor")
    return(meltedCormat)
  })
  ## Find outlying plots based on PCA angles.
  annotatePlantsPcaAngle <- lapply(X = names(plantPcas), FUN = function(geno) {
    ## Calculate the pairwise difference of coordinates on the 2nd axis and
    ## annotate plant with average diff larger than threshold.
    plantPcaPlot <- factoextra::fviz_pca_var(plantPcas[[geno]])
    PcVecs <- as.matrix(plantPcaPlot$data[, 2:3])
    PcAngles <- matrix(nrow = nrow(PcVecs), ncol = nrow(PcVecs),
                       dimnames = list(rownames(PcVecs), rownames(PcVecs)))
    PcAngles[lower.tri(PcAngles)] <-
      combn(x = rownames(PcVecs), m = 2, FUN = function(plants) {
        angle(PcVecs[plants, ])
      }, simplify = TRUE)
    PcAngles[upper.tri(PcAngles)] <- t(PcAngles)[upper.tri(PcAngles)]
    meanPcAngles <- rowMeans(PcAngles, na.rm = TRUE)
    if (!is.null(geno.decomp)) {
      thrPcaPlant <- thrPca[attr(plantPcas[[geno]], which = "genoDecomp")]
    } else {
      thrPcaPlant <- thrPca
    }
    if (any(meanPcAngles >= thrPcaPlant)) {
      ## Create data.frame with info on plants with average difference
      ## above threshold.
      annPlotsPcAngle <- meanPcAngles[meanPcAngles >= thrPcaPlant]
      return(data.frame(plotId = names(annPlotsPcAngle), reason = "angle",
                        value = annPlotsPcAngle))
    } else {
      return(NULL)
    }
  })
  ## Find outlying plots based on slopes.
  annotatePlantsSlope <- lapply(X = names(slopemats), FUN = function(geno) {
    if (!is.null(slopemats[[geno]])) {
      ## Compute mean slope per plot on all but the worst slope per plant.
      ## This prevents one bad correlation causing all plants to be outliers.
      meanSlope <- apply(slopemats[[geno]], MARGIN = 1,
                         FUN = function(plant) {
                           mean(sort(plant)[-1])
                         })
      thrSlopePlant <- attr(slopemats[[geno]], which = "thrSlope")
    } else {
      meanSlope <- NULL
    }
    if (any(meanSlope < thrSlopePlant)) {
      ## Create data.frame with info on plants with average difference
      ## above threshold.
      annPlotsSlope <- meanSlope[meanSlope < thrSlopePlant]
      return(data.frame(plotId = names(annPlotsSlope), reason = "slope",
                        value = annPlotsSlope))
    } else {
      return(NULL)
    }
  })
  ## Convert slopemats to format used by ggplot.
  slopemats <- lapply(X = slopemats, FUN = function(slopemat) {
    ## Set lower part of cormat to NA for plotting.
    slopemat[lower.tri(slopemat)] <- NA
    ## Melt to format used by ggplot.
    meltedSlopemat <- reshape2::melt(slopemat, na.rm = TRUE)
    ## Convert Var1 and Var2 to factor needed for plotting.
    if (!is.factor(meltedSlopemat[["Var1"]])) {
      meltedSlopemat[["Var1"]] <- as.factor(meltedSlopemat[["Var1"]])
      meltedSlopemat[["Var2"]] <- as.factor(meltedSlopemat[["Var2"]])
    }
    attr(meltedSlopemat, which = "thrSlope") <-
      attr(slopemat, which = "thrSlope")
    return(meltedSlopemat)
  })
  ## Create full data.frame with annotated plants.
  annotatePlants <- do.call(rbind, c(annotatePlantsCor,
                                     annotatePlantsPcaAngle,
                                     annotatePlantsSlope))
  if (!is.null(annotatePlants)) {
    ## Merge genotype and geno.decomp to annotated plants.
    annotatePlants <- merge(unique(corrDatPred[c("genotype", geno.decomp,
                                                 "plotId")]),
                            annotatePlants)
    ## Order by genotype, geno.decomp and plotId.
    if (!is.null(geno.decomp)) {
      annOrd <- order(annotatePlants[["genotype"]],
                      annotatePlants[[geno.decomp]], annotatePlants[["plotId"]])
    } else {
      annOrd <- order(annotatePlants[["genotype"]], annotatePlants[["plotId"]])
    }
    annotatePlants <- droplevels(annotatePlants[annOrd, ])
  } else {
    annotatePlants <- data.frame()
  }
  plotInfo <- unique(corrDatPred[c("genotype", geno.decomp)])
  plotInfo <- plotInfo[!interaction(plotInfo) %in% genoPlotIdLim, , drop = FALSE]
  class(annotatePlants) <- c("serieOut", class(annotatePlants))
  attr(x = annotatePlants, which = "thrCor") <- thrCor
  attr(x = annotatePlants, which = "thrPca") <- thrPca
  attr(x = annotatePlants, which = "thrSlope") <- thrSlope
  attr(x = annotatePlants, which = "trait") <- trait
  attr(x = annotatePlants, which = "geno.decomp") <- geno.decomp
  attr(x = annotatePlants, which = "plotInfo") <- plotInfo
  attr(x = annotatePlants, which = "cormats") <- cormats
  attr(x = annotatePlants, which = "plantPcas") <- plantPcas
  attr(x = annotatePlants, which = "slopemats") <- slopemats
  attr(x = annotatePlants, which = "genoPreds") <- genoPreds
  attr(x = annotatePlants, which = "genoDats") <- genoDats
  return(annotatePlants)
}

#' Plot outliers for series of observations
#'
#' Plot the fitted spline, correlation matrix and PCA biplot for each of the
#' genotypes. Outlying series of observations are shown as filled dots in the
#' fitted spline plot, other observations are shown as open dots.
#'
#' @inheritParams plot.TP
#'
#' @param x An object of class \code{serieOut}.
#' @param ... Ignored.
#' @param reason A character vector indicating which types of outliers should
#' be plotted.
#' @param genotypes A character vector indicating which genotypes should be
#' plotted. If \code{NULL} all genotypes are plotted.
#' @param geno.decomp A character vector indicating which levels of
#' \code{geno.decomp} should be plotted. If \code{NULL} all levels are plotted.
#' Ignored if \code{geno.decomp} was not used when fitting models.
#' @param useTimeNumber Should the timeNumber be used instead of the timePoint
#' in the labels on the x-axis?
#' @param timeNumber If \code{useTimeNumber = TRUE}, a character vector
#' indicating the column containing the numerical time to use.
#'
#' @returns A list of ggplot objects is invisibly returned.
#'
#' @examples
#' \donttest{
#' ## The data from the Phenovator platform have been corrected for spatial
#' ## trends and outliers for single observations have been removed.
#'
#' ## Fit P-splines on a subset of genotypes
#' subGenoVator <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGenoVator,
#'                         knots = 50)
#'
#' ## Extract the data.frames with predicted values and P-Spline coefficients.
#' predDat <- fit.spline$predDat
#' coefDat <- fit.spline$coefDat
#'
#' ## The coefficients are then used to tag suspect time courses.
#' outVator <- detectSerieOut(corrDat = spatCorrectedVator,
#'                            predDat = predDat,
#'                            coefDat = coefDat,
#'                            trait = "EffpsII_corr",
#'                            genotypes = subGenoVator,
#'                            thrCor = 0.9,
#'                            thrPca = 30,
#'                            thrSlope = 0.7)
#'
#' ## The `outVator` can be visualized for selected genotypes.
#' plot(outVator, genotypes = "G151")
#'
#' ## Only visualize outliers tagged because of low correlation
#' ## between slopes of the regression.
#' plot(outVator, genotypes = "G151", reason = "slope")
#' }
#'
#' @family functions for detecting outliers for series of observations
#'
#' @export
plot.serieOut <- function(x,
                          ...,
                          reason = c("mean corr", "angle", "slope"),
                          genotypes = NULL,
                          geno.decomp = NULL,
                          useTimeNumber = FALSE,
                          timeNumber = NULL,
                          title = NULL,
                          output = TRUE) {
  if (useTimeNumber && (is.null(timeNumber) || !is.character(timeNumber) ||
                        length(timeNumber) > 1)) {
    stop("timeNumber should be a character string of length 1.\n")
  }
  thrCor <- attr(x = x, which = "thrCor")
  thrPca <- attr(x = x, which = "thrPca")
  trait <- attr(x = x, which = "trait")
  geno.decompVar <- attr(x = x, which = "geno.decomp")
  plotInfo <- attr(x = x, which = "plotInfo")
  ## Restrict x to selected reasons.
  reason <- match.arg(reason, several.ok = TRUE)
  x <- x[x[["reason"]] %in% reason, ]
  ## Restrict to selected levels of geno.decomp.
  if (!is.null(geno.decomp) && !is.null(geno.decompVar)) {
    if (!all(geno.decomp %in% plotInfo[[geno.decompVar]])) {
      stop("All selected geno.decomp levels should be in the data.\n")
    }
    plotInfo <- plotInfo[plotInfo[[geno.decompVar]] %in% geno.decomp, ]
  }
  if (!is.null(genotypes) && (!is.character(genotypes) ||
                              !all(genotypes %in% plotInfo[["genotype"]]))) {
    stop("genotypes should be a character vector of genotypes used for ",
         "outlier detection.\n")
  }
  if (is.null(genotypes)) {
    genotypes <- interaction(plotInfo[c("genotype", geno.decompVar)])
  } else {
    genotypes <- interaction(plotInfo[plotInfo[["genotype"]] %in% genotypes,
                                      c("genotype", geno.decompVar)])
  }
  genotypes <- as.character(genotypes)
  cormats <- attr(x = x, which = "cormats")[genotypes]
  plantPcas <- attr(x = x, which = "plantPcas")[genotypes]
  slopemats <- attr(x = x, which = "slopemats")[genotypes]
  genoPreds <- attr(x = x, which = "genoPreds")[genotypes]
  genoDats <- attr(x = x, which = "genoDats")[genotypes]
  ## Get minimum correlation. Cannot be higher than 0.8.
  minCor <- min(c(unlist(cormats, use.names = FALSE), 0.8), na.rm = TRUE)
  ## Get minimum slope. Cannot be higher than 0.8.
  minSlope <- min(c(unlist(slopemats, use.names = FALSE), 0.8), na.rm = TRUE)
  ## Compute the number of breaks for the time scale based on all plants.
  ## If there are less than 3 time points use the number of time points.
  ## Otherwise use 3.
  useTimePoint <- hasName(x = genoPreds[[1]], name = "timePoint") &&
    !useTimeNumber
  timeVar <- if (useTimePoint) "timePoint" else "timeNumber"
  timeVar2 <- if (useTimeNumber) timeNumber else timeVar
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
    plotShapes[names(plotShapes) %in% x[["plotId"]]] <- 21
    ## Reorder according to levels in genoDats.
    ## This prevents a duplicate legend in specific cases.
    plotShapes <- na.omit(plotShapes[levels(genoDats[[genotype]]$plotId)])
    ## Plot of time course per genotype: corrected data + spline per plant.
    kinetic <- ggplot2::ggplot(genoDats[[genotype]],
                               ggplot2::aes(x = .data[[timeVar2]],
                                            y = .data[[trait]],
                                            color = .data[["plotId"]])) +
      ggplot2::geom_point(ggplot2::aes(shape = .data[["plotId"]]), size = 2,
                          na.rm = TRUE, fill = "red") +
      ggplot2::geom_line(data = genoPreds[[genotype]],
                         ggplot2::aes(x = .data[[timeVar]],
                                      y = .data[["pred.value"]]),
                         linewidth = 0.5, na.rm = TRUE) +
      ggplot2::scale_shape_manual(values = plotShapes) +
      ggplot2::theme_light() +
      ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                     axis.title = ggplot2::element_text(size = 13))
    if (useTimePoint) {
      ## Format the time scale to Month + day.
      kinetic <- kinetic +
        ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                  labels = scales::date_format("%B %d"))
    }
    ## Correlation plot.
    if (any(c("mean corr", "slope") %in% reason)) {
      thrCorGeno <- attr(cormats[[genotype]], which = "thrCor")
      thrSlopeGeno <- attr(slopemats[[genotype]], which = "thrSlope")
      correl <- ggplot2::ggplot(data = cormats[[genotype]],
                                ggplot2::aes(x = .data[["Var2"]],
                                             y = .data[["Var1"]],
                                             fill = .data[["value"]])) +
        ## Move y-axis to the right.
        ggplot2::scale_y_discrete(position = "right") +
        ## Use coord fixed to create a square shaped output.
        ggplot2::coord_fixed() +
        ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
                       panel.border = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_line(color = "grey92"),
                       plot.title = ggplot2::element_text(hjust = 0.5),
                       axis.ticks = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                           hjust = 1),
                       legend.position = "left",
                       legend.box = "horizontal") +
        ggplot2::labs(title = "Correlations", x = NULL, y = NULL)
      if ("mean corr" %in% reason) {
        correl <- correl +
          ggplot2::geom_tile(color = "white") +
          ggplot2::scale_fill_gradientn(colors = c("red", "white", "blue"),
                                        values = scales::rescale(c(minCor,
                                                                   thrCorGeno, 1)),
                                        limits = c(minCor, 1),
                                        name = "Pearson\nCorrelation")
      }
      if ("slope" %in% reason) {
        correl <- correl +
          ggnewscale::new_scale_fill() +
          ## Add slope to upper left.
          ggplot2::geom_tile(data = slopemats[[genotype]],
                             ggplot2::aes(x = .data[["Var1"]],
                                          y = .data[["Var2"]],
                                          fill = .data[["value"]]),
                             color = "white") +
          ggplot2::scale_fill_gradientn(colors = c("cyan", "white", "darkgreen"),
                                        values = scales::rescale(c(minSlope,
                                                                   thrSlopeGeno, 1)),
                                        limits = c(minSlope, 1),
                                        name = "Slope\nCorrelation")
      }
    } else {
      correl <- grid::nullGrob()
    }
    ## PCA biplot.
    if ("angle" %in% reason) {
      pcaplot <- factoextra::fviz_pca_var(plantPcas[[genotype]])
    } else {
      pcaplot <- grid::nullGrob()
    }
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
    pGeno <-
      gridExtra::arrangeGrob(kinetic, correl, pcaplot,
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

#' Replace outliers for series of observations by NA
#'
#' Function for replacing outliers for series of observations in the data by NA.
#' The input can either be a data.frame, specified in \code{dat}, or the output
#' of the \code{fitSpline} function, specified in \code{fitSpline}. Exactly one
#' of these should be provided as input for the function.
#'
#' @param dat A \code{data.frame}.
#' @param fitSpline An object of class \code{HTPSpline}, the output of the
#' \code{\link{fitSpline}} function.
#' @param serieOut A data.frame with at least the column plotId with
#' values corresponding to those in dat/fitSpline.
#' @param reason A character vector indicating which types of outliers should
#' be replaced by NA.
#' @param traits The traits that should be replaced by NA. When using the
#' output of \code{detectSerieOut} as input for \code{serieOut} this defaults
#' to the trait used for when detecting the outliers.
#'
#' @returns Depending on the input either a \code{data.frame} or an object of
#' class \code{HTPSpline} for which the outliers specified in \code{serieOut}
#' are replaced by NA.
#'
#' @examples
#' ## Run the function to fit P-splines on a subset of genotypes.
#' subGenoVator <- c("G160", "G151")
#' fit.spline <- fitSpline(inDat = spatCorrectedVator,
#'                         trait = "EffpsII_corr",
#'                         genotypes = subGenoVator,
#'                         knots = 50)
#'
#' ## Extract the tables of predicted values and P-spline coefficients.
#' predDat <- fit.spline$predDat
#' coefDat <- fit.spline$coefDat
#'
#' ## The coefficients are then used to tag suspect time courses
#' outVator <- detectSerieOut(corrDat = spatCorrectedVator,
#'                            predDat = predDat,
#'                            coefDat = coefDat,
#'                            trait = "EffpsII_corr",
#'                            genotypes = subGenoVator,
#'                            thrCor = 0.9,
#'                            thrPca = 30,
#'                            thrSlope = 0.7)
#'
#' ## Replace the outliers by NA in the corrected data.
#' spatCorrectedVatorOut <- removeSerieOut(dat = spatCorrectedVator,
#'                                         serieOut = outVator)
#'
#' ## Only replace the slope outliers by NA in the corrected data.
#' spatCorrectedVatorOut2 <- removeSerieOut(dat = spatCorrectedVator,
#'                                         serieOut = outVator,
#'                                         reason = "slope")
#'
#' ## Replace the outliers by NA in the corrected data.
#' ## Replace both the corrected value and the raw trait value by NA.
#' spatCorrectedVatorOut3 <-
#'   removeSerieOut(dat = spatCorrectedVator,
#'                  serieOut = outVator,
#'                  traits = c("EffpsII", "EffpsII_corr"))
#'
#' @family functions for detecting outliers for series of observations
#'
#' @export
removeSerieOut <- function(dat = NULL,
                           fitSpline = NULL,
                           serieOut,
                           reason = c("mean corr", "angle", "slope"),
                           traits = attr(x = serieOut,
                                         which = "trait")) {
  reason <- match.arg(reason, several.ok = TRUE)
  ## Check that one of dat and fitSpline are specified.
  if ((is.null(dat) && is.null(fitSpline)) || (
    !is.null(dat) && !is.null(fitSpline))) {
    stop("Specify exactly one of dat and fitSpline as inputs.\n")
  }
  if (!is.null(dat)) {
    if (!inherits(dat, "data.frame")) {
      stop("dat should be a data.frame.\n")
    }
    if (!hasName(dat, "plotId")) {
      stop("dat should at least contain the column plotId.\n")
    }
  }
  if (!is.null(fitSpline)) {
    if (!inherits(fitSpline, "HTPSpline")) {
      stop("fitSpline should be an object of class HTPSpline.\n")
    }
  }
  if (!inherits(serieOut, "data.frame")) {
    stop("serieOut should be a data.frame.\n")
  }
  if (nrow(serieOut) > 0) {
    if (!hasName(serieOut, "plotId")) {
      stop("serieOut should at least contain the column plotId.\n")
    }
    if (length(reason) < 3)
      if (!hasName(serieOut, "reason")) {
        stop("serieOut should contain a column reason.\n")
      } else {
        serieOut <- serieOut[serieOut[["reason"]] %in% reason, ]
      }
    if (!is.null(dat)) {
      ## Check if all traits are present in dat.
      traitMiss <- traits[!traits %in% colnames(dat)]
      if (length(traitMiss) > 0) {
        stop("All traits should be in dat. The following traits are missing:\n",
             paste(traitMiss, collapse = ", "))
      }
      ## Remove plots that are in serieOut.
      for (trait in traits) {
        dat[dat[["plotId"]] %in% serieOut[["plotId"]], trait] <- NA
      }
    } else if (!is.null(fitSpline)) {
      fitSpline$coefDat[fitSpline$coefDat[["plotId"]] %in%
                          serieOut[["plotId"]], "obj.coefficients"] <- NA
      fitSpline$predDat[fitSpline$predDat[["plotId"]] %in% serieOut[["plotId"]],
                        c("pred.value", "deriv", "deriv2")] <- NA
      modDat <- attr(x = fitSpline, which = "modDat")
      for (trait in traits) {
        modDat[modDat[["plotId"]] %in% serieOut[["plotId"]], trait] <- NA
      }
      attr(x = fitSpline, which = "modDat") <- modDat
    }
  }
  res <- if (!is.null(dat)) dat else fitSpline
  return(res)
}
