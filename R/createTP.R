#' @export
createTimePoints <- function(dat,
                             genotype,
                             timePoint,
                             plotId = NULL,
                             repId = NULL,
                             rowNum = NULL,
                             colNum = NULL,
                             rowId = rowNum,
                             colId = colNum,
                             addCheck = FALSE,
                             checkGenotypes = NULL) {
  ## Save name of original dat for naming output.
  datName <- deparse(substitute(dat))
  if (length(datName) > 1) {
    datName <- "dat"
  }
  ## Checks.
  if (missing(dat) || !is.data.frame(dat)) {
    stop("dat has to be a data.frame.\n")
  }
  ## Convert input to data.frame. This needs to be done to be able to handle
  ## tibbles and possibly other data structures in the future.
  dat <- as.data.frame(dat)
  cols <- colnames(dat)
  for (param in c(genotype, timePoint, plotId, repId, rowId, colId,
                  rowNum, colNum)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !hasName(dat, param))) {
      stop(paste(deparse(param), "has to be NULL or a column in dat.\n"))
    }
  }
  ## Create list of reserved column names for renaming columns.
  renameCols <- c("genotype", "timePoint", "plotId", "repId", "rowId", "colId",
                  "rowNum", "colNum")
  ## First rename duplicate colums and add duplicated columns to dat
  renameFrom <- as.character(sapply(X = renameCols, FUN = function(x) {
    get(x)
  }))
  ## Create a data.frame with renamed cols to add to TP as an attribute.
  renamed <- data.frame(orig = renameFrom[renameFrom != "NULL"],
                        new = renameCols[renameFrom != "NULL"],
                        stringsAsFactors = FALSE)
  ## Get duplicate columns.
  dupCols <- which(duplicated(renameFrom) & renameFrom != "NULL")
  for (dupCol in dupCols) {
    ## Copy original column as extra column in dat for each duplicate.
    tempName <- paste0(".temp", dupCol)
    dat[tempName] <- dat[, colnames(dat) == renameFrom[dupCol]]
    ## Add new replacementname to cols and renameFrom.
    cols[length(cols) + 1] <- tempName
    renameFrom[dupCol] <- tempName
  }
  ## Rename columns.
  for (i in 1:length(renameCols)) {
    cols[cols == renameFrom[i]] <- renameCols[i]
  }
  colnames(dat) <- cols
  ## Convert timepoint to nice format.
  dat[["timePoint"]] <- lubridate::as_datetime(dat[["timePoint"]])
  ## Add check genotype.
  if (addCheck) {
    dat[["check"]] <- ifelse(dat[["genotype"]] %in% checkGenotypes,
                             dat[["genotype"]], "noCheck")
    dat[["genoCheck"]] <- dat[["genotype"]]
    dat[["genoCheck"]][dat[["check"]] != "noCheck"] <- NA
  }
  ## Convert columns to factor if neccessary.
  factorCols <-  c("genotype", "plotId", "repId", "rowId", "colId", "check",
                   "genoCheck")
  for (factorCol in factorCols) {
    if (hasName(dat, factorCol)) {
      dat[[factorCol]] <- as.factor(dat[[factorCol]])
    }
  }
  ## Convert columns to numeric if neccessary.
  numCols <- c("rowNum", "colNum")
  for (numCol in numCols) {
    if (hasName(dat, numCol) && !is.numeric(dat[cols == numCol])) {
      dat[cols == numCol] <- as.numeric(dat[, cols == numCol])
    }
  }
  ## Check plotId for uniqueness.
  if (is.null(plotId)) {
    warning("No unique plot identifier supplied.\n")
  } else {
    if (max(by(data = dat[["plotId"]], INDICES = dat[["timePoint"]],
               FUN = anyDuplicated) > 0)) {
      warning("plotId should be unique per time point.\n")
    }
  }
  listData <- split(x = dat, f = dat[["timePoint"]])
  ## Set meta for all timePoints in dat.
  for (tr in names(listData)) {
    ## Add a list of columns that have been renamed as attribute to TP.
    attr(x = listData[[tr]], which = "renamedCols") <-
      if (nrow(renamed) > 0) renamed else NULL
  }
  TP <- structure(listData,
                  class = c("TP", "list"))
  return(TP)
}

#' Plot function for class TP
#'
#' Plotting function for objects of class TP. Plots either the layout of the
#' different timePoints within the TP object or locates the timePoints on a map. Also a
#' boxplot can be made for selected traits and timePoints and a plot of correlations
#' between timePoints. A detailed description and optional extra parameters of the
#' different plots is given in the sections below.
#'
#' @section Layout Plot:
#' Plots the layout of the selected timePoints (all available timePoints by default).
#' This plot can only be made for timePoints that contain both row (\code{rowNum})
#' and column (\code{colNum}) information. If either one of those is missing
#' the timePoint is skipped with a warning. If blocks (\code{subBlock}) are
#' available for a timePoint these are indicated in different colors per block,
#' otherwise all plots are colored in grey. If replicates (\code{repId}) are
#' available a black line is plotted between diffent replicates. Missing plots
#' are indicated in white. This can either be single plots in a timePoint or
#' complete missing columns or rows.\cr
#' Extra parameter options:
#' \describe{
#' \item{showGeno}{Should individual genotypes be indicated in the plot?
#' Defaults to \code{FALSE}}
#' \item{highlight}{A character vector of genotypes to be highlighted in the
#' plot.}
#' \item{colorSubBlock}{Should subBlocks be colored with a different color per
#' subBlock? Defaults to \code{FALSE}. \code{colorSubBlock} is ignored when
#' highlight is used to highlight genotypes.}
#' }
#'
#' @section Box Plot:
#' Creates a boxplot per selected trait grouped by timePoint. Extra parameter
#' options:
#' \describe{
#' \item{groupBy}{A character string indicating a column in \code{TP} by which
#' the boxes in the plot should be grouped. By default the boxes are grouped
#' per timePoint.}
#' \item{colorBy}{A character string indicating a column in \code{TP} by which
#' the boxes are colored. Coloring will be done within the groups indicated by
#' the \code{groupBy} parameter.}
#' \item{orderBy}{A character string indicating the way the boxes should be
#' ordered. Either "alphabetic" for alphabetical ordering of the groups,
#' "ascending" for ordering by ascending mean, or "descending" for ordering by
#' descending mean. Default boxes are ordered alphabetically.}
#' }
#'
#' @section Correlation Plot:
#' Draws a heatmap of correlations between timePoints per selected trait. If
#' genotypes are replicated within timePoints genotypic means are taken before
#' computing correlations. The order of the timePoints in the heatmap is determined
#' by clustering them.
#'
#' @param x An object of class TP.
#' @param ... Extra plot options. Described per plotType in their respective
#' section.
#' @param plotType A single character string indicating which plot should be
#' made. See the sections below for a detailed explanation of the plots.
#' @param timePoints A character vector indicating the timePoints to be plotted.
#' Only used if \code{plotType} = "layout" or "box".
#' @param traits A character vector indicating the traits to be plotted in
#' a boxplot. Only used if \code{plotType} = "box" or "cor".
#' @param genotypes A character vector indicating the genotypes to be plotted.
#' Only used if \code{plotType} = "raw".
#' @param geno.decomp A character vector indicating the grouping of the
#' genotypes to be plotted. Only used if \code{plotType} = "raw".
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @family functions for TP objects
#'
#' @export
plot.TP <- function(x,
                    ...,
                    plotType = c("layout", "box", "cor", "raw"),
                    timePoints = names(x),
                    traits = NULL,
                    genotypes = NULL,
                    geno.decomp = NULL,
                    output = TRUE) {
  ## Checks.
  if (!is.character(timePoints) || !all(hasName(x = x, name = timePoints))) {
    stop(paste0("All timePoints should be in ", deparse(substitute(x)), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (plotType == "layout") {
    showGeno <- isTRUE(dotArgs$showGeno)
    highlight <- dotArgs$highlight
    colorSubBlock <- isTRUE(dotArgs$colorSubBlock)
    if (!is.null(highlight) && !is.character(highlight)) {
      stop("highlight should be a character vector.\n")
    }
    p <- setNames(vector(mode = "list", length = length(timePoints)), timePoints)
    for (timePoint in timePoints) {
      trDat <- x[[timePoint]]
      if (!"rowNum" %in% colnames(trDat)) {
        warning(paste0("rowNum should be a column in ", timePoint, ".\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (!"colNum" %in% colnames(trDat)) {
        warning(paste0("colNum should be a column in ", timePoint, ".\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (length(highlight) > 0) {
        trDat$highlight. <- ifelse(trDat$genotype %in% highlight,
                                   as.character(trDat$genotype), NA)
      }
      trLoc <- attr(trDat, "trLocation")
      ylen <- attr(trDat, "trPlLength")
      xlen <- attr(trDat, "trPlWidth")
      ## Compute aspect for proper depiction of field size. If no information
      ## is available plots are assumed to be square.
      if (is.null(ylen) || is.null(xlen)) {
        aspect <- length(unique(trDat$colNum)) /
          length(unique(trDat$rowNum))
      } else {
        aspect <- ylen / xlen
      }
      plotRep <- hasName(x = trDat, name = "repId")
      plotSubBlock <- hasName(x = trDat, name = "subBlock")
      ## Create data for lines between replicates.
      if (plotRep) {
        repBord <- calcPlotBorders(trDat = trDat, bordVar = "repId")
      }
      ## Create base plot.
      pTr <- ggplot2::ggplot(data = trDat,
                             ggplot2::aes_string(x = "colNum",
                                                 y = "rowNum")) +
        ggplot2::coord_fixed(ratio = aspect,
                             xlim = range(trDat$colNum) + c(-0.5, 0.5),
                             ylim = range(trDat$rowNum) + c(-0.5, 0.5),
                             clip = "off") +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ## Move ticks to edge of the plot.
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(),
                                    expand = c(0, 0)) +
        ggplot2::ggtitle(trLoc)
      if (sum(!is.na(trDat$highlight.)) > 0) {
        ## Genotypes to be highlighted get a color.
        ## Everything else the NA color.
        pTr <- pTr + ggplot2::geom_tile(
          ggplot2::aes_string(fill = "highlight."), color = "grey75") +
          ggplot2::labs(fill = "Highlighted") +
          ## Remove NA from scale.
          ggplot2::scale_fill_discrete(na.translate = FALSE)
      } else if (plotSubBlock && colorSubBlock) {
        ## Color tiles by subblock.
        pTr <- pTr + ggplot2::geom_tile(
          ggplot2::aes_string(fill = "subBlock"), color = "grey75") +
          ggplot2::guides(fill = ggplot2::guide_legend(ncol = 3))
      } else {
        ## No subblocks and no hightlights so just a single fill color.
        pTr <- pTr + ggplot2::geom_tile(fill = "white", color = "grey75")
      }
      ## Create data for lines between subBlocks.
      if (plotSubBlock) {
        subBlockBord <- calcPlotBorders(trDat = trDat, bordVar = "subBlock")
        pTr <- pTr +
          ## Add verical lines as segment.
          ## adding/subtracting 0.5 assures plotting at the borders of
          ## the tiles.
          ggplot2::geom_segment(
            ggplot2::aes_string(x = "x - 0.5", xend = "x - 0.5",
                                y = "y - 0.5", yend = "y + 0.5",
                                linetype = "'subBlocks'"),
            data = subBlockBord$vertW, size = 0.4) +
          ggplot2::geom_segment(
            ggplot2::aes_string(x = "x - 0.5", xend = "x + 0.5",
                                y = "y - 0.5", yend = "y - 0.5"),
            data = subBlockBord$horW, size = 0.4)
      }
      if (showGeno) {
        ## Add names of genotypes to the center of the tiles.
        pTr <- pTr + ggplot2::geom_text(ggplot2::aes_string(label = "genotype"),
                                        size = 2, check_overlap = TRUE)
      }
      if (plotRep) {
        ## Add lines for replicates.
        pTr <- pTr +
          ## Add verical lines as segment.
          ## adding/subtracting 0.5 assures plotting at the borders of
          ## the tiles.
          ggplot2::geom_segment(
            ggplot2::aes_string(x = "x - 0.5", xend = "x - 0.5",
                                y = "y - 0.5", yend = "y + 0.5",
                                linetype = "'replicates'"),
            data = repBord$vertW, size = 1) +
          ggplot2::geom_segment(
            ggplot2::aes_string(x = "x - 0.5", xend = "x + 0.5",
                                y = "y - 0.5", yend = "y - 0.5"),
            data = repBord$horW, size = 1)
      }
      if (plotSubBlock || plotRep) {
        shwVals <- c(plotRep, plotSubBlock)
        pTr <- pTr +
          ## Add a legend entry for replicates and subBlocks.
          ggplot2::scale_linetype_manual(c("replicates", "subBlocks")[shwVals],
                                         values = c("replicates" = "solid",
                                                    "subBlocks" = "solid")[shwVals],
                                         name = ggplot2::element_blank()) +
          ggplot2::guides(linetype = ggplot2::guide_legend(override.aes =
                                                             list(size = c(1, 0.4)[shwVals])))
      }
      p[[timePoint]] <- pTr
      if (output) {
        plot(pTr)
      }
    }
  } else if (plotType == "box") {
    if (is.null(traits) || !is.character(traits)) {
      stop("traits should be a character vector.\n")
    }
    groupBy <- dotArgs$groupBy
    if (!is.null(groupBy) && (!is.character(groupBy) || length(groupBy) > 1)) {
      stop("groupBy should be a single character string.\n")
    }
    if (!is.null(groupBy) && !all(sapply(X = x, FUN = function(timePoint) {
      hasName(x = timePoint, name = groupBy)
    }))) {
      stop("groupBy should be a column in TP.\n")
    }
    colorBy <- dotArgs$colorBy
    if (!is.null(colorBy) && (!is.character(colorBy) || length(colorBy) > 1)) {
      stop("colorBy should be a single character string.\n")
    }
    if (!is.null(colorBy) && !all(sapply(X = x, FUN = function(timePoint) {
      hasName(x = timePoint, name = colorBy)
    }))) {
      stop("colorBy should be a column in TP.\n")
    }
    orderBy <- dotArgs$orderBy
    if (!is.null(orderBy)) {
      orderBy <- match.arg(orderBy, choices = c("alphabetic", "ascending",
                                                "descending"))
    } else {
      orderBy <- "alphabetic"
    }
    p <- setNames(vector(mode = "list", length = length(traits)), traits)
    for (trait in traits) {
      ## Create a single data.frame from x with only columns timePoint, trait and
      ## genotype. Genotype is needed to be able to display hovering info.
      ## timePoints where trait is not measured/available are removed by setting
      ## them to NULL.
      xVar <- if (is.null(groupBy)) "timePoint" else groupBy
      plotDat <- Reduce(f = rbind,
                        x = lapply(X = x[timePoints],
                                   function(timePoint) {
                                     if (!hasName(x = timePoint, name = trait)) {
                                       NULL
                                     } else {
                                       if (!hasName(x = timePoint,
                                                    name = "timePoint")) {
                                         timePoint[["timePoint"]] <- names(x)
                                       }
                                       timePoint[c(trait, "genotype", xVar,
                                                   if (!is.null(colorBy)) colorBy)]
                                     }
                                   }))
      if (is.null(plotDat)) {
        warning(paste0(trait, " isn't a column in any of the timePoints.\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (orderBy != "alphabetic") {
        ## Reorder levels in timePoint so plotting is done according to orderBy.
        levNw <- reorder(x = plotDat[[xVar]], X = plotDat[[trait]],
                         FUN = mean, na.rm = TRUE, order = TRUE)
        if (orderBy == "ascending") {
          plotDat[xVar] <- factor(plotDat[[xVar]], levels = levels(levNw))
        } else {
          plotDat[xVar] <- factor(plotDat[[xVar]], levels = rev(levels(levNw)))
        }
      } else {
        ## Always convert timepoint to factor.
        plotDat[xVar] <- factor(plotDat[[xVar]])
      }
      ## Create boxplot.
      ## Back ticks around variable names are needed to handle spaces in names.
      pTr <- ggplot2::ggplot(plotDat,
                             ggplot2::aes_string(x = paste0("`", xVar, "`"),
                                                 y = paste0("`", trait, "`"),
                                                 fill = if (is.null(colorBy)) NULL else
                                                   paste0("`", colorBy, "`"))) +
        ggplot2::geom_boxplot(na.rm = TRUE) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5,
                                                           hjust = 1)) +
        ggplot2::labs(x = xVar, y = trait)
      p[[trait]] <- pTr
      if (output) {
        plot(pTr)
      }
    }
  } else if (plotType == "cor") {
    if (length(timePoints) == 1) {
      stop("At least two timePoints requiered for a correlation plot.\n")
    }
    if (is.null(traits) || !is.character(traits)) {
      stop("traits should be a character vector.\n")
    }
    p <- setNames(vector(mode = "list", length = length(traits)), traits)
    for (trait in traits) {
      ## Create a single data.frame from x with only columns timePoint and trait.
      ## timePoints where trait is not measured/available are removed by setting
      ## them to NULL.
      plotDat <- Reduce(f = rbind, x = lapply(X = x, FUN = function(timePoint) {
        if (!hasName(x = timePoint, name = trait)) {
          NULL
        } else {
          timePoint[c("genotype", "timePoint", trait)]
        }
      }))
      if (is.null(plotDat)) {
        warning(paste0(trait, " isn't a column in any of the timePoints.\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      ## Create table with values trait per genotype per timePoint.
      ## If TP already contains BLUEs/BLUPs taking means doesn't do anything
      ## but it is needed for raw data where there can be replicates.
      plotTab <- tapply(plotDat[[trait]],
                        INDEX = list(plotDat[["genotype"]],
                                     plotDat[["timePoint"]]),
                        FUN = mean, na.rm = TRUE)
      ## Create a correlation matrix.
      corMat <- cor(plotTab, use = "pairwise.complete.obs")
      ## Remove rows and columns with only NA.
      corKeep <- sapply(X = 1:ncol(corMat), FUN = function(i) {
        any(!is.na(corMat[, i]))
      })
      corMat <- corMat[corKeep, corKeep, drop = FALSE]
      ## Melt to get the proper format for ggplot.
      meltedCorMat <- reshape2::melt(corMat)
      ## If timePoint names consist of only numbers melt converts them to numeric.
      ## This gives problems with plotting, so reconvert them to factor.
      if (is.numeric(meltedCorMat[["Var1"]])) {
        meltedCorMat[["Var1"]] <- factor(meltedCorMat[["Var1"]],
                                         levels = rownames(corMat))
        meltedCorMat[["Var2"]] <- factor(meltedCorMat[["Var2"]],
                                         levels = rownames(corMat))
      }
      ## Remove top left of the plot. Only plotting a bottom right triangle.
      ## Diagonal is removed as well.
      meltedCorMat <- meltedCorMat[as.numeric(meltedCorMat$Var1) >
                                     as.numeric(meltedCorMat$Var2), ]
      ## Create plot.
      pTr <- ggplot2::ggplot(data = meltedCorMat,
                             ggplot2::aes_string(x = "Var1", y = "Var2",
                                                 fill = "value")) +
        ggplot2::geom_tile(color = "grey50") +
        ## Create a gradient scale.
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                      na.value = "grey", limit = c(-1, 1)) +
        ## Move y-axis to the right for easier reading.
        ggplot2::scale_y_discrete(position = "right") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           vjust = 1, size = 6,
                                                           hjust = 1)) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
        ## Center title.
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ## Remove grid behind empty bit of triangle.
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank()) +
        ## No axis and legend titles.
        ggplot2::labs(x = "", y = "", color = "") +
        ggplot2::ggtitle(paste("Correlations of timepoints for", trait)) +
        ggplot2::guides(size = FALSE) +
        ## Fix coordinates to get a square sized plot.
        ggplot2::coord_fixed()
      p[[trait]] <- pTr
      if (output) {
        plot(pTr)
      }
    }
  } else if (plotType == "raw") {
    if (is.null(traits) || !is.character(traits)) {
      stop("traits should be a character vector.\n")
    }
    if (!is.null(geno.decomp) && !all(sapply(X = x, FUN = function(timePoint) {
      hasName(x = timePoint, name = geno.decomp)
    }))) {
      stop("geno.decomp should be a column in TP.\n")
    }
    p <- setNames(vector(mode = "list", length = length(traits)), traits)
    for (trait in traits) {
      ## Create a single data.frame from x with only columns timePoint and
      ## trait. timePoints where trait is not measured/available are removed
      ## by setting them to NULL.
      plotDat <- Reduce(f = rbind, x = lapply(X = x, FUN = function(timePoint) {
        if (!hasName(x = timePoint, name = trait)) {
          NULL
        } else {
          timePoint[c("genotype", "timePoint", "plotId", trait, geno.decomp)]
        }
      }))
      if (is.null(plotDat)) {
        warning(paste0(trait, " isn't a column in any of the timePoints.\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      ## Restrict genotypes.
      if (!is.null(genotypes)) {
        plotDat <- plotDat[plotDat[["genotype"]] %in% genotypes, ]
        plotDat <- droplevels(plotDat)
      }
      ## Add interaction with geno.decomp to genotypes.
      if (!is.null(geno.decomp)) {
        plotDat[["genotype"]] <- interaction(plotDat[[geno.decomp]],
                                             plotDat[["genotype"]], sep = "_")
      }
      ## Add combinations missing in data to plotDat.
      plotDat <- addMissVals(dat = plotDat, trait = trait)
      ## Create actual plots.
      xyFacetPlot(baseDat = plotDat, yVal = trait,
                  title = "Phenovator platform - Raw data", yLab = trait,
                  output = output)
    }
  }
  invisible(p)
}
