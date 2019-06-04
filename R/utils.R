#' Helper function for constructing two data.frames containing the coordinates
#' that can be used for plotting a border around parts of a raster plot using
#' geom_path in ggplot2. This can be used to construct an outside border
#' around each replicate in a plot. ggplot2 itself doesn't have this
#' functionality.
#'
#' @noRd
#' @keywords internal
calcPlotBorders <- function(trDat,
                            bordVar) {
  yMin <- min(trDat$rowNum)
  yMax <- max(trDat$rowNum)
  xMin <- min(trDat$colNum)
  xMax <- max(trDat$colNum)
  ## Create matrix containing replicates.
  ## First create an empty matrix containing all row/column values
  ## between min and max to assure complete missing rows/columns
  ## are added.
  M <- matrix(nrow = yMax - yMin + 1, ncol = xMax - xMin + 1,
              dimnames = list(yMin:yMax, xMin:xMax))
  for (i in 1:nrow(trDat)) {
    M[as.character(trDat[i, "rowNum"]),
      as.character(trDat[i, "colNum"])] <- trDat[i, bordVar]
  }
  ## Create an imputed version of M for plotting borders around NA values.
  MImp <- M
  MImp[is.na(MImp)] <- nlevels(trDat[[bordVar]]) + 1
  has.breaks <- function(x) {
    ncol(x) == 2 & nrow(x) > 0
  }
  ## Create a data.frame with positions where the value of rep in the
  ## data changes in vertical direction.
  vertW <- do.call(rbind.data.frame,
                   Filter(f = has.breaks, x = Map(function(i, x) {
                     cbind(y = i, x = which(diff(c(0, x, 0)) != 0))
                   }, 1:nrow(MImp), split(MImp, 1:nrow(MImp)))))
  ## Remove vertical walls that are on the outside bordering an NA value
  ## to prevent drawing of unneeded lines.
  vertW <- vertW[!(vertW$x == 1 & is.na(M[vertW$y, 1])) &
                   !(vertW$x == ncol(M) + 1 &
                       is.na(M[vertW$y, ncol(M)])), ]
  ## Add min row value for plotting in the correct position.
  vertW$y <- vertW$y + yMin - 1
  vertW$x <- vertW$x + xMin - 1
  ## For horizontal walls follow the same procedure as above.
  horW <- do.call(rbind.data.frame,
                  Filter(f = has.breaks, x = Map(function(i, y) {
                    cbind(x = i, y = which(diff(c(0, y, 0)) != 0))
                  }, 1:ncol(MImp), as.data.frame(MImp))))
  horW <- horW[!(horW$y == 1 & is.na(M[1, horW$x])) &
                 !(horW$y == nrow(M) + 1 & is.na(M[nrow(M), horW$x])), ]
  horW$y <- horW$y + yMin - 1
  horW$x <- horW$x + xMin - 1
  return(list(horW = horW, vertW = vertW))
}

#' @noRd
#' @keywords internal
addMissVals <- function(dat,
                        trait) {
  ## Create lhs formula for dcast using all columns in dat except
  ## timePoint (rhs) and trait (value var).
  castCols <- setdiff(colnames(dat), c("timePoint", trait))
  castForm <- formula(paste(paste(castCols, collapse = "+"), "~ timePoint"))
  ## Melt and reshape with default settings adds missing combinations to the
  ## data table.
  datOut <- reshape2::melt(data = reshape2::dcast(data = dat,
                                                  formula = castForm,
                                                  value.var = trait),
                           id.vars = castCols, variable.name = "timePoint",
                           value.name = trait)
  ## Melt loses date format for timePoint so resetting it here.
  datOut[["timePoint"]] <- lubridate::as_datetime(datOut[["timePoint"]])
  return(datOut)
}

#' @noRd
#' @keywords internal
xyFacetPlot <- function(baseDat,
                        overlayDat = NULL,
                        xVal = "timePoint",
                        yVal,
                        yValOverlay = yVal,
                        groupVal = "plotId",
                        colVal = groupVal,
                        facetVal = "genotype",
                        title,
                        xLab = "Time",
                        yLab = "Trait",
                        output = TRUE) {
  p <- ggplot2::ggplot(baseDat,
                       ggplot2::aes_string(x = xVal, y = yVal)) +
    ggplot2::geom_line(ggplot2::aes_string(group = groupVal, color = colVal),
                       show.legend = FALSE, na.rm = TRUE) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0, "cm"),
                   panel.border = ggplot2::element_rect(color = "black",
                                                        fill = "transparent"),
                   strip.background = ggplot2::element_rect(color = "black",
                                                            fill = "bisque"),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(title = title, x = xLab, y = yLab)
  if (!is.null(overlayDat)) {
    p <- p + ggplot2::geom_line(ggplot2::aes_string(x = "timePoint",
                                                    y = yValOverlay),
                                data = overlayDat, color = "black", size = 1,
                                show.legend = FALSE, na.rm = TRUE)
  }
  for (i in 1:ceiling(nlevels(baseDat[[facetVal]]) / 25)) {
    plot(p + ggforce::facet_wrap_paginate(facets = facetVal, nrow = 5, ncol = 5,
                                          page = i))
  }
}
