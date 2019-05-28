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
  yMin <- min(trDat$rowCoord)
  yMax <- max(trDat$rowCoord)
  xMin <- min(trDat$colCoord)
  xMax <- max(trDat$colCoord)
  ## Create matrix containing replicates.
  ## First create an empty matrix containing all row/column values
  ## between min and max to assure complete missing rows/columns
  ## are added.
  M <- matrix(nrow = yMax - yMin + 1, ncol = xMax - xMin + 1,
              dimnames = list(yMin:yMax, xMin:xMax))
  for (i in 1:nrow(trDat)) {
    M[as.character(trDat[i, "rowCoord"]),
      as.character(trDat[i, "colCoord"])] <- trDat[i, bordVar]
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
