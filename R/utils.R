#' Custom tryCatch to return result, errors and warnings.
#' Copied from http://stackoverflow.com/a/24569739/2271856.
#'
#' @noRd
#' @keywords internal
tryCatchExt <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }), warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warning = warn, error = err)
}

#' Helper function for suppressing a single warning message.
#'
#' @noRd
#' @keywords internal
supprWarn <- function(expression,
                      message) {
  withCallingHandlers(expression,
                      warning = function(w) {
                        if (grepl(message, w$message)) {
                          invokeRestart("muffleWarning")
                        }
                      })
}

#' Helper function for checking whether error message about 1% change on
#' last iteration for asreml is worth mentioning as a warning.
#' If the corresponding parameter is close to zero then changes of 1%
#' or more can be expected and are ok.
#'
#' @noRd
#' @keywords internal
chkLastIter <- function(model) {
  wrnMsg <- "changed by more than 1%"
  if (any(grepl(pattern = wrnMsg, x = model$warning))) {
    ## EXtract monitor df from model object.
    mon <- model$value$monitor
    ## Extract values for parameters for last 2 iterations.
    ## First 3 rows give general model info. Last col a summary.
    lastIt <- mon[-(1:3), c(ncol(mon) - 2, ncol(mon) - 1)]
    ## Compute change of parameters in last iteration.
    change <- ifelse(lastIt[, 1] == 0, 0, abs((lastIt[, 2] - lastIt[, 1]) /
                                                lastIt[, 1]) * 100)
    ## Suppress warning if the change was less than 5% or the param value less
    ## than 0.1.
    if (all(change <= 5) || all(lastIt[change > 5, 1] < 0.1)) {
      model$warning <- model$warning[!grepl(pattern = wrnMsg,
                                            x = model$warning)]
    }
  }
  return(model)
}

#' Helper function for converting certain asreml warnings to errors.
#'
#' @noRd
#' @keywords internal
wrnToErr <- function(model) {
  wrns <- c("Abnormal termination", "returning -Inf")
  for (wrn in wrns) {
    if (any(grepl(pattern = wrn, x = model$warning))) {
      ## Remove from warnings and add to errors
      model$error <- c(model$error, model$warning[grepl(pattern = wrn,
                                                        x = model$warning)])
      model$warning <- model$warning[!grepl(pattern = wrn,
                                            x = model$warning)]
    }
  }
  return(model)
}

#' Extended version of asreml.predict
#'
#' Asreml has a bug that may throw a warning message:
#' Abnormal termination
#' Insufficient workspace - (reset workspace or pworkspace arguments)
#' This may be avoided by increasing pworkspace, but this doesn't
#' always work.
#' If this happens pworkspace is increased in 'small' steps.
#'
#' @noRd
#' @keywords internal
predictAsreml <- function(model,
                          classify = "genotype",
                          associate = as.formula("~ NULL"),
                          vcov = TRUE,
                          ...) {
  wrnMsg <- "reset workspace or pworkspace arguments"
  ## Predict using default settings, i.e. pworkspace = 8e6
  modelP <- tryCatchExt(predict(model, classify = classify,
                                vcov = vcov, associate = associate,
                                maxiter = 20, trace = FALSE, ...))
  pWorkSpace <- 8e6
  ## While there is a warning, increase pWorkSpace and predict again.
  while (!is.null(modelP$warning) &&
         any(grepl(pattern = wrnMsg, x = modelP$warning))
         && pWorkSpace < 160e6) {
    pWorkSpace <- pWorkSpace + 8e6
    modelP <- tryCatchExt(predict(model, classify = classify,
                                  vcov = vcov, associate = associate,
                                  maxiter = 20, pworkspace = pWorkSpace,
                                  trace = FALSE, ...))
  }
  if (!is.null(modelP$warning) && !all(grepl(pattern = wrnMsg,
                                             x = modelP$warning))) {
    modelP <- chkLastIter(modelP)
    if (length(modelP$warning) != 0) {
      warning(modelP$warning, call. = FALSE)
    }
  }
  if ((length(modelP$warning) == 0 ||
       !all(grepl(pattern = wrnMsg, x = modelP$warning))) &&
      is.null(modelP$error)) {
    return(modelP$value)
  } else {
    stop(paste("Error in asreml when running predict. Asreml message:\n",
               modelP$error, "\n",
               modelP$warning, "\n"), call. = FALSE)
  }
}

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

#' @noRd
#' @keywords internal
chkFile <- function(outFile) {
  if (!is.character(outFile) || length(outFile) > 1 ||
      tools::file_ext(outFile) != "csv") {
    stop("outFile should be a single character string ending in .csv.\n")
  }
  if (file.access(dirname(outFile), 2)) {
    stop("no permission to write to ", outFile, ".\n")
  }
}

#' Helper function for checking timepoints
#'
#' Helper function that checks if the timePoints provided are in the object
#' (this can be an object of class TP or fitMod). If the timePoints provided are
#' in a numeric format they are converted to their corresponding character
#' values.
#'
#' @param x An R object.
#' @param timePoints a character or numeric vector containing timePoints.
#'
#' @return A character vector containing timePoints.
#'
#' @noRd
#' @keywords internal
chkTimePoints <- function(x,
                          timePoints) {
  if (!inherits(x, "TP") && !inherits(x, "fitMod")) {
    stop(deparse(substitute(x)), " should be an object of class TP or fitMod.\n")
  }
  timePointsX <- attr(x, which = "timePoints")
  if (is.character(timePoints)) {
    if (!all(timePoints %in% timePointsX[["timePoint"]])) {
      stop("All timePoints should be in ", deparse(substitute(x)), ".\n")
    }
  } else if (is.numeric(timePoints)) {
    if (!all(timePoints %in% timePointsX[["timeNumber"]])) {
      stop("All timePoints should be in ", deparse(substitute(x)), ".\n")
    }
    timePoints <- timePointsX[timePoints, "timePoint"]
  } else {
    stop("timePoints should be a character or numeric vector.\n")
  }
  return(timePoints)
}


