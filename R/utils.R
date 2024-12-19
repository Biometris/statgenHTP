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
calcPlotBorders <- function(tpDat,
                            bordVar) {
  yMin <- min(tpDat$rowNum)
  yMax <- max(tpDat$rowNum)
  xMin <- min(tpDat$colNum)
  xMax <- max(tpDat$colNum)
  ## Create matrix containing replicates.
  ## First create an empty matrix containing all row/column values
  ## between min and max to assure complete missing rows/columns
  ## are added.
  M <- matrix(nrow = yMax - yMin + 1, ncol = xMax - xMin + 1,
              dimnames = list(yMin:yMax, xMin:xMax))
  for (i in seq_len(nrow(tpDat))) {
    M[as.character(tpDat[i, "rowNum"]),
      as.character(tpDat[i, "colNum"])] <- tpDat[i, bordVar]
  }
  ## Create an imputed version of M for plotting borders around NA values.
  MImp <- M
  MImp[is.na(MImp)] <- nlevels(tpDat[[bordVar]]) + 1
  has.breaks <- function(x) {
    ncol(x) == 2 & nrow(x) > 0
  }
  ## Create a data.frame with positions where the value of rep in the
  ## data changes in vertical direction.
  vertW <- do.call(rbind.data.frame,
                   Filter(f = has.breaks, x = Map(function(i, x) {
                     cbind(y = i, x = which(diff(c(0, x, 0)) != 0))
                   }, seq_len(nrow(MImp)), split(MImp, seq_len(nrow(MImp))))))
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
                  }, seq_len(ncol(MImp)), as.data.frame(MImp))))
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
  ## timePoint (rhs), timeNumber (rhs) and trait (value var).
  castCols <- setdiff(colnames(dat), c("timePoint", "timeNumber",  trait))
  rhsCols <- intersect(c("timePoint", "timeNumber"), colnames(dat))
  castForm <- formula(paste(paste(castCols, collapse = "+"), "~ " ,
                            paste(rhsCols, collapse = "+")))
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
                        output = TRUE,
                        plotLine = FALSE) {
  ## Create plot.
  p <- ggplot2::ggplot(baseDat, ggplot2::aes(x = .data[[xVal]],
                                             y = .data[[yVal]])) +
    ## Format the time scale to Month + day.
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0, "cm"),
                   panel.border = ggplot2::element_rect(color = "black",
                                                        fill = "transparent"),
                   strip.background = ggplot2::element_rect(color = "black",
                                                            fill = "bisque"),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 25, vjust = 1,
                                                       hjust = 1)) +
    ggplot2::labs(title = title, x = xLab, y = yLab)
  nTp <- length(unique(baseDat[["timePoint"]]))
  if (nTp < 5) {
    p <- p + ggplot2::scale_x_datetime(breaks = unique(baseDat[["timePoint"]]),
                                       labels = scales::date_format("%B %d"))
  } else {
    ## Format the time scale to Month + day.
    p <- p + ggplot2::scale_x_datetime(breaks = prettier(n = 3),
                                       labels = scales::date_format("%B %d"))
  }

  if (!plotLine || length(unique(baseDat[[xVal]])) == 1) {
    ## Multiple time points in data. Display a line.
    p <- p + ggplot2::geom_point(ggplot2::aes(group = .data[[groupVal]],
                                              color = .data[[colVal]]),
                                 show.legend = FALSE, na.rm = TRUE, size = 1)
  } else {
    ## Only one time point. Makes geom_line crash. Display as point.
    p <- p + ggplot2::geom_line(ggplot2::aes(group = .data[[groupVal]],
                                             color = .data[[colVal]]),
                                show.legend = FALSE, na.rm = TRUE)
  }
  if (!is.null(overlayDat)) {
    ## Add a second data set as overlay over the first plot.
    if (!plotLine || length(unique(baseDat[[xVal]])) == 1) {
      ## Multiple time points in data. Display a line.
      p <- p + ggplot2::geom_point(ggplot2::aes(x = .data[[xVal]],
                                                y = .data[[yValOverlay]]),
                                   data = overlayDat, color = "black", size = 1,
                                   show.legend = FALSE, na.rm = TRUE)
    } else {
      ## Only one time point. Makes geom_line crash. Display as point.
      p <- p +
        ggplot2::geom_line(ggplot2::aes(x = .data[[xVal]],
                                        y = .data[[yValOverlay]]),
                           data = overlayDat, color = "black", linewidth = 1,
                           show.legend = FALSE, na.rm = TRUE)
    }
  }
  ## Calculate the total number of plots.
  ## facetVal is a vector so use interaction to get the number of levels.
  nPlots <- nlevels(interaction(baseDat[facetVal], drop = TRUE))
  ## 25 Plots per page.
  nPag <- ceiling(nPlots / 25)
  if (nPlots >= 25) {
    ## More than 25 plots.
    ## For identical layout on all pages use 5 x 5 plots throughout.
    rowPag <- colPag <- rep(x = 5, times = nPag)
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
      ggforce::facet_wrap_paginate(facets = facetVal,
                                   nrow = rowPag[i], ncol = colPag[i],
                                   labeller = ggplot2::label_wrap_gen(multi_line = FALSE),
                                   page = i)
    if (output) {
      suppressMessages(plot(pPag[[i]]))
    }
  }
  invisible(pPag)
}

#' @noRd
#' @keywords internal
chkFile <- function(outFile,
                    fileType = "csv") {
  if (!is.character(outFile) || length(outFile) > 1 ||
      tools::file_ext(outFile) != fileType) {
    stop("outFile should be a single character string ending in .",
         fileType, ".\n")
  }
  if (file.access(dirname(outFile), 2)) {
    stop("No permission to write to ", outFile, ".\n")
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
#' @returns A character vector containing timePoints.
#'
#' @noRd
#' @keywords internal
chkTimePoints <- function(x,
                          timePoints) {
  if (!inherits(x, "TP") && !inherits(x, "fitMod")) {
    stop(deparse(substitute(x)),
         " should be an object of class TP or fitMod.\n")
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
    timePoints <- timePointsX[timePointsX[["timeNumber"]] %in% timePoints,
                              "timePoint"]
  } else {
    stop("timePoints should be a character or numeric vector.\n")
  }
  return(timePoints)
}

#' Helper function for checking structure of smooth arguments to
#' fitSplineHDM.
#'
#' @noRd
#' @keywords internal
chkSmooth <- function(x) {
  xName <- deparse(substitute(x))
  if (!is.list(x) || length(x) != 3 ||
      !(setequal(names(x), c("nseg", "bdeg", "pord")))) {
    stop(xName, " should be a named list of length 3.\n")
  }
  if (!is.numeric(x$nseg) || length(x$nseg) > 1 || x$nseg < 1) {
    stop("nseg in ", xName, " should be a positive numerical value.\n")
  }
  if (!is.numeric(x$bdeg) || length(x$bdeg) > 1 || x$bdeg < 1) {
    stop("bdeg in ", xName, " should be a positive numerical value.\n")
  }
  if (!is.numeric(x$pord) || length(x$pord) > 1 || x$pord < 1) {
    stop("pord in ", xName, " should be a positive numerical value.\n")
  }
}

#' Helper function for minimal plot theme.
#'
#' @noRd
#' @keywords internal
plotTheme <- function() {
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 legend.key = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(color = "black"),
                 plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Helper function for creating a decent looking time scale on the x-axis of
#' time series plots.
#'
#' @noRd
#' @keywords internal
prettier <- function(n = 3) {
  function(x) {
    ## Get first time point.
    intStart <- x[1]
    ## Compute interval length in seconds.
    sec <- x[2] - x[1]
    if (n == 1) {
      ## Just one time point. Label in the middle.
      intStart + sec / 2
    } else if (n == 2) {
      ## Two time points, labels at 1/5 and 4/5.
      intStart + sec / 5 * c(1, 4)
    } else {
      ## Three time points, labels at 1/8, 1/2 and 7/8.
      intStart + sec / 8 * c(1, 4, 7)
    }
  }
}

#' Helper function for checking if asreml 4.0 or higher is installed and licence
#' activated.
#'
#' @importFrom utils packageVersion
#'
#' @noRd
#' @keywords internal
checkAsreml <- function() {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("No valid installation of asreml found.\n")
  }
  asremlVersion <- packageVersion("asreml")
  if (asremlVersion[1] < "4") {
    stop("asreml version 4.0 or higher is requiered.\n")
  }
  licenceStatus <- asreml::asreml.license.status(quiet = TRUE)
  if (licenceStatus$status != 0) {
    stop("Error checking asreml licence status:\n", licenceStatus$statusMessage)
  }
  invisible(TRUE)
}

#' Helper function for checking row and column information
#'
#' Check that row and column information is available in data.frame.
#' Columns should be present and cannot contain NA values.
#'
#' @noRd
#' @keywords internal
chkRowCol <- function(dat) {
  timePoint <- dat[["timePoint"]][1]
  rowCol <- TRUE
  if (!"rowNum" %in% colnames(dat)) {
    rowCol <- FALSE
    warning("rowNum should be a column in ", timePoint, ".\n",
            "Plot skipped.\n", call. = FALSE)
  } else if (sum(is.na(dat[["rowNum"]])) > 0) {
    rowCol <- FALSE
    warning("rowNum contains missing values for ", timePoint, ".\n",
            "Plot skipped.\n", call. = FALSE)
  }
  if (!"colNum" %in% colnames(dat)) {
    rowCol <- FALSE
    warning("colNum should be a column in ", timePoint, ".\n",
            "Plot skipped.\n", call. = FALSE)
  } else if (sum(is.na(dat[["colNum"]])) > 0) {
    rowCol <- FALSE
    warning("colNum contains missing values for ", timePoint, ".\n",
            "Plot skipped.\n", call. = FALSE)
  }
  return(rowCol)
}

#' Helper function for row binding data.frames with different columns.
#'
#' @param dfList A list of data.frames.
#'
#' @noRd
#' @keywords internal
dfBind <- function(dfList) {
  ## Remove empty data.frames from dfList
  for (i in rev(seq_along(dfList))) {
    if (nrow(dfList[[i]]) == 0) {
      dfList[[i]] <- NULL
    }
  }
  if (length(dfList) == 0) {
    return(data.frame())
  }
  ## Get variable names from all data.frames.
  allNms <- unique(unlist(lapply(dfList, names)))
  ## rbind all data.frames setting values for missing columns to NA.
  do.call(rbind,
          c(lapply(X = dfList, FUN = function(x) {
            nwDat <- sapply(X = setdiff(allNms, names(x)), FUN = function(y) {
              NA
            })
            data.frame(c(x, nwDat), check.names = FALSE,
                       stringsAsFactors = FALSE)
          }), make.row.names = FALSE)
  )
}

#' Helper function for computing angles
#'
#' @noRd
#' @keywords internal
angle <- function(M) {
  dotProd <- M[1, ] %*% M[2, ]
  norm1 <- norm(M[1, ], type = "2")
  norm2 <- norm(M[2, ], type = "2")
  theta <- acos(dotProd / (norm1 * norm2))
  return(as.numeric(theta) * 180 / pi)
}

#' Count valid observations per time point for a given trait
#'
#' Count valid observations per time point for a given trait.
#'
#' @param TP An object of class TP.
#' @param trait A character string indicating the trait for which valid
#' observations should be counted.
#'
#' @returns A named numerical vector with he number of valid observations per
#' time point .
#'
#' @examples
#' ## Create a TP object containing the data from the Phenovator.
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#' ## Count valid observations for EffpsII per time point.
#' validPheno <- countValid(phenoTP, trait = "EffpsII")
#' head(validPheno)
#'
#' @export
countValid <- function(TP,
                       trait) {
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length one.\n")
  }
  sapply(X = TP, FUN = function(timepoint) {
    sum(!is.na(timepoint[[trait]]))
  })
}

#' Count valid observations per plotId for a given trait
#'
#' Count valid observations per plotId for a given trait.
#'
#' @param TP An object of class TP.
#' @param trait A character string indicating the trait for which valid
#' observations should be counted.
#' @param plotIds A character vector indicating the plotIds for which valid
#' observations should be checked. If \code{NULL} valid observations are
#' counted for all plotIds in TP.
#'
#' @returns A named numerical vector with he number of valid observations per
#' plotId.
#'
#' @examples
#' ## Create a TP object containing the data from the Phenovator.
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3", "check4"))
#' ## Count valid observations for EffpsII for a subset of plots.
#' countValidPlot(phenoTP,
#'                trait = "EffpsII",
#'                plotIds = c("c12r22", "c24r41", "c14r32"))
#'
#' @export
countValidPlot <- function(TP,
                           trait,
                           plotIds = NULL) {
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  TPTot <- as.data.frame(TP)
  if (!is.null(plotIds) && !is.character(plotIds)) {
    stop("plotIds should be NULL or a character vector.\n")
  }
  if (is.null(plotIds)) {
    plotIds <- as.character(unique(TPTot[["plotId"]]))
  } else {
    if (!all(plotIds %in% TPTot[["plotId"]])) {
      stop("All plotIds should be in TP.\n")
    }
  }
  if (!is.character(trait) || length(trait) > 1) {
    stop("trait should be a character string of length one.\n")
  }
  if (!hasName(x = TPTot, name = trait)) {
    stop(trait, " should be a column in TP.\n")
  }
  plotIdsValid <- as.character(TPTot[!is.na(TPTot[[trait]]), "plotId"])
  sapply(X = plotIds, FUN = function(plotId) {
    sum(plotIdsValid == plotId)
  })
}

#' @noRd
#' @keywords internal
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  N <- length(lmat)
  k <- dim(lmat[[1]])[1]
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M, M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M - 1L), nrow = k)[, rep(seq_len(N), each = k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive = FALSE, use.names = FALSE)))
}
