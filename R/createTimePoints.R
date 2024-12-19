#' Create an object of class TP
#'
#' Convert a data.frame to an object of class TP (Time Points).
#' The function converts a data.frame to an object of class TP in the following
#' steps:
#' \itemize{
#' \item{Quality control on the input data. For example, warnings will be given
#' when more than 50% of observations are missing for a plant.}
#' \item{Rename columns to default column names used by the functions in the
#' statgenHTP package. For example, the column in the data containing
#' variety/accession/genotype is renamed to “genotype.” Original column names
#' are stored as an attribute of the individual data.frames in the TP object.}
#' \item{Convert column types to the default column types. For example, the
#' column “genotype” is converted to a factor and “rowNum” to a numeric column.}
#' \item{Convert the column containing time into time format. If needed, the
#' time format can be provided in \code{timeFormat}. For example, with a
#' date/time input of the form “day/month/year hour:minute”, use
#' "%d/%m/%Y %H:%M". For a full list of abbreviations see the R package
#' strptime. When the input time is a numeric value, the function will
#' convert it to time from 01-01-1970.}
#' \item{If \code{addCheck = TRUE}, the genotypes listed in
#' \code{checkGenotypes} are reference genotypes (or check). It will add a
#' column check with a value "noCheck" for the genotypes that are not in
#' \code{checkGenotypes} and the name of the genotypes for the
#' \code{checkGenotypes}. A column genoCheck is also added with the names of
#' the genotypes that are not in \code{checkGenotypes} and \code{NA} for the
#' \code{checkGenotypes}. These columns are necessary for fitting models on
#' data that includes check genotypes, e.g. reference genotypes that are
#' highly replicated or in case of augmented design.}
#' \item{Split the data into separate data.frames by time point. A TP object is
#' a list of data.frames where each data.frame contains the data for a single
#' time point. If there is only one time point the output will be a list with
#' only one item.}
#' \item{Add a data.frame with columns timeNumber and timePoint as attribute
#' “timePoints” to the TP object. This data.frame can be used for referencing
#' time points by a unique number.}
#' }
#' Note that \code{plotId} needs to be a unique identifier for a plot or a
#' plant. It cannot occur more than once per time point.
#'
#' @param dat A data.frame.
#' @param experimentName A character string, the name of the experiment. Stored
#' with the data and used in default plot titles.
#' @param genotype A character string indicating the column in dat containing
#' the genotypes.
#' @param timePoint A character string indicating the column in dat containing
#' the time points.
#' @param timeFormat A character string indicating the input format of the time
#' points. E.g. for a date/time input of the form day/month/year hour:minute,
#' use "%d/%m/%Y %H:%M". For a full list of abbreviations see
#' \code{\link[base]{strptime}}. If \code{NULL}, a best guess is done based on
#' the input.
#' @param plotId A character string indicating the column in dat containing
#' the plotId. This has to be a unique identifier per plot/plant per time point.
#' @param repId A character string indicating the column in dat containing
#' the replicates.
#' @param rowNum A character string indicating the column in dat containing
#' the row number of the plot.
#' @param colNum A character string indicating the column in dat containing
#' the column number of the plot.
#' @param addCheck Should a column check be added to the output? If \code{TRUE},
#' checkGenotypes cannot be \code{NULL}.
#' @param checkGenotypes A character vector containing the genotypes used as
#' checks in the experiment.
#'
#' @returns An object of class \code{TP}. A list with, per time point in the
#' input, a data.frame containing the data for that time point. A data.frame
#' with columns timeNumber and timePoint is added as attribute timePoints to
#' the data. This data.frame can be used for referencing timePoints by their
#' number.
#'
#' @examples ## Create a TP object containing the data from the Phenovator.
#' phenoTP <- createTimePoints(dat = PhenovatorDat1,
#'                             experimentName = "Phenovator",
#'                             genotype = "Genotype",
#'                             timePoint = "timepoints",
#'                             repId = "Replicate",
#'                             plotId = "pos",
#'                             rowNum = "y", colNum = "x",
#'                             addCheck = TRUE,
#'                             checkGenotypes = c("check1", "check2",
#'                                                "check3","check4"))
#' summary(phenoTP)
#'
#' @family functions for data preparation
#'
#' @export
createTimePoints <- function(dat,
                             experimentName,
                             genotype,
                             timePoint,
                             timeFormat = NULL,
                             plotId,
                             repId = NULL,
                             rowNum = NULL,
                             colNum = NULL,
                             addCheck = FALSE,
                             checkGenotypes = NULL) {
  ## Checks.
  if (missing(dat) || !is.data.frame(dat)) {
    stop("dat has to be a data.frame.\n")
  }
  if (missing(experimentName) || !is.character(experimentName) ||
      length(experimentName) > 1) {
    stop("experimentName should be a character string of length one.\n")
  }
  ## Convert input to data.frame. This needs to be done to be able to handle
  ## tibbles and possibly other data structures in the future.
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  cols <- colnames(dat)
  ## Check that all columns to be renamed are columns in dat.
  for (param in c(genotype, timePoint, plotId, repId, rowNum, colNum)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !hasName(dat, param))) {
      stop(substitute(param), " has to be NULL or a column in dat.\n")
    }
  }
  ## Check uniqueness of plotId.
  if (max(table(dat[[plotId]], dat[[timePoint]])) > 1) {
    stop("plotId has to be unique within each time point.\n")
  }
  ## Check for plotIds that have a limited amount of observations.
  plotTab <- table(dat[[plotId]])
  nTimePoints <- length(unique(dat[[timePoint]]))
  plotLimObs <- names(plotTab[plotTab < 0.5 * nTimePoints])
  if (length(plotLimObs) > 5) {
    warning("More than 5 plotIds have observations for less than 50% of ",
            "the time points. The first 5 are printed, to see them all run
              attr(..., 'plotLimObs') on the output\n",
            paste(plotLimObs[1:5], collapse = ", "), "\n", call. = FALSE)
  } else if (length(plotLimObs) > 0) {
    warning("The following plotIds have observations for less than 50% of ",
            "the time points:\n", paste(plotLimObs, collapse = ", "), "\n",
            call. = FALSE)
  }
  ## Check presence of check genotypes.
  if (addCheck) {
    if (is.null(checkGenotypes) || !is.character(checkGenotypes)) {
      stop("checkGenotypes should be a character vector.\n")
    }
    if (!all(checkGenotypes %in% dat[[genotype]])) {
      stop("All checkGenotypes should be genotypes in dat.\n")
    }
  }
  ## Set rowId and colId to rowNum and colNum.
  ## rowId and colId will be factor with identical information as numeric
  ## rowNum and colNum needed by both SpATS and asreml in modeling.
  rowId <- rowNum
  colId <- colNum
  ## Create list of reserved column names for renaming columns.
  renameCols <- c("genotype", "timePoint", "plotId", "repId", "rowId", "colId",
                  "rowNum", "colNum")
  ## Get the value associated with the column to be renamed.
  ## These values are derived from the input parameters.
  renameFrom <- as.character(sapply(X = renameCols, FUN = function(x) {
    get(x)
  }))
  ## Create a data.frame with renamed cols to add to TP as an attribute.
  renamed <- data.frame(orig = renameFrom[renameFrom != "NULL"],
                        new = renameCols[renameFrom != "NULL"],
                        stringsAsFactors = FALSE)
  ## First rename duplicate colums and add duplicated columns to dat.
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
  for (i in seq_along(renameCols)) {
    cols[cols == renameFrom[i]] <- renameCols[i]
  }
  colnames(dat) <- cols
  ## Convert timepoint to nice format.
  if (is.numeric(dat[["timePoint"]])) {
    warning("timePoint is a numeric column. It will be converted using ",
            "1970-01-01 as 0.\n")
    dateTime <- try(lubridate::as_datetime(dat[["timePoint"]]))
  } else {
    ## When using format an the output is of class POSIXlt.
    ## This doesn't work with split later on and therefore is converted to
    ## POSIXct.
    dateTime <- try({
      as.POSIXct(lubridate::as_datetime(as.character(dat[["timePoint"]]),
                                        format = timeFormat))
    })
  }
  if (inherits(dateTime, "try-error") || all(is.na(dateTime))) {
    stop("Error when converting timePoints to Date format. Please check the
         input format of the timePoints and use the timeFormat argument.\n")
  } else {
    dat[["timePoint"]] <- dateTime
  }
  ## Add check genotype.
  if (addCheck) {
    ## Add column check with genotype for check genotypes and noCheck for
    ## all other genotypes. To be used as covariate in model.
    dat[["check"]] <- ifelse(dat[["genotype"]] %in% checkGenotypes,
                             as.character(dat[["genotype"]]), "noCheck")
    ## Add column genoCheck with NA for check genotypes and genotype for
    ## all other genotypes. Preserve original column so models without check
    ## can still be fitted.
    dat[["genoCheck"]] <- ifelse(dat[["genotype"]] %in% checkGenotypes,
                                 NA, as.character(dat[["genotype"]]))
  }
  ## Convert columns to factor if necessary.
  factorCols <-  c("genotype", "plotId", "repId", "rowId", "colId", "check",
                   "genoCheck")
  for (factorCol in factorCols) {
    if (hasName(dat, factorCol) && !is.factor(dat[[factorCol]])) {
      dat[[factorCol]] <- factor(dat[[factorCol]],
                                 levels = sort(unique(dat[[factorCol]])))
    }
  }
  ## Convert columns to numeric if necessary.
  numCols <- c("rowNum", "colNum")
  for (numCol in numCols) {
    if (hasName(dat, numCol) && !is.numeric(dat[cols == numCol])) {
      dat[cols == numCol] <- as.numeric(dat[, cols == numCol])
    }
  }
  if (all(hasName(dat, c("rowNum", "colNum")))) {
    ## Check that row column combinations are unique within time points.
    rowColTab <- table(dat[["timePoint"]], dat[["rowNum"]], dat[["colNum"]])
    if (any(rowColTab > 1)) {
      stop("Combinations of row and column coordinates should be unique ",
           "within time points.\n")
    }
  }
  listData <- split(x = dat, f = dat[["timePoint"]])
  ## Set meta for all timePoints in dat.
  for (tr in names(listData)) {
    ## Add a list of columns that have been renamed as attribute to TP.
    attr(x = listData[[tr]], which = "renamedCols") <-
      if (nrow(renamed) > 0) renamed else NULL
  }
  ## Create a data.frame with coding of time points for reference.
  timePoints <- data.frame(timeNumber = seq_along(listData),
                           timePoint = names(listData),
                           stringsAsFactors = FALSE)
  TP <- structure(listData,
                  experimentName = experimentName,
                  timePoints = timePoints,
                  checkGenotypes = checkGenotypes,
                  plotLimObs = plotLimObs,
                  class = c("TP", "list"))
  return(TP)
}

#' Remove time points from an object of class TP
#'
#' Function for removing selected time points from an object of class TP.
#'
#' @param TP An object of class TP.
#' @param timePoints A character or numeric vector indicating the time points
#' to be removed. When using a character string to reference a time point, the
#' value has to be an exact match to one of the existing timePoints. When using
#' a number it will be matched by its number ("timeNumber") in the timePoints
#' attribute of the TP object.
#'
#' @returns An object of class TP, the input with the selected time points
#' removed.
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
#'                                                "check3","check4"))
#' ## Remove the first and last time point from the TP object.
#' phenoTPNew <- removeTimePoints(phenoTP,
#'                                timePoints = c(1, 73))
#'
#' ## Compare by looking at summaries.
#' summary(phenoTP)
#' summary(phenoTPNew)
#'
#' @family functions for data preparation
#'
#' @export
removeTimePoints <- function(TP,
                             timePoints) {
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  timePoints <- chkTimePoints(TP, timePoints)
  TPOut <- TP[setdiff(names(TP), timePoints)]
  return(TPOut)
}

#' Summary function for TP objects
#'
#' Function for creating a short summary of the contents of a \code{TP} object.
#' The summary consists of the name of the experiment, the number of time
#' points, the first and last time point and the genotypes defined as checks.
#'
#' @param object An object of class TP.
#' @param ... Ignored.
#'
#' @returns No return value, a summary is printed.
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
#'                                                "check3","check4"))
#' ## Create a summary.
#' summary(phenoTP)
#'
#' @family functions for data preparation
#'
#' @export
summary.TP <- function(object,
                       ...)  {
  experimentName <- attr(x = object, which = "experimentName")
  noTP <- length(object)
  firstTP <- names(object)[1]
  lastTP <- names(object)[noTP]
  checkGenotypes <- attr(x = object, which = "checkGenotypes")

  cat(deparse(substitute(object)), " contains data for experiment ",
      experimentName, ".\n\n", sep = "")
  cat("It contains", noTP, "time points.\n")
  cat("First time point:", firstTP, "\n")
  cat("Last time point:", lastTP, "\n\n")
  if (!is.null(checkGenotypes)) {
    cat("The following genotypes are defined as check genotypes: " ,
        paste(checkGenotypes, collapse = ", "), ".\n", sep = "")
  } else {
    cat("No check genotypes are defined.\n")
  }
}

#' Plot function for class TP
#'
#' Plotting function for objects of class TP. Plots the layout of the platform
#' for different time points within the TP object. Also a boxplot can be made
#' for selected traits and time points and a plot of correlations
#' between time points. Finally the raw data can be displayed per
#' genotype. A detailed description and optional extra parameters for the
#' different plots are given in the sections below.
#'
#' @section Layout Plot:
#' Plots the layout of the platform for selected time points (all available time
#' points by  default). This plot can only be made for time points that contain
#' both row (\code{rowNum}) and column (\code{colNum}) information. If either
#' one of those is missing the timePoint is skipped with a warning.
#' If replicates (\code{repId}) are available, a black line is plotted between
#' diffent replicates. Missing plots are indicated in white. This can either be
#' single plots in a time point or complete missing columns or rows.\cr
#' Extra parameter options:
#' \describe{
#' \item{showGeno}{Should individual genotypes be labeled in the plot?
#' Defaults to \code{FALSE}}
#' \item{highlight}{A character vector of genotypes to be highlighted in the
#' plot.}
#' }
#'
#' @section Box Plot:
#' Creates a boxplot per selected trait grouped by time point (all available
#' time points by  default). Extra parameter options:
#' \describe{
#' \item{groupBy}{A character string indicating a column in \code{TP} by which
#' the boxes in the plot should be grouped.
#' By default the boxes are grouped per time point.}
#' \item{colorBy}{A character string indicating a column in \code{TP} by which
#' the boxes are colored. Coloring will be done within the groups indicated by
#' the \code{groupBy} parameter, e.g. per replicate within each time point
#' using \code{repId}.}
#' \item{orderBy}{A character string indicating the way the boxes should be
#' ordered. Either "alphabetic" for alphabetical ordering of the groups,
#' "ascending" for ordering by ascending mean, or "descending" for ordering by
#' descending mean. By default boxes are ordered alphabetically.}
#' }
#'
#' @section Correlation Plot:
#' Draws a heatmap of correlations of raw data between time points per selected
#' trait for selected time points (all available time points by default).
#'
#' @section Raw data plot:
#' Create a plot of the raw data of the selected trait over time for selected
#' time points (all available time points by default). Plots are grouped by
#' genotype, or by genotype x treatment when the \code{geno.decomp} option is
#' specified. By default, all the genotypes will be plotted which might take
#' time and memory when the output is not saved in a file (see parameter
#' \code{outFile}). Extra parameter options:
#'\describe{
#' \item{genotypes}{A character vector indicating the genotypes to be plotted.}
#' \item{geno.decomp}{A character vector indicating the grouping of the
#' genotypes to be plotted.}
#' \item{plotLine}{Should the data be displayed as lines? Default is FALSE.}
#' }
#'
#' @param x An object of class TP.
#' @param ... Extra plot options. Described per plotType in their respective
#' section.
#' @param plotType A single character string indicating which plot should be
#' made. See the sections below for a detailed explanation of the plots.
#' @param timePoints A character or numeric vector indicating the time points
#' to be plotted. When using a character string to reference a time point, the
#' value has to be an exact match to one of the existing timePoints. When using
#' a number it will be matched by its number ("timeNumber") in the timePoints
#' attribute of the TP object.
#' @param traits A character vector indicating the traits to be plotted.
#' If \code{plotType} = "layout" only a single trait may be plotted. For the
#' other plotTypes, providing multiple traits will create multiple plots.
#' @param title A character string used as title for the plot. If \code{NULL} a
#' default title is added to the plot depending on \code{plotType}.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a (list of) ggplot object(s) is invisibly returned. Ignored if
#' \code{outFile} is specified.
#' @param outFile A character string indicating the .pdf file to which the
#' plots should be written. If \code{NULL}, no file is written.
#' @param outFileOpts A named list of extra options for the pdf outfile, e.g.
#' width and height. See \code{\link[grDevices]{pdf}} for all possible options.
#'
#' @returns Depending on the plot type, either a ggplot object or a list of
#' ggplot objects is invisibly returned.
#'
#' @examples
#' \donttest{
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
#'
#' ## Plot the layout for the third time point with the check genotypes
#' ## highlighted
#' plot(phenoTP,
#'      plotType = "layout",
#'      timePoints = 3,
#'      highlight = c("check1", "check2", "check3", "check4"))
#'
#' ## Create a boxplot for "EffpsII" with 5 time points and boxes colored
#' ## by "repId" within time point.
#' plot(phenoTP,
#'      plotType = "box",
#'      traits = "EffpsII",
#'      timePoints = 1:5,
#'      colorBy = "repId")
#'
#' ## Create a correlation plot for "EffpsII" for a selection of time points.
#' plot(phenoTP,
#'      plotType = "cor",
#'      traits = "EffpsII",
#'      timePoints = seq(from=1, to=73, by=5))
#'
#' ## Plot the raw data of four genotypes for the trait "EffpsII":
#' plot(phenoTP,
#'      traits = "EffpsII",
#'      plotType = "raw",
#'      genotypes = c("G001","G002","check1","check2"))
#' }
#'
#' @family functions for data preparation
#'
#' @export
plot.TP <- function(x,
                    ...,
                    plotType = c("layout", "box", "cor", "raw"),
                    timePoints = names(x),
                    title = NULL,
                    traits = NULL,
                    output = TRUE,
                    outFile = NULL,
                    outFileOpts = NULL) {
  ## Checks.
  timePoints <- chkTimePoints(x, timePoints)
  plotType <- match.arg(plotType)
  experimentName <- attr(x = x, which = "experimentName")
  dotArgs <- list(...)
  if (!is.null(outFile)) {
    chkFile(outFile, fileType = "pdf")
    output <- TRUE
    outFileOpts <- c(list(file = outFile), outFileOpts)
    on.exit(dev.off(), add = TRUE)
    do.call(pdf, args = outFileOpts)
  }
  if (plotType == "layout") {
    showGeno <- isTRUE(dotArgs$showGeno)
    highlight <- dotArgs$highlight
    if (!is.null(highlight) && !is.character(highlight)) {
      stop("highlight should be a character vector.\n")
    }
    if (!is.null(traits) && (!is.character(traits) || length(traits) > 1)) {
      stop("traits should be NULL or a character string.\n")
    }
    plotTitle <- title
    p <- setNames(vector(mode = "list", length = length(timePoints)),
                  timePoints)
    for (timePoint in timePoints) {
      if (is.null(title)) {
        plotTitle <- paste(experimentName, "-", timePoint)
      }
      tpDat <- x[[timePoint]]
      if (!is.null(traits) && !hasName(x = tpDat, name = traits)) {
        stop("traits should be a column in TP.\n")
      }
      if (!chkRowCol(tpDat)) next
      if (length(highlight) > 0) {
        tpDat$highlight. <- ifelse(tpDat$genotype %in% highlight,
                                   as.character(tpDat$genotype), NA)
      }
      ## Compute aspect for proper depiction of field size.
      ## Plots are assumed to be square.
      aspect <- length(unique(tpDat$colNum)) / length(unique(tpDat$rowNum))
      plotRep <- hasName(x = tpDat, name = "repId")
      ## Create data for lines between replicates.
      if (plotRep) {
        repBord <- calcPlotBorders(tpDat = tpDat, bordVar = "repId")
      }
      ## Create base plot.
      pTp <- ggplot2::ggplot(data = tpDat,
                             ggplot2::aes(x = .data[["colNum"]],
                                          y = .data[["rowNum"]])) +
        ggplot2::coord_fixed(ratio = aspect,
                             xlim = range(tpDat$colNum) + c(-0.5, 0.5),
                             ylim = range(tpDat$rowNum) + c(-0.5, 0.5),
                             clip = "off") +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ## Move ticks to edge of the plot.
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(),
                                    expand = c(0, 0)) +
        ggplot2::ggtitle(plotTitle)
      if (sum(!is.na(tpDat$highlight.)) > 0) {
        ## Genotypes to be highlighted get a color.
        ## Everything else the NA color.
        pTp <- pTp +
          ggplot2::geom_tile(ggplot2::aes(fill = .data[["highlight."]]),
                             color = "grey75") +
          ggplot2::labs(fill = "Highlighted") +
          ## Remove NA from scale.
          ggplot2::scale_fill_discrete(na.translate = FALSE)
      } else if (!is.null(traits)) {
        pTp <- pTp +
          ggplot2::geom_tile(ggplot2::aes(fill = .data[[traits]]),
                             color = "grey75") +
          ## Adjust plot colors.
          ggplot2::scale_fill_gradientn(colors = topo.colors(100))
      } else {
        ## No hightlights so just a single fill color.
        pTp <- pTp +
          ggplot2::geom_tile(fill = "white", color = "grey75")
      }
      if (showGeno) {
        ## Add names of genotypes to the center of the tiles.
        pTp <- pTp +
          ggplot2::geom_text(ggplot2::aes(label = .data[["genotype"]]),
                             size = 2, check_overlap = TRUE)
      }
      if (plotRep) {
        ## Add lines for replicates.
        pTp <- pTp +
          ## Add vertical lines as segment.
          ## adding/subtracting 0.5 assures plotting at the borders of
          ## the tiles.
          ggplot2::geom_segment(ggplot2::aes(x = .data[["x"]] - 0.5,
                                             xend = .data[["x"]] - 0.5,
                                             y = .data[["y"]] - 0.5,
                                             yend = .data[["y"]] + 0.5,
                                             linetype = "replicates"),
                                data = repBord$vertW, linewidth = 1) +
          ggplot2::geom_segment(ggplot2::aes(x = .data[["x"]] - 0.5,
                                             xend = .data[["x"]] + 0.5,
                                             y = .data[["y"]] - 0.5,
                                             yend = .data[["y"]] - 0.5),
                                data = repBord$horW, linewidth = 1)
      }
      if (plotRep) {
        pTp <- pTp +
          ## Add a legend entry for replicates and subBlocks.
          ggplot2::scale_linetype_manual("replicates",
                                         values = c("replicates" = "solid"),
                                         name = ggplot2::element_blank()) +
          ggplot2::guides(linetype = ggplot2::guide_legend(override.aes =
                                                             list(size = 1)))
      }
      p[[timePoint]] <- pTp
      if (output) {
        plot(pTp)
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
    plotTitle <- title
    for (trait in traits) {
      ## Create a single data.frame from x with only columns timePoint, trait
      ## and genotype. Genotype is needed to be able to display hovering info.
      ## timePoints where trait is not measured/available are removed by setting
      ## them to NULL.
      if (is.null(title)) {
        plotTitle <- paste(experimentName, "-", trait)
      }
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
        warning(trait, " isn't a column in any of the timePoints.\n",
                "Plot skipped.\n", call. = FALSE)
        next
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
      ## Colorby is ignored in plot if it is not a factor.
      if (!is.null(colorBy) && !is.factor(plotDat[colorBy])) {
        plotDat[colorBy] <- factor(plotDat[[colorBy]])
      }
      ## Create boxplot.
      ## Back ticks around variable names are needed to handle spaces in names.
      if (is.null(colorBy)) {
        colorByTr <- ".colorBy"
        plotDat[[".colorBy"]] <- factor(1)
      } else {
        colorByTr <- colorBy
      }
      pTp <- ggplot2::ggplot(plotDat,
                             ggplot2::aes(x = .data[[xVar]],
                                          y = .data[[trait]],
                                          fill = .data[[colorByTr]])) +
        ggplot2::geom_boxplot(na.rm = TRUE,
                              show.legend = colorByTr != ".colorBy") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5,
                                                           hjust = 1)) +
        ggplot2::labs(title = plotTitle, x = xVar, y = trait)
      p[[trait]] <- pTp
      if (output) {
        plot(pTp)
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
    plotTitle <- title
    for (trait in traits) {
      if (is.null(title)) {
        plotTitle <- paste(experimentName, "- Correlations of timepoints for",
                           trait)
      }
      ## Create a single data.frame from x with only columns timePoint and trait
      ## timePoints where trait is not measured/available are removed by setting
      ## them to NULL.
      plotDat <- Reduce(f = rbind,
                        x = lapply(X = x[timePoints],
                                   FUN = function(timePoint) {
                                     if (!hasName(x = timePoint, name = trait)) {
                                       NULL
                                     } else {
                                       timePoint[c("plotId", "timePoint",
                                                   trait)]
                                     }
                                   }))
      if (is.null(plotDat)) {
        warning(trait, " isn't a column in any of the timePoints.\n",
                "Plot skipped.\n", call. = FALSE)
        next
      }
      ## Create table with the value of trait per plotId per timePoint.
      plotTab <- tapply(plotDat[[trait]],
                        INDEX = list(plotDat[["plotId"]],
                                     plotDat[["timePoint"]]), FUN = I)
      ## Create a correlation matrix.
      corMat <- cor(plotTab, use = "pairwise.complete.obs")
      ## Remove rows and columns with only NA.
      corKeep <- sapply(X = seq_len(ncol(corMat)), FUN = function(i) {
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
      pTp <- ggplot2::ggplot(data = meltedCorMat,
                             ggplot2::aes(x = .data[["Var1"]],
                                          y = .data[["Var2"]],
                                          fill = .data[["value"]])) +
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
        ggplot2::labs(title = plotTitle, x = "", y = "", color = "") +
        ggplot2::guides(size = "none") +
        ## Fix coordinates to get a square sized plot.
        ggplot2::coord_fixed()
      p[[trait]] <- pTp
      if (output) {
        plot(pTp)
      }
    }
  } else if (plotType == "raw") {
    if (is.null(traits) || !is.character(traits)) {
      stop("traits should be a character vector.\n")
    }
    genotypes <- dotArgs$genotypes
    plotLine  <- isTRUE(dotArgs$plotLine)
    geno.decomp <- dotArgs$geno.decomp
    if (!is.null(geno.decomp) && !all(sapply(X = x, FUN = function(timePoint) {
      hasName(x = timePoint, name = geno.decomp)
    }))) {
      stop("geno.decomp should be a column in TP.\n")
    }
    p <- setNames(vector(mode = "list", length = length(traits)), traits)
    plotTitle <- title
    for (trait in traits) {
      if (is.null(title)) {
        plotTitle <- paste(experimentName, "-", trait, "- raw data")
      }
      ## Create a single data.frame from x with only columns timePoint and
      ## trait. timePoints where trait is not measured/available are removed
      ## by setting them to NULL.
      plotDat <- Reduce(f = rbind,
                        x = lapply(X = x[timePoints],
                                   FUN = function(timePoint) {
                                     if (!hasName(x = timePoint, name = trait)) {
                                       NULL
                                     } else {
                                       timePoint[c("genotype", "timePoint",
                                                   "plotId", trait,
                                                   geno.decomp)]
                                     }
                                   }))
      if (is.null(plotDat)) {
        warning(trait, " isn't a column in any of the timePoints.\n",
                "Plot skipped.\n", call. = FALSE)
        next
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
                  title = plotTitle, yLab = trait, output = output,
                  plotLine = plotLine)
    }
  }
  invisible(p)
}

#' Function for extracting for objects of class TP that keeps class and
#' attributes.
#'
#' @param x An object of class TP.
#' @param i An index specifying the element to extract or replace.
#' @param ... Ignored.
#'
#' @noRd
#' @export
`[.TP` <- function(x,
                   i,
                   ...) {
  timePoints <- chkTimePoints(x, i)
  timePointsX <- attr(x, which = "timePoints")
  timePointsR <- timePointsX[timePointsX[["timePoint"]] %in% timePoints, ]
  if (nrow(timePointsR) > 0) {
    class(x) <- "list"
    r <- x[timePointsR[["timePoint"]]]
    attr(r, "timePoints") <- timePointsR
    attr(r, "experimentName") <- attr(x, "experimentName")
    attr(r, "class") <- c("TP", "list")
  } else {
    r <- NULL
  }
  return(r)
}

#' Coerce TP object to data.frame
#'
#' Function for converting an object of class TP to a data.frame.
#'
#' @param x An object of class TP.
#' @param ... Ignored.
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
#' ## Convert phenoTP to data.frame.
#' phenoDat <- as.data.frame(phenoTP)
#'
#' @returns A data.frame containing the data.frames for all time points in the
#' TP object bound together.
#'
#' @family functions for data preparation
#'
#' @export
as.data.frame.TP <- function (x,
                              ...) {
  res <- dfBind(x)
  return(res)
}
