#' @export
createTD <- function(dat,
                     genotype,
                     timePoint,
                     plotId = NULL,
                     rowNum = NULL,
                     colNum = NULL,
                     rowId = rowNum,
                     colId = colNum,
                     checkId = NULL) {
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
  for (param in c(genotype, timePoint, plotId,
                  rowId, colId, rowNum, colNum, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !hasName(dat, param))) {
      stop(paste(deparse(param), "has to be NULL or a column in dat.\n"))
    }
  }
  ## Create list of reserved column names for renaming columns.
  renameCols <- c("genotype", "timePoint", "plotId", "rowId", "colId",
                  "rowNum", "colNum", "checkId")
  ## First rename duplicate colums and add duplicated columns to dat
  renameFrom <- as.character(sapply(X = renameCols, FUN = function(x) {
    get(x)
  }))
  ## Create a data.frame with renamed cols to add to TD as an attribute.
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
  dat[["timePoint"]] <- lubridate::ymd_hms(dat[["timePoint"]])
  ## Add time
  dat[["time"]] <- as.factor(dat[["timePoint"]])
  ## Convert columns to factor if neccessary.
  factorCols <-  c("genotype", "plotId", "rowId", "colId", "checkId")
  for (factorCol in factorCols) {
    if (hasName(dat, factorCol)) {
      dat[[factorCol]] <- as.factor(dat[[factorCol]])
    }
  }
  ## Convert columns to numeric if neccessary.
  numCols <- c("rowCoord", "colCoord")
  for (numCol in numCols) {
    if (hasName(dat, numCol) && !is.numeric(dat[cols == numCol])) {
      dat[cols == numCol] <- as.numeric(dat[, cols == numCol])
    }
  }
  listData <- split(x = dat, f = dat$time)
  ## Set meta for all trials in dat.
  for (tr in names(listData)) {
    ## Add a list of columns that have been renamed as attribute to TD.
    attr(x = listData[[tr]], which = "renamedCols") <-
      if (nrow(renamed) > 0) renamed else NULL
  }
  TD <- structure(listData,
                  class = c("TD", "list"))
  return(TD)
}
