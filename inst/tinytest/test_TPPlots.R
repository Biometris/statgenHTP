## Test plot.TP

## Testing the exact plot output is difficult but since also the ggplot
## objects on which the plots are based are invisibly returned at least some
## checking can be done.

## Some problems only occur when the actual plotting is done.
## Therefore setting an outfile for all checks.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x")

## Check that general checks in plot.TP function correctly.
expect_error(plot(testTP, plotType = "test"), "should be one of")

## Check that an outfile is actually created.
tmpFile <- tempfile(fileext = ".pdf")
plot(testTP, outFile = tmpFile)
# Only checking that tmpFile has some content.
expect_true(file.size(tmpFile) > 0)

### Check layout plot.

p0 <- plot(testTP, plotType = "layout", outFile = tmpFile)
expect_true(inherits(p0, "list"))
expect_equal(length(p0), 5)
expect_true(inherits(p0[[1]], "ggplot"))

## Create copy for testing row/column condition.
testTP1a <- testTP1b <- testTP
testTP1a[[1]][["colNum"]] <- NULL
expect_warning(plot(testTP1a, plotType = "layout"), "colNum should be")
testTP1b[[1]][["rowNum"]] <- NULL
expect_warning(plot(testTP1b, plotType = "layout"), "rowNum should be")

## Check option showGeno in layout plot.

p1 <- plot(testTP, plotType = "layout", showGeno = TRUE, outFile = tmpFile)
## Difference with default plot p0 should be the extra GeomText layer.
geoms0 <- sapply(p0[[1]]$layers, function(x) class(x$geom)[1])
geoms1 <- sapply(p1[[1]]$layers, function(x) class(x$geom)[1])
expect_equal(setdiff(geoms1, geoms0), "GeomText")


## Check option highlight in layout plot.

expect_error(plot(testTP, plotType = "layout", highlight = 1),
             "highlight should be a character vector")
p1 <- plot(testTP, plotType = "layout", highlight = "check1", outFile = tmpFile)
geoms1 <- sapply(p1[[1]]$layers, function(x) class(x$geom)[1])
## Two plots should be highlighted as defined in variable highlight..
expect_equal(as.character(p1[[1]]$layers[geoms1 == "GeomTile"][[1]]$mapping),
             "~highlight.")
expect_equal(sum(!is.na(p1[[1]]$data$highlight.)), 2)

### Check box plot.

expect_error(plot(testTP, plotType = "box", traits = 1),
             "traits should be a character vector")
expect_warning(plot(testTP, plotType = "box", traits = "trait"),
               "trait isn't a column in any of the timePoints")
p <- plot(testTP, plotType = "box", traits = "t1", outFile = tmpFile)
expect_true(inherits(p, "list"))
expect_equal(length(p), 1)
expect_true(inherits(p[[1]], "ggplot"))

## Check option groupBy for box plot.

expect_error(plot(testTP, plotType = "box", traits = "t1", groupBy = 1),
             "groupBy should be a single character string")
expect_error(plot(testTP, plotType = "box", traits = "t1", groupBy = "grp"),
             "groupBy should be a column in TP")
p <- plot(testTP, plotType = "box", traits = "t1", groupBy = "repId",
          outFile = tmpFile)
expect_true("~repId" %in% as.character(p$t1$mapping))

## Check option colorBy for box plot.

expect_error(plot(testTP, plotType = "box", traits = "t1", colorBy = 1),
             "colorBy should be a single character string")
expect_error(plot(testTP, plotType = "box", traits = "t1", colorBy = "grp"),
             "colorBy should be a column in TP")
p <- plot(testTP, plotType = "box", traits = "t1", colorBy = "repId",
          outFile = tmpFile)
expect_true(all(c("~repId", "~timePoint") %in% as.character(p$t1$mapping)))

## Check option orderBy for box plot.

p0 <- plot(testTP, plotType = "box", traits = "t1", outFile = tmpFile)
p1 <- plot(testTP, plotType = "box", traits = "t1",
           orderBy = "ascending", outFile = tmpFile)
p2 <- plot(testTP, plotType = "box", traits = "t1", orderBy = "descending",
           outFile = tmpFile)
## This basically only checks that releveling took place.
expect_equal(setdiff(names(p1$t1$plot_env), names(p0$t1$plot_env)),
             "levNw")
expect_equal(setdiff(names(p2$t1$plot_env), names(p0$t1$plot_env)),
             "levNw")

### Check correlation plot.

expect_error(plot(testTP[1], plotType = "cor", traits = "trait"),
             "At least two timePoints requiered for a correlation plot")
expect_error(plot(testTP, plotType = "cor", traits = 1),
             "traits should be a character vector")
expect_warning(plot(testTP, plotType = "cor", traits = "trait"),
               "trait isn't a column in any of the timePoints")
p <- plot(testTP, plotType = "cor", traits = "t1", outFile = tmpFile)
expect_true(inherits(p, "list"))
expect_equal(length(p), 1)
expect_true(inherits(p[[1]], "ggplot"))

### Check raw plot.

expect_error(plot(testTP, plotType = "raw", traits = 1),
             "traits should be a character vector")
expect_warning(plot(testTP, plotType = "raw", traits = "trait"),
               "trait isn't a column in any of the timePoints")
p <- plot(testTP, plotType = "raw", traits = "t1", outFile = tmpFile)
expect_true(inherits(p, "list"))
expect_equal(length(p), 1)

## Check raw plot for single time point.
expect_silent(plot(testTP[1], plotType = "raw", traits = "t1",
                   outFile = tmpFile))

## Check option genotypes for raw plot.

expect_silent(plot(testTP, plotType = "raw", traits = "t1",
                   genotypes = "G12", outFile = tmpFile))

## Check option geno.decomp for raw plot.

expect_error(plot(testTP, plotType = "raw", traits = "t1",
                  geno.decomp = "gd"),
             "geno.decomp should be a column in TP")
expect_silent(plot(testTP, plotType = "raw", traits = "t1",
                   geno.decomp = "Basin", outFile = tmpFile))

## Remove tmpFile
unlink(tmpFile)

