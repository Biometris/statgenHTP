## Test plot.fitMod

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
                           colNum = "x", addCheck = TRUE,
                           checkGenotypes = "check1")

## Create fitMod object.
testFitMod <- fitModels(testTP, trait = "t1", quiet = TRUE)
## Create another fitMod object for testing with geno.decomp.
testFitMod2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)
## Create another fitMod object for testing with check.
#testFitMod3 <- fitModels(testTP, trait = "t1", useCheck = TRUE, quiet = TRUE)

if (at_home() && FALSE) {
  ## Create fitMods for additional testing with asreml.
  testFitModAs <- fitModels(testTP, trait = "t1", engine = "asreml",
                            quiet = TRUE)
  testFitModAs2 <- fitModels(testTP, trait = "t1", engine = "asreml",
                             spatial = TRUE, quiet = TRUE)
}

## Create a temporary outfile for writing plots.
tmpFile <- tempfile(fileext = ".pdf")

## Check that general checks in plot.TP function correctly.
expect_error(plot(testFitMod, plotType = "test"), "should be one of")
expect_error(plot(testFitMod, title = 1), "title should be NULL or a character")

### Check rawPred plot.

expect_error(plot(testFitMod, plotType = "rawPred", genotypes = 1),
             "genotypes should be NULL or a character vector")
expect_error(plot(testFitMod, plotType = "rawPred", genotypes = "g1"),
             "All genotypes should be in testFitMod")

expect_silent(p0 <- plot(testFitMod, plotType = "rawPred", outFile = tmpFile))
expect_true(inherits(p0, "list"))
expect_equal(length(p0), 1)
expect_true(inherits(p0[[1]], "ggplot"))

geoms0 <- sapply(p0[[1]]$layers, function(x) class(x$geom)[1])
expect_equal(geoms0, c("GeomPoint", "GeomPoint"))

## Check that rawPred plots function correctly for single timePoints.
expect_silent(p1 <- plot(testFitMod[1], plotType = "rawPred",
                         outFile = tmpFile))
geoms1 <- sapply(p1[[1]]$layers, function(x) class(x$geom)[1])
expect_equal(geoms1, c("GeomPoint", "GeomPoint"))

## Check option genotypes in rawpred plots.
expect_silent(p2 <- plot(testFitMod, plotType = "rawPred", genotypes = "G12",
                         outFile = tmpFile))
nCol <- ggplot2::ggplot_build(p2[[1]])$layout$facet$params$ncol
nRow <- ggplot2::ggplot_build(p2[[1]])$layout$facet$params$nrow

## Grid should be modified since only one genotype left.
expect_equal(nRow, 1)
expect_equal(nCol, 1)

## Check that rawPred plot functions for model with geno.decomp.
expect_silent(p3 <- plot(testFitMod2, plotType = "rawPred", outFile = tmpFile))

## Check that rawPred plot functions for model with check.
# expect_silent(p4 <- plot(testFitMod3, plotType = "rawPred", outFile = tmpFile))
# expect_silent(p5 <- plot(testFitMod3, plotType = "rawPred", plotChecks = TRUE,
#                          outFile = tmpFile))
# expect_equal(nrow(p4[[1]]$data), 105)
# expect_equal(nrow(p5[[1]]$data), 125)

if (at_home() && FALSE) {
  ## Check that rawPred plot functions for asreml.
  expect_silent(plot(testFitModAs, plotType = "rawPred"))
}

### Check corrPred plot.

expect_error(plot(testFitMod, plotType = "corrPred", genotypes = 1),
             "genotypes should be NULL or a character vector")
expect_error(plot(testFitMod, plotType = "corrPred", genotypes = "g1"),
             "All genotypes should be in testFitMod")

expect_silent(p0 <- plot(testFitMod, plotType = "corrPred", outFile = tmpFile))
expect_true(inherits(p0, "list"))
expect_equal(length(p0), 1)
expect_true(inherits(p0[[1]], "ggplot"))

geoms0 <- sapply(p0[[1]]$layers, function(x) class(x$geom)[1])
expect_equal(geoms0, c("GeomPoint", "GeomPoint"))

## Check that corrPred plots function correctly for single timePoints.
expect_silent(p1 <- plot(testFitMod[1], plotType = "corrPred",
                         outFile = tmpFile))
geoms1 <- sapply(p1[[1]]$layers, function(x) class(x$geom)[1])
expect_equal(geoms1, c("GeomPoint", "GeomPoint"))

## Check option genotypes in corrPred plots.
expect_silent(p2 <- plot(testFitMod, plotType = "corrPred", genotypes = "G12",
                         outFile = tmpFile))
nCol <- ggplot2::ggplot_build(p2[[1]])$layout$facet$params$ncol
nRow <- ggplot2::ggplot_build(p2[[1]])$layout$facet$params$nrow

## Grid should be modified since only one genotype left.
expect_equal(nRow, 1)
expect_equal(nCol, 1)

## Check that corrPred plot functions for model with geno.decomp.
expect_silent(p3 <- plot(testFitMod2, plotType = "corrPred", outFile = tmpFile))

## Check that rawPred plot functions for model with check.
# expect_silent(p4 <- plot(testFitMod3, plotType = "corrPred", outFile = tmpFile))
# expect_silent(p5 <- plot(testFitMod3, plotType = "corrPred", plotChecks = TRUE,
#                          outFile = tmpFile))
# expect_equal(nrow(p4[[1]]$data), 105)
# expect_equal(nrow(p5[[1]]$data), 125)

if (at_home() && FALSE) {
  ## Check that rawPred plot functions for asreml.
  expect_silent(plot(testFitModAs, plotType = "corrPred"))
}

### Check heritability plot.

expect_silent(p0 <- plot(testFitMod, plotType = "herit", outFile = tmpFile))
expect_true(inherits(p0, "ggplot"))

## Output should be a combination of points and lines.
geoms0 <- sapply(p0$layers, function(x) class(x$geom)[1])
expect_equal(geoms0, c("GeomPoint", "GeomLine"))

## If there is only one timepoint output should be only points.
expect_silent(p1 <- plot(testFitMod[1], plotType = "herit", outFile = tmpFile))
geoms1 <- sapply(p1$layers, function(x) class(x$geom)[1])
expect_equal(geoms1, c("GeomPoint"))

## Check option yLim in heritability plot.
expect_silent(p2 <- plot(testFitMod, plotType = "herit", yLim = c(0, 1),
                         outFile = tmpFile))
expect_equal(as.list(p2$scales$get_scales("y"))$limits, c(0, 1))

## Check that heritability plot functions for model with geno.decomp.
expect_silent(p3 <- plot(testFitMod2, plotType = "herit", outFile = tmpFile))

### Check effective dimensions plot.

expect_silent(p0 <- plot(testFitMod, plotType = "effDim", outFile = tmpFile))
expect_true(inherits(p0, "ggplot"))

## Output should be a combination of points and lines.
geoms0 <- sapply(p0$layers, function(x) class(x$geom)[1])
expect_equal(geoms0, c("GeomPoint", "GeomLine"))

## If there is only one timepoint output should be only points.
expect_silent(p1 <- plot(testFitMod[1], plotType = "effDim", outFile = tmpFile))
geoms1 <- sapply(p1$layers, function(x) class(x$geom)[1])
expect_equal(geoms1, c("GeomPoint"))

## Check option yLim in effDim plot.

expect_silent(p2 <- plot(testFitMod, plotType = "effDim", yLim = c(0, 100),
                         outFile = tmpFile))
expect_equal(as.list(p2$scales$get_scales("y"))$limits, c(0, 100))

## Check option EDType in effDim plot.

expect_error(plot(testFitMod, plotType = "effDim", EDType = "ED"),
             "should be one of")
expect_silent(p3 <- plot(testFitMod, plotType = "effDim", EDType = "ratio",
                         outFile = tmpFile))
expect_equal(as.list(p3$scales$get_scales("y"))$limits, c(0, 0.777319493997852))

## Check option which in effDim plot.
expect_silent(p4 <- plot(testFitMod, plotType = "effDim", whichED = "colId",
                         outFile = tmpFile))

## Check that effDim plot functions for model with geno.decomp.
expect_silent(p5 <- plot(testFitMod2, plotType = "effDim", outFile = tmpFile))

## Check that plotting is not possible with models fitted with asreml.
if (at_home() && FALSE) {
  expect_error(plot(testFitModAs, plotType = "effDim"),
               "only be plotted for models fitted with SpATS")
}

### Check variance plot.

expect_silent(p0 <- plot(testFitMod, plotType = "variance", outFile = tmpFile))
expect_true(inherits(p0, "ggplot"))

## Output should be a combination of points and lines.
geoms0 <- sapply(p0$layers, function(x) class(x$geom)[1])
expect_equal(geoms0, c("GeomPoint", "GeomLine"))

## If there is only one timepoint output should be only points.
expect_silent(p1 <- plot(testFitMod[1], plotType = "variance",
                         outFile = tmpFile))
geoms1 <- sapply(p1$layers, function(x) class(x$geom)[1])
expect_equal(geoms1, c("GeomPoint"))

## Check option yLim in variance plots.

expect_silent(p2 <- plot(testFitMod, plotType = "variance", yLim = c(0, 1e-3),
                         outFile = tmpFile))
expect_equal(as.list(p2$scales$get_scales("y"))$limits, c(0, 1e-3))

## Check that effDim plot functions for model with geno.decomp.
expect_silent(p3 <- plot(testFitMod2, plotType = "variance", outFile = tmpFile))

### Check spatial plots.

expect_silent(p0 <- plot(testFitMod, plotType = "spatial", outFile = tmpFile))
expect_true(inherits(p0, "list"))
expect_equal(length(p0), 5)
expect_true(inherits(p0[[1]], "list"))
expect_equal(length(p0[[1]]), 6)
expect_true(inherits(p0[[1]][[1]], "ggplot"))

## Check option spaTrend in spatial plots.

expect_error(plot(testFitMod, plotType = "spatial", spaTrend = "sTr"),
             "should be one of")
expect_silent(p1 <- plot(testFitMod, plotType = "spatial",
                         spaTrend = "percentage", outFile = tmpFile))

## Check that spatial plot functions for model with check.
# expect_silent(plot(testFitMod3, plotType = "spatial", outFile = tmpFile))

## Check that effDim plot functions for model with geno.decomp.
expect_silent(p2 <- plot(testFitMod2, plotType = "spatial", outFile = tmpFile))

if (at_home() && FALSE) {
  ## Check that spatial plots cannot be made for asreml when spatial = FALSE.
  expect_error(plot(testFitModAs, plotType = "spatial"),
               "when setting spatial = TRUE when fitting the asreml models")

  p3 <- plot(testFitModAs2, plotType = "spatial")
  expect_equal(length(p3), 5)
}

### Check time lapse plots.

## Create a second outfile with .gif extension
tmpFile2 <- tempfile(fileext = ".gif")

expect_silent(p0 <- plot(testFitMod, plotType = "timeLapse",
                         outFile = tmpFile2))

## Check that effDim plot functions for model with geno.decomp.
expect_silent(p1 <- plot(testFitMod2, plotType = "timeLapse",
                         outFile = tmpFile2))

## Remove tmpFiles
unlink(tmpFile)
unlink(tmpFile2)

