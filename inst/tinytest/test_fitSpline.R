### Test fitSpline

Sys.setlocale("LC_COLLATE", "C")

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)

## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", rowNum = "y", colNum = "x")

## Fit model.
testFitMod <- fitModels(testTP, trait = "t1", quiet = TRUE)

## Get corrected values.
corr <- getCorrected(testFitMod)

## Check that general checks in fitSpline function correctly.
expect_error(fitSpline(1, trait = "t1"),
             "inDat should be a data.frame")
expect_error(fitSpline(corr, trait = 1),
             "trait should be a character string of length 1")
expect_error(fitSpline(corr, trait = "t2"),
             paste("inDat should at least contain the following columns:",
                   "genotype, t2, plotId, timePoint"))

## Fit a spline on the corrected values.
expect_warning(splineRes <- fitSpline(inDat = corr,
                                      trait = "t1_corr"),
               "for less than the minimum number of time points, which is 4")

## Check that general structure of the output is correct.
expect_inherits(splineRes, c("HTPSpline", "list"))
expect_equal(length(splineRes), 2)
expect_equal(names(splineRes), c("coefDat", "predDat"))

expect_inherits(splineRes[["coefDat"]], "data.frame")
expect_inherits(splineRes[["predDat"]], "data.frame")

## get coefDat and predDat.
coefDat <- splineRes[["coefDat"]]
predDat <- splineRes[["predDat"]]

## Check that full coefDat and predDat results are correct.
expect_equal_to_reference(coefDat, file = "coefDat", tolerance = 1e-6)
expect_equal_to_reference(predDat, file = "predDat", tolerance = 1e-4)

## Check that option timeNumber functions correctly.
expect_error(fitSpline(corr, trait = "t1_corr", useTimeNumber = TRUE),
             "timeNumber should be a character string of length 1")
expect_error(fitSpline(corr, trait = "t1_corr", useTimeNumber = TRUE,
                       timeNumber = "genotype"),
             "timeNumber should be a numerical column")

expect_warning(splineRes2 <- fitSpline(inDat = corr,
                                       trait = "t1_corr",
                                       useTimeNumber = TRUE,
                                       timeNumber = "timeNumber"),
               "for less than the minimum number of time points, which is 4")

## Check that option genotypes functions correctly.
expect_error(fitSpline(inDat = corr, trait = "t1_corr", genotypes = 1),
             "genotypes should be a character vector of genotypes in inDat")
expect_error(fitSpline(inDat = corr, trait = "t1_corr", genotypes = "G0"),
             "genotypes should be a character vector of genotypes in inDat")
splineRes3 <- fitSpline(inDat = corr, trait = "t1_corr", genotypes = "G12")
expect_true(all(splineRes3[["coefDat"]][["genotype"]] == "G12"))
expect_true(all(splineRes3[["predDat"]][["genotype"]] == "G12"))

## Check that option plotIds functions correctly.
expect_error(fitSpline(inDat = corr, trait = "t1_corr", plotIds = 1),
             "plotIds should be a character vector of plotIds in inDat")
expect_error(fitSpline(inDat = corr, trait = "t1_corr", plotIds = "p1"),
             "plotIds should be a character vector of plotIds in inDat")
expect_error(fitSpline(inDat = corr, trait = "t1_corr", genotypes = "G12",
                        plotIds = "c13r2"),
             "At least one valid combination of genotype and plotId should be selected")
splineRes4 <- fitSpline(inDat = corr, trait = "t1_corr", plotIds = "c13r2")
expect_true(all(splineRes4[["coefDat"]][["plotId"]] == "c13r2"))
expect_true(all(splineRes4[["predDat"]][["plotId"]] == "c13r2"))

## Check that option knots functions correctly.
expect_error(fitSpline(inDat = corr, trait = "t1_corr", knots = "a"),
             "knots should be a positive numerical value")
expect_error(fitSpline(inDat = corr, trait = "t1_corr", knots = 3),
             "Number of knots should be at least 4 for proper spline fitting")
expect_silent(fitSpline(inDat = corr, trait = "t1_corr",
                        genotypes = "G12", knots = 20))

## Check that option minNoTP functions correctly.
expect_error(fitSpline(inDat = corr, trait = "t1_corr", minNoTP = "a"),
             "minNoTP should be a numerical value")
expect_error(fitSpline(inDat = corr, trait = "t1_corr", minNoTP = 50),
             "minNoTP should be a number bewtween 0 and 5")
expect_silent(fitSpline(inDat = corr, trait = "t1_corr",
                        genotypes = "G12", minNoTP = 5))

# Add some extra NA to corr.
corr2 <- corr
corr2[["t1_corr"]][1:10] <- NA

expect_warning(splineRes5 <- fitSpline(inDat = corr2, trait = "t1_corr",
                                       minNoTP = 5),
               "More than 5 plotIds have observations for less than the minimum")
expect_equal(attr(splineRes5, which = "plotLimObs"),
             c("c10r1", "c10r2", "c10r3", "c10r5", "c11r2", "c12r3", "c12r4",
               "c13r1", "c13r5", "c14r2", "c14r3"))

## Check that splines are fitted correctly when plotId is absent.
corr3 <- corr[, !colnames(corr) == "plotId"]
splineRes6 <- fitSpline(inDat = corr3, trait = "t1_corr")

expect_equal(ncol(splineRes6[["coefDat"]]), 3)
expect_equal(ncol(splineRes6[["predDat"]]), 6)

### Check plotting of fitSpline results.

# Create temporary file for exporting plots.
tmpFile <- tempfile(fileext = ".pdf")

## Check that general checks in plot function correctly.
expect_error(plot(splineRes, plotType = "a"),
             "should be one of")
expect_error(plot(splineRes, genotypes = "a"),
             "genotypes should be a character vector of genotypes in predDat")
expect_error(plot(splineRes, plotIds = "a"),
             "plotIds should be a character vector of plotIds in predDat")
expect_error(plot(splineRes, genotypes = "G12", plotIds = "c13r2"),
             "At least one valid combination of genotype and plotId should be selected")

## Check that general output structure is correct.
expect_silent(p <- plot(splineRes))
expect_inherits(p, "list")
expect_equal(length(p), 1)
expect_inherits(p[[1]], "ggplot")

## Check that option title functions correctly.
p1 <- plot(splineRes, plotIds = "c13r2", title = "bla")
expect_equal(p1[[1]]$labels$title, "bla")

# Changing plotType should change the default title.
expect_equal(p[[1]]$labels$title, "Corrected data and P-spline prediction")
p2 <- plot(splineRes, plotIds = "c13r2", plotType = "derivatives")
expect_equal(p2[[1]]$labels$title, "P-spline first derivatives")
p3 <- plot(splineRes, plotIds = "c13r2", plotType = "derivatives2")
expect_equal(p3[[1]]$labels$title, "P-spline second derivatives")

## Check that plotting fitted splines without plotId functions correctly.
expect_silent(plot(splineRes6, genotypes = "G12"))

## Check that plotting to pdf functions correctly.
expect_silent(plot(splineRes, plotIds = "c13r2", output = TRUE,
                   outFile = tmpFile))

## Remove tmpFile
unlink(tmpFile)


