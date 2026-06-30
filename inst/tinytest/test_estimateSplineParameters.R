### Test estimateSplineParameters

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

## Fit a spline on the corrected values.
splineRes <- fitSpline(inDat = corr, trait = "t1_corr")

## Fit a spline with plotId absent.
corr2 <- corr[, !colnames(corr) == "plotId"]
splineRes2 <- fitSpline(inDat = corr2, trait = "t1_corr")

### Check extracting parameters from fitted splines.

## Check that general checks in estimateSplineParameters function correctly.
expect_error(estimateSplineParameters(splineRes, what = "a"),
             "should be one of")
expect_error(estimateSplineParameters(splineRes, what = "pp"),
             "A percentile should be give as pN, with N between 0 and 100")
expect_error(estimateSplineParameters(splineRes, what = "min", genotypes = "a"),
             "genotypes should be a character vector of genotypes in predDat")
expect_error(estimateSplineParameters(splineRes, what = "min", plotIds = "a"),
             "plotIds should be a character vector of plotIds in predDat")
expect_error(estimateSplineParameters(splineRes, what = "min",
                                      genotypes = "G12", plotIds = "c13r2"),
             "At least one valid combination of genotype and plotId should be selected")

expect_silent(est1 <- estimateSplineParameters(splineRes, what = "min"))

expect_inherits(est1, "splineEst")

expect_equal_to_reference(est1, "splineEst", tolerance = 1e-6)

## Check that options timeMin and timeMax function correctly.
# Get first and last timePoint from data.
startTime <- min(corr[["timePoint"]])
endTime <- max(corr[["timePoint"]])

expect_error(estimateSplineParameters(splineRes, what = "min",
                                      timeMin = startTime - 1),
             "timeMin should be within the time interval in the data")
expect_error(estimateSplineParameters(splineRes, what = "min",
                                      timeMax = endTime + 1),
             "timeMax should be within the time interval in the data")
expect_error(estimateSplineParameters(splineRes, what = "min",
                                      timeMin = endTime, timeMax = startTime),
             "timeMax should be larger than timeMin")

est2 <- estimateSplineParameters(splineRes, what = "min",
                                 timeMin = startTime, timeMax = endTime)
expect_equal(est1, est2)

## Check that estimates are made correctly when plotId is absent.
est3 <- estimateSplineParameters(splineRes2, what = "min")

expect_equal(ncol(est3), 4)

## Check that estimates are made correctly when spline where fitted on timeNumber.
expect_silent(estimateSplineParameters(splineRes2, what = "min"))

## Check that estimates for multiple parameters are made correctly.
expect_silent(est4 <- estimateSplineParameters(splineRes,
                                               what = c("min", "max")))

expect_equal(est4[["min_predictions"]], est1[["min_predictions"]])


### HDM splines

## Create TP object.
testTP2 <- createTimePoints(dat = testDat, experimentName = "testExp",
                            genotype = "Genotype", timePoint = "timepoints",
                            repId = "Replicate", plotId = "pos",
                            rowNum = "y", colNum = "x")

## Fit model.
testFitMod2 <- fitModels(testTP2, trait = "t1", geno.decomp = "repId",
                         extraFixedFactors = "Basin", quiet = TRUE)

## Get corrected values.
corr2 <- getCorrected(testFitMod2)

## Move all check1s to Basin 1 for a consistent population structure.
corr2[corr2[["genotype"]] == "check1", "Basin"] <- 1

## Fit a HDM spline on the corrected values.
splineRes3 <- fitSplineHDM(inDat = corr2, trait = "t1_corr",
                           pop = "Basin", trace = FALSE)

## Check that estimates work correctly for HDM splines.
est5 <- estimateSplineParameters(splineRes3,
                                 what = c("min", "max", "AUC", "p10"),
                                 AUCScale = "hour")

expect_equal_to_reference(est5, "splineEstHDM", tolerance = 1e-6)


## Check that option fitLevel works correctly.

est6 <- estimateSplineParameters(splineRes3, what = "mean", fitLevel = "plot")
est7 <- estimateSplineParameters(splineRes3, what = "mean", fitLevel = "genoDev")

expect_equal(nrow(est6), 25)
expect_equal(nrow(est7), 22)


## Check that plotting of estimates works correctly.

## Create temporary file for output.
tmpFile <- tempfile(fileext = ".pdf")

p1 <- plot(est1, plotType = "box", outFile = tmpFile)
expect_inherits(p1, "list")
expect_inherits(p1[[1]], "ggplot")

p2 <- plot(est5, plotType = "hist", what = "p10", outFile = tmpFile)
expect_inherits(p2, "list")
expect_inherits(p2[[1]], "ggplot")

