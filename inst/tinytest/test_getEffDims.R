### Test getCorrected.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x", addCheck = TRUE,
                           checkGenotypes = "check1")

## Fit some models with different options that might influence output
## of getEffDims function.
testFitMod1 <- fitModels(testTP, trait = "t1", quiet = TRUE)
testFitMod2 <- fitModels(testTP, trait = "t1", useRepId = TRUE, quiet = TRUE)

### Check input.

expect_error(getEffDims("fitMod"),
             "fitMod should be an object of class fitMod")

effDims1 <- getEffDims(testFitMod1)
effDims2 <- getEffDims(testFitMod2)
effRats1 <- getEffDims(testFitMod1, EDType = "ratio")
effRats2 <- getEffDims(testFitMod2, EDType = "ratio")

## Check output structure.

expect_inherits(effDims1, "data.frame")
expect_equal(dim(effDims1), c(5, 10))
expect_equal(colnames(effDims1), c("timeNumber", "timePoint", "colId", "rowId",
                                   "fCol", "fRow", "fColRow", "colfRow",
                                   "fColfRow", "surface"))

# useRepId should give different colnames.
expect_equal(colnames(effDims2), c("timeNumber", "timePoint", "repId:colId",
                                   "repId:rowId", "fCol", "fRow", "fColRow",
                                   "colfRow", "fColfRow", "surface"))

# Compare results.

expect_equal_to_reference(effDims1, file = "effDims1Comp")
expect_equal_to_reference(effDims2, file = "effDims2Comp")
expect_equal_to_reference(effRats1, file = "effRats1Comp")
expect_equal_to_reference(effRats2, file = "effRats2Comp")

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getEffDims(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
effDimsOut <- getEffDims(testFitMod1, outFile = tmpFile)
effDimsIn <- read.csv(tmpFile)

# Ignore timePoint
# It is imported as factor, but is a date in the original data.
expect_equal(effDimsOut[, -2], effDimsIn[, -2])

## Test that eff dims cannot be computed for asreml models.

if (at_home()) {
  testFitModAs1 <- fitModels(testTP, trait = "t1", engine = "asreml",
                             quiet = TRUE)
  expect_error(getEffDims(testFitModAs1),
               "Models in testFitModAs1 should be fitted using SpATS")
}
