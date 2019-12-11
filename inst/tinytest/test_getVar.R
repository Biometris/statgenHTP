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
## of getVar function.
testFitMod1 <- fitModels(testTP, trait = "t1", quiet = TRUE)
testFitMod2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)
testFitMod3 <- fitModels(testTP, trait = "t1", useCheck = TRUE, quiet = TRUE)
testFitMod4 <- fitModels(testTP, trait = "t1", useRepId = TRUE, quiet = TRUE)

### Check input.

expect_error(getVar("fitMod"),
             "fitMod should be an object of class fitMod")

var1 <- getVar(testFitMod1)
var2 <- getVar(testFitMod2)
var3 <- getVar(testFitMod3)
var4 <- getVar(testFitMod4)

## Check output structure.

expect_true(inherits(var1, "data.frame"))
expect_equal(dim(var1), c(5, 6))
expect_equal(colnames(var1), c("timeNumber", "timePoint", "varGen", "varRes",
                               "varCol", "varRow"))

# Read expected results.
var1Orig <- read.csv("var1")
var2Orig <- read.csv("var2")
var3Orig <- read.csv("var3")
var4Orig <- read.csv("var4")

# Compare results.
expect_equal(var1[, 3:6], var1Orig[, 3:6])
#expect_equal(var2[, 3:6], var2Orig[, 3:6])
expect_equal(var3[, 3:6], var3Orig[, 3:6])
expect_equal(var4[, 3:6], var4Orig[, 3:6])

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getVar(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
varOut <- getVar(testFitMod1, outFile = tmpFile)
varIn <- read.csv(tmpFile)

# Ignore timePoint
# It is imported as factor, but is a date in the original data.
expect_equal(varOut[, -2], varIn[, -2])

