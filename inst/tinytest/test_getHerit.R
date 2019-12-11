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
## of getHerit function.
testFitMod1 <- fitModels(testTP, trait = "t1", quiet = TRUE)
testFitMod2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)

### Check input.

expect_error(getHerit("fitMod"),
             "fitMod should be an object of class fitMod")

herit1 <- getHerit(testFitMod1)
herit2 <- getHerit(testFitMod2)

## Check output structure.

expect_true(inherits(herit1, "data.frame"))
expect_equal(dim(herit1), c(5, 3))
expect_equal(colnames(herit1), c("timeNumber", "timePoint", "h2"))

# Geno.decomp -> 1 column per level of geno.decomp.
expect_equal(dim(herit2), c(5, 4))
expect_equal(colnames(herit2), c("timeNumber", "timePoint", "1", "2"))

# Compare results.
expect_equal(herit1[["h2"]], c(0, 0.21, 0.01, 0.91, 0.42))
expect_equal(herit2[["1"]], c(0, 0.7, 0.02, 0.49, 0))
expect_equal(herit2[["2"]], c(0.3, 0.77, 0.58, 0.78, 0.71))

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getHerit(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
heritOut <- getHerit(testFitMod1, outFile = tmpFile)
heritIn <- read.csv(tmpFile)

# Ignore timePoint
# It is imported as factor, but is a date in the original data.
expect_equal(heritOut[, -2], heritIn[, -2])
