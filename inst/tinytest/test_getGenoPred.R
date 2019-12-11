### Test getGenoPred.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x", addCheck = TRUE,
                           checkGenotypes = "check1")

## Fit some models with different options that influence output
## of getGenoPred function.
testFitMod1 <- fitModels(testTP, trait = "t1", quiet = TRUE)
testFitMod2 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                         quiet = TRUE)
testFitMod3 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)
testFitMod4 <- fitModels(testTP, trait = "t1", useCheck = TRUE, quiet = TRUE)
testFitMod5 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                         useCheck = TRUE, quiet = TRUE)
testFitMod6 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         useCheck = TRUE, quiet = TRUE)

### Check input.

expect_error(getGenoPred("fitMod"),
             "fitMod should be an object of class fitMod")

gp1 <- getGenoPred(testFitMod1)
gp2 <- getGenoPred(testFitMod2)
gp3 <- getGenoPred(testFitMod3)
gp4a <- getGenoPred(testFitMod4)
gp4b <- getGenoPred(testFitMod4, predictChecks = TRUE)
gp5 <- getGenoPred(testFitMod5)
gp6a <- getGenoPred(testFitMod6)
gp6b <- getGenoPred(testFitMod6, predictChecks = TRUE)

## Check output structure

expect_true(inherits(gp1, "list"))
expect_equal(names(gp1), c("genoPred", "checkPred"))
expect_null(gp1$checkPred)
expect_true(inherits(gp1$genoPred, "data.frame"))
expect_equal(dim(gp1$genoPred), c(110, 5))
expect_equal(colnames(gp1$genoPred), c("timeNumber", "timePoint", "genotype",
                                       "predicted.values", "standard.errors"))
# Extra column for geno.decomp.
# Check1 is in both replicates, so one extra prediction per time point.
expect_equal(dim(gp3$genoPred), c(115, 6))
expect_equal(setdiff(colnames(gp3$genoPred),
                     colnames(gp1$genoPred)), "geno.decomp")

# checkPred not empty for useCheck.
expect_equal(gp4a$genoPred, gp4b$genoPred)
expect_true(inherits(gp4b$checkPred, "data.frame"))
expect_equal(dim(gp4b$checkPred), c(5, 5))
expect_equal(colnames(gp4b$checkPred), c("timeNumber", "timePoint", "check",
                                         "predicted.values", "standard.errors"))

# Extra column for geno.decomp in check.
# Check1 is in both replicates, so one extra prediction per time point.
expect_equal(dim(gp6b$checkPred), c(10, 6))
expect_equal(setdiff(colnames(gp6b$genoPred),
                     colnames(gp4a$genoPred)), "geno.decomp")

## Check that results are as expected.
## Function is complicated and not always obvious so a thorough check is needed.

# Extract genotype and check predictions.
gp1 <- gp1$genoPred
gp2 <- gp2$genoPred
gp3 <- gp3$genoPred
gp4 <- gp4a$genoPred
gp4Check <- gp4b$checkPred
gp5 <- gp5$genoPred
gp6 <- gp6a$genoPred
gp6Check <- gp6b$checkPred

# Read expected results.
gp1Orig <- read.csv("gp1")
gp2Orig <- read.csv("gp2")
gp3Orig <- read.csv("gp3")
gp4Orig <- read.csv("gp4")
gp4CheckOrig <- read.csv("gp4c")
gp5Orig <- read.csv("gp5")
gp6Orig <- read.csv("gp6")
gp6CheckOrig <- read.csv("gp6c")

# timePoint, genotype and geno.decomp are imported from .csv as factors.
# Therefore checking if character values match with new results.
expect_equal(gp1[["timeNumber"]], gp1Orig[["timeNumber"]])
expect_equal(as.character(gp1[["timePoint"]]),
             as.character(gp1Orig[["timePoint"]]))
expect_equal(as.character(gp1[["genotype"]]),
             as.character(gp1Orig[["genotype"]]))
expect_equal(gp1[["predicted.values"]], gp1Orig[["predicted.values"]])
expect_equal(gp1[["standard.errors"]], gp1Orig[["standard.errors"]])

expect_equal(gp2[["timeNumber"]], gp2Orig[["timeNumber"]])
expect_equal(as.character(gp2[["timePoint"]]),
             as.character(gp2Orig[["timePoint"]]))
expect_equal(as.character(gp2[["genotype"]]),
             as.character(gp2Orig[["genotype"]]))
expect_equal(gp2[["predicted.values"]], gp2Orig[["predicted.values"]])
expect_equal(gp2[["standard.errors"]], gp2Orig[["standard.errors"]])

expect_equal(gp3[["timeNumber"]], gp3Orig[["timeNumber"]])
expect_equal(as.character(gp3[["timePoint"]]),
             as.character(gp3Orig[["timePoint"]]))
expect_equal(as.character(gp3[["geno.decomp"]]),
             as.character(gp3Orig[["geno.decomp"]]))
expect_equal(as.character(gp3[["genotype"]]),
             as.character(gp3Orig[["genotype"]]))
expect_equal(gp3[["predicted.values"]], gp3Orig[["predicted.values"]])
expect_equal(gp3[["standard.errors"]], gp3Orig[["standard.errors"]])

expect_equal(gp4[["timeNumber"]], gp4Orig[["timeNumber"]])
expect_equal(as.character(gp4[["timePoint"]]),
             as.character(gp4Orig[["timePoint"]]))
expect_equal(as.character(gp4[["genotype"]]),
             as.character(gp4Orig[["genotype"]]))
expect_equal(gp4[["predicted.values"]], gp4Orig[["predicted.values"]])
expect_equal(gp4[["standard.errors"]], gp4Orig[["standard.errors"]])

expect_equal(gp4Check[["timeNumber"]], gp4CheckOrig[["timeNumber"]])
expect_equal(as.character(gp4Check[["timePoint"]]),
             as.character(gp4CheckOrig[["timePoint"]]))
expect_equal(as.character(gp4Check[["check"]]),
             as.character(gp4CheckOrig[["check"]]))
expect_equal(gp4Check[["predicted.values"]], gp4CheckOrig[["predicted.values"]])
expect_equal(gp4Check[["standard.errors"]], gp4CheckOrig[["standard.errors"]])

expect_equal(gp5[["timeNumber"]], gp5Orig[["timeNumber"]])
expect_equal(as.character(gp5[["timePoint"]]),
             as.character(gp5Orig[["timePoint"]]))
expect_equal(as.character(gp5[["genotype"]]),
             as.character(gp5Orig[["genotype"]]))
expect_equal(gp5[["predicted.values"]], gp5Orig[["predicted.values"]])
expect_equal(gp5[["standard.errors"]], gp5Orig[["standard.errors"]])

expect_equal(gp6[["timeNumber"]], gp6Orig[["timeNumber"]])
expect_equal(as.character(gp6[["timePoint"]]),
             as.character(gp6Orig[["timePoint"]]))
expect_equal(as.character(gp6[["geno.decomp"]]),
             as.character(gp6Orig[["geno.decomp"]]))
expect_equal(as.character(gp6[["genotype"]]),
             as.character(gp6Orig[["genotype"]]))
expect_equal(gp6[["predicted.values"]], gp6Orig[["predicted.values"]])
expect_equal(gp6[["standard.errors"]], gp6Orig[["standard.errors"]])

expect_equal(gp6Check[["timeNumber"]], gp6CheckOrig[["timeNumber"]])
expect_equal(as.character(gp6Check[["timePoint"]]),
             as.character(gp6CheckOrig[["timePoint"]]))
expect_equal(as.character(gp6Check[["check"]]),
             as.character(gp6CheckOrig[["check"]]))
expect_equal(gp6Check[["predicted.values"]], gp6CheckOrig[["predicted.values"]])
expect_equal(gp6Check[["standard.errors"]], gp6CheckOrig[["standard.errors"]])

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getGenoPred(testFitMod4, outFile = "outfile"),
             "a single character string ending in .csv")
gpOut <- getGenoPred(testFitMod4, outFile = tmpFile)$genoPred
gpIn <- read.csv(tmpFile)
expect_equal(gpOut[["predicted.values"]], gpIn[["predicted.values"]])

# Checks should be written to a separate file.
tmpFile2 <- tempfile(fileext = ".csv")

gpOut2 <- getGenoPred(testFitMod4, predictChecks = TRUE, outFile = tmpFile2)
gpIn2 <- read.csv(tmpFile)
expect_equal(gpIn2, gpIn)
cpIn <- read.csv(gsub(pattern = ".csv", replacement = "Check.csv",
                      x = tmpFile2))
expect_equal(gpOut2$checkPred[["predicted.values"]], cpIn[["predicted.values"]])
