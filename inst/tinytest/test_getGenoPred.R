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
gp4 <- getGenoPred(testFitMod4)
gp5 <- getGenoPred(testFitMod5)
gp6 <- getGenoPred(testFitMod6)

## Check output structure

expect_true(inherits(gp1, "data.frame"))
expect_equal(dim(gp1), c(120, 5))
expect_equal(colnames(gp1), c("timeNumber", "timePoint", "genotype",
                              "predicted.values", "standard.errors"))
# Extra column for geno.decomp.
expect_equal(dim(gp3), c(120, 6))
expect_equal(setdiff(colnames(gp3), colnames(gp1)), "geno.decomp")

## Check that results are as expected.
## Function is complicated and not always obvious so a thorough check is needed.

# Read expected results.
gp1Orig <- read.csv("gp1")
gp2Orig <- read.csv("gp2")
gp3Orig <- read.csv("gp3")
gp4Orig <- read.csv("gp4")
gp5Orig <- read.csv("gp5")
gp6Orig <- read.csv("gp6")

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
