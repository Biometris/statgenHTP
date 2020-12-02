### Test getCorrected.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x", addCheck = TRUE,
                           checkGenotypes = "check1")

## Fit some models with different options that influence output
## of getCorrected function.
testFitMod1 <- fitModels(testTP, trait = "t1", quiet = TRUE)
testFitMod2 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                         quiet = TRUE)
testFitMod3 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)
# testFitMod4 <- fitModels(testTP, trait = "t1", useCheck = TRUE, quiet = TRUE)
# testFitMod5 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
#                          useCheck = TRUE, quiet = TRUE)
# testFitMod6 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
#                          useCheck = TRUE, quiet = TRUE)

### Check input.

expect_error(getCorrected("fitMod"),
             "fitMod should be an object of class fitMod")

corr1 <- getCorrected(testFitMod1)
corr2 <- getCorrected(testFitMod2)
corr3 <- getCorrected(testFitMod3)
# corr4 <- getCorrected(testFitMod4)
# corr5 <- getCorrected(testFitMod5)
# corr6 <- getCorrected(testFitMod6)

## Check output structure

expect_true(inherits(corr1, "data.frame"))
expect_equal(dim(corr1), c(125, 9))
expect_equal(colnames(corr1), c("timeNumber", "timePoint", "t1_corr", "t1",
                                "wt", "genotype", "rowId", "colId", "plotId"))

# Extra column for fixed factors and geno.decomp.
expect_equal(setdiff(colnames(corr2), colnames(corr1)), "Basin")
expect_equal(setdiff(colnames(corr3), colnames(corr1)), "geno.decomp")

## Check that results are as expected.
## Function is complicated and not always obvious so a thorough check is needed.

# Read expected results.
corr1Orig <- read.csv("corr1")
corr2Orig <- read.csv("corr2")
corr3Orig <- read.csv("corr3")
# corr4Orig <- read.csv("corr4")
# corr5Orig <- read.csv("corr5")
# corr6Orig <- read.csv("corr6")

if (FALSE) {


# timePoint, genotype, geno.decomp and plotId are imported from .csv as factors.
# Therefore checking if character values match with new results.
expect_equal(corr1[["timeNumber"]], corr1Orig[["timeNumber"]])
expect_equal(as.character(corr1[["timePoint"]]),
             as.character(corr1Orig[["timePoint"]]))
expect_equal(as.character(corr1[["genotype"]]),
             as.character(corr1Orig[["genotype"]]))
expect_equal(corr1[["t1"]], corr1Orig[["t1"]])
expect_equal(corr1[["t1_corr"]], corr1Orig[["t1_corr"]])
expect_equal(as.character(corr1[["plotId"]]),
             as.character(corr1Orig[["plotId"]]))

expect_equal(corr2[["timeNumber"]], corr2Orig[["timeNumber"]])
expect_equal(as.character(corr2[["timePoint"]]),
             as.character(corr2Orig[["timePoint"]]))
expect_equal(as.character(corr2[["genotype"]]),
             as.character(corr2Orig[["genotype"]]))
expect_equal(corr2[["t1"]], corr2Orig[["t1"]])
expect_equal(corr2[["t1_corr"]], corr2Orig[["t1_corr"]])
expect_equal(as.character(corr2[["plotId"]]),
             as.character(corr2Orig[["plotId"]]))
expect_equal(as.character(corr2[["Basin"]]),
             as.character(corr2Orig[["Basin"]]))

expect_equal(corr3[["timeNumber"]], corr3Orig[["timeNumber"]])
expect_equal(as.character(corr3[["timePoint"]]),
             as.character(corr3Orig[["timePoint"]]))
expect_equal(as.character(corr3[["geno.decomp"]]),
             as.character(corr3Orig[["geno.decomp"]]))
expect_equal(as.character(corr3[["genotype"]]),
             as.character(corr3Orig[["genotype"]]))
expect_equal(corr3[["t1"]], corr3Orig[["t1"]])
expect_equal(corr3[["t1_corr"]], corr3Orig[["t1_corr"]])
expect_equal(as.character(corr3[["plotId"]]),
             as.character(corr3Orig[["plotId"]]))

expect_equal(corr4[["timeNumber"]], corr4Orig[["timeNumber"]])
expect_equal(as.character(corr4[["timePoint"]]),
             as.character(corr4Orig[["timePoint"]]))
expect_equal(as.character(corr4[["genotype"]]),
             as.character(corr4Orig[["genotype"]]))
expect_equal(corr4[["t1"]], corr4Orig[["t1"]])
expect_equal(corr4[["t1_corr"]], corr4Orig[["t1_corr"]])
expect_equal(as.character(corr4[["plotId"]]),
             as.character(corr4Orig[["plotId"]]))

expect_equal(corr5[["timeNumber"]], corr5Orig[["timeNumber"]])
expect_equal(as.character(corr5[["timePoint"]]),
             as.character(corr5Orig[["timePoint"]]))
expect_equal(as.character(corr5[["genotype"]]),
             as.character(corr5Orig[["genotype"]]))
expect_equal(corr5[["t1"]], corr5Orig[["t1"]])
expect_equal(corr5[["t1_corr"]], corr5Orig[["t1_corr"]])
expect_equal(as.character(corr5[["plotId"]]),
             as.character(corr5Orig[["plotId"]]))
expect_equal(as.character(corr5[["Basin"]]),
             as.character(corr5Orig[["Basin"]]))

expect_equal(corr6[["timeNumber"]], corr6Orig[["timeNumber"]])
expect_equal(as.character(corr6[["timePoint"]]),
             as.character(corr6Orig[["timePoint"]]))
expect_equal(as.character(corr6[["geno.decomp"]]),
             as.character(corr6Orig[["geno.decomp"]]))
expect_equal(as.character(corr6[["genotype"]]),
             as.character(corr6Orig[["genotype"]]))
expect_equal(corr6[["t1"]], corr6Orig[["t1"]])
expect_equal(corr6[["t1_corr"]], corr6Orig[["t1_corr"]], tol = 1e-6)
expect_equal(as.character(corr6[["plotId"]]),
             as.character(corr6Orig[["plotId"]]))

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getCorrected(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
corrOut <- getCorrected(testFitMod1, outFile = tmpFile)
corrIn <- read.csv(tmpFile)
expect_equal(corrOut[["t1_corr"]], corrIn[["t1_corr"]])

## Test for models fitted using asreml.
## Limited set of tests, focused on diferences between SpATS and asreml.

## Fit some models with different options that influence output
## of getCorrected function.
testFitModAs1 <- fitModels(testTP, trait = "t1", engine = "asreml",
                           quiet = TRUE)
testFitModAs2 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                           engine = "asreml", quiet = TRUE)
testFitModAs3 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                           engine = "asreml", quiet = TRUE)
testFitModAs4 <- fitModels(testTP, trait = "t1", useCheck = TRUE,
                           engine = "asreml", quiet = TRUE)
testFitModAs5 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                           useCheck = TRUE, engine = "asreml", quiet = TRUE)
testFitModAs6 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                           useCheck = TRUE, engine = "asreml", quiet = TRUE)

expect_warning(corrAs1 <- getCorrected(testFitModAs1),
               "No spatial or fixed effects to correct for")
corrAs2 <- getCorrected(testFitModAs2)
expect_warning(corrAs3 <- getCorrected(testFitModAs3))
corrAs4 <- getCorrected(testFitModAs4)
corrAs5 <- getCorrected(testFitModAs5)
corrAs6 <- getCorrected(testFitModAs6)

# Output structure should be similar to SpATS output.
# No row + column correction, so rowId and colId not in output.
expect_equal(dim(corrAs1), dim(corr1) - c(0, 2))
expect_equal(dim(corrAs2), dim(corr2) - c(0, 2))
expect_equal(dim(corrAs3), dim(corr3) - c(0, 2))
expect_equal(dim(corrAs4), dim(corr4) - c(0, 2))
expect_equal(dim(corrAs5), dim(corr5) - c(0, 2))
expect_equal(dim(corrAs6), dim(corr6) - c(0, 2))

expect_equal(setdiff(colnames(corr1), colnames(corrAs1)), c("colId", "rowId"))
expect_equal(setdiff(colnames(corr6), colnames(corrAs6)), c("colId", "rowId"))

## Check that results are as expected.
## Function is complicated and not always obvious so a thorough check is needed.

# Read expected results.
corrAs1Orig <- read.csv("corrAs1")
corrAs2Orig <- read.csv("corrAs2")
corrAs3Orig <- read.csv("corrAs3")
corrAs4Orig <- read.csv("corrAs4")
corrAs5Orig <- read.csv("corrAs5")
corrAs6Orig <- read.csv("corrAs6")

# Ignoring timePoint and timeNumber since there is no difference with SpATS.
expect_equal(as.character(corrAs1[["genotype"]]),
             as.character(corrAs1Orig[["genotype"]]))
expect_equal(corrAs1[["t1"]], corrAs1Orig[["t1"]])
expect_equal(corrAs1[["t1_corrAs"]], corrAs1Orig[["t1_corrAs"]])
expect_equal(as.character(corrAs1[["plotId"]]),
             as.character(corrAs1Orig[["plotId"]]))

expect_equal(as.character(corrAs2[["genotype"]]),
             as.character(corrAs2Orig[["genotype"]]))
expect_equal(corrAs2[["t1"]], corrAs2Orig[["t1"]])
expect_equal(corrAs2[["t1_corrAs"]], corrAs2Orig[["t1_corrAs"]])
expect_equal(as.character(corrAs2[["plotId"]]),
             as.character(corrAs2Orig[["plotId"]]))
expect_equal(as.character(corrAs2[["Basin"]]),
             as.character(corrAs2Orig[["Basin"]]))

expect_equal(as.character(corrAs3[["geno.decomp"]]),
             as.character(corrAs3Orig[["geno.decomp"]]))
expect_equal(as.character(corrAs3[["genotype"]]),
             as.character(corrAs3Orig[["genotype"]]))
expect_equal(corrAs3[["t1"]], corrAs3Orig[["t1"]])
expect_equal(corrAs3[["t1_corrAs"]], corrAs3Orig[["t1_corrAs"]])
expect_equal(as.character(corrAs3[["plotId"]]),
             as.character(corrAs3Orig[["plotId"]]))

expect_equal(as.character(corrAs4[["genotype"]]),
             as.character(corrAs4Orig[["genotype"]]))
expect_equal(corrAs4[["t1"]], corrAs4Orig[["t1"]])
expect_equal(corrAs4[["t1_corrAs"]], corrAs4Orig[["t1_corrAs"]])
expect_equal(as.character(corrAs4[["plotId"]]),
             as.character(corrAs4Orig[["plotId"]]))

expect_equal(as.character(corrAs5[["genotype"]]),
             as.character(corrAs5Orig[["genotype"]]))
expect_equal(corrAs5[["t1"]], corrAs5Orig[["t1"]])
expect_equal(corrAs5[["t1_corrAs"]], corrAs5Orig[["t1_corrAs"]])
expect_equal(as.character(corrAs5[["plotId"]]),
             as.character(corrAs5Orig[["plotId"]]))
expect_equal(as.character(corrAs5[["Basin"]]),
             as.character(corrAs5Orig[["Basin"]]))

expect_equal(as.character(corrAs6[["geno.decomp"]]),
             as.character(corrAs6Orig[["geno.decomp"]]))
expect_equal(as.character(corrAs6[["genotype"]]),
             as.character(corrAs6Orig[["genotype"]]))
expect_equal(corrAs6[["t1"]], corrAs6Orig[["t1"]])
expect_equal(corrAs6[["t1_corrAs"]], corrAs6Orig[["t1_corrAs"]])
expect_equal(as.character(corrAs6[["plotId"]]),
             as.character(corrAs6Orig[["plotId"]]))



}
