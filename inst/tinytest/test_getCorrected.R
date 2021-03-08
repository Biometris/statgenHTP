### Test getCorrected.

Sys.setlocale("LC_COLLATE", "C")

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
testFitMod4 <- fitModels(testTP, trait = "t1", useCheck = TRUE, quiet = TRUE)
testFitMod5 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                         useCheck = TRUE, quiet = TRUE)
testFitMod6 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         useCheck = TRUE, quiet = TRUE)

### Check input.

expect_error(getCorrected("fitMod"),
             "fitMod should be an object of class fitMod")

corr1 <- getCorrected(testFitMod1)
corr2 <- getCorrected(testFitMod2)
corr3 <- getCorrected(testFitMod3)
corr4 <- getCorrected(testFitMod4)
corr5 <- getCorrected(testFitMod5)
corr6 <- getCorrected(testFitMod6)

## Check output structure

expect_inherits(corr1, "data.frame")
expect_equal(dim(corr1), c(125, 9))
expect_equal(colnames(corr1), c("timeNumber", "timePoint", "t1_corr", "t1",
                                "wt", "genotype", "rowId", "colId", "plotId"))

# Extra column for fixed factors and geno.decomp.
expect_equal(setdiff(colnames(corr2), colnames(corr1)), "Basin")
expect_equal(setdiff(colnames(corr3), colnames(corr1)), "geno.decomp")

## Check that results are as expected.
## Function is complicated and not always obvious so a thorough check is needed.

expect_equal_to_reference(corr1, file = "corr1Comp", tolerance = 1e-5)
expect_equal_to_reference(corr2, file = "corr2Comp", tolerance = 1e-5)
expect_equal_to_reference(corr3, file = "corr3Comp", tolerance = 1e-5)
expect_equal_to_reference(corr4, file = "corr4Comp", tolerance = 1e-5)
expect_equal_to_reference(corr5, file = "corr5Comp", tolerance = 1e-5)
expect_equal_to_reference(corr6, file = "corr6Comp", tolerance = 1e-5)

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getCorrected(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
corrOut <- getCorrected(testFitMod1, outFile = tmpFile)
corrIn <- read.csv(tmpFile)
expect_equal(corrOut[["t1_corr"]], corrIn[["t1_corr"]])

## Test for models fitted using asreml.
## Limited set of tests, focused on differences between SpATS and asreml.

if (at_home()) {
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
  expect_warning(corrAs3 <- getCorrected(testFitModAs3),
                 "No spatial or fixed effects to correct for")
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

  expect_equal(setdiff(colnames(corr1), colnames(corrAs1)), c("rowId", "colId"))
  expect_equal(setdiff(colnames(corr6), colnames(corrAs6)), c("rowId", "colId"))

  ## Check that results are as expected.
  ## Function is complicated and not always obvious so a thorough check is needed.

  expect_equal_to_reference(corrAs1, file = "corrAs1Comp")
  expect_equal_to_reference(corrAs2, file = "corrAs2Comp")
  expect_equal_to_reference(corrAs3, file = "corrAs3Comp")
  expect_equal_to_reference(corrAs4, file = "corrAs4Comp")
  expect_equal_to_reference(corrAs5, file = "corrAs5Comp")
  expect_equal_to_reference(corrAs6, file = "corrAs6Comp")
}
