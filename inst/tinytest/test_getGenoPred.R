### Test getGenoPred.

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
expect_warning(getGenoPred(testFitMod1, predictChecks = TRUE),
               "check was not used when fitting the model")

gp1 <- getGenoPred(testFitMod1)
gp2 <- getGenoPred(testFitMod2)
gp3 <- getGenoPred(testFitMod3)
gp4a <- getGenoPred(testFitMod4)
gp4b <- getGenoPred(testFitMod4, predictChecks = TRUE)
gp5 <- getGenoPred(testFitMod5)
gp6a <- getGenoPred(testFitMod6)
gp6b <- getGenoPred(testFitMod6, predictChecks = TRUE)

## Check output structure

expect_inherits(gp1, "list")
expect_equal(names(gp1), c("genoPred", "checkPred"))
expect_null(gp1$checkPred)
expect_inherits(gp1$genoPred, "data.frame")
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
expect_inherits(gp4b$checkPred, "data.frame")
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

expect_equal_to_reference(gp1, file = "gp1Comp")
expect_equal_to_reference(gp2, file = "gp2Comp")
expect_equal_to_reference(gp3, file = "gp3Comp")
expect_equal_to_reference(gp4, file = "gp4Comp")
expect_equal_to_reference(gp4Check, file = "gp4CheckComp")
expect_equal_to_reference(gp5, file = "gp5Comp")
expect_equal_to_reference(gp6, file = "gp6Comp")
expect_equal_to_reference(gp6Check, file = "gp6CheckComp")

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


## Test for models fitted using asreml.
## Limited set of tests, focused on differences between SpATS and asreml.

if (at_home()) {
  ## Fit some models with different options that influence output
  ## of getGenoPred function.
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

  gpAs1 <- getGenoPred(testFitModAs1)
  gpAs2 <- getGenoPred(testFitModAs2)
  gpAs3 <- getGenoPred(testFitModAs3)
  gpAs4a <- getGenoPred(testFitModAs4)
  gpAs4b <- getGenoPred(testFitModAs4, predictChecks = TRUE)
  gpAs5 <- getGenoPred(testFitModAs5)
  gpAs6a <- getGenoPred(testFitModAs6)
  gpAs6b <- getGenoPred(testFitModAs6, predictChecks = TRUE)

  ## Check that results are as expected.
  ## Function is complicated and not always obvious so a thorough check is needed.
  gpAs1 <- gpAs1$genoPred
  gpAs2 <- gpAs2$genoPred
  gpAs3 <- gpAs3$genoPred
  gpAs4 <- gpAs4a$genoPred
  gpAs4Check <- gpAs4b$checkPred
  gpAs5 <- gpAs5$genoPred
  gpAs6 <- gpAs6a$genoPred
  gpAs6Check <- gpAs6b$checkPred

  # Output structure should be identical to SpATS output.
  expect_equal(dim(gpAs1), dim(gp1))
  expect_equal(dim(gpAs2), dim(gp2))
  expect_equal(dim(gpAs3), dim(gp3))
  expect_equal(dim(gpAs4), dim(gp4))
  expect_equal(dim(gpAs4Check), dim(gp4Check))
  expect_equal(dim(gpAs5), dim(gp5))
  expect_equal(dim(gpAs6), dim(gp6))
  expect_equal(dim(gpAs6Check), dim(gp6Check))

  expect_equal(colnames(gpAs1), colnames(gp1))
  expect_equal(colnames(gpAs6), colnames(gp6))
  expect_equal(colnames(gpAs6Check), colnames(gp6Check))

  expect_equal_to_reference(gpAs1, file = "gpAs1Comp")
  expect_equal_to_reference(gpAs2, file = "gpAs2Comp")
  expect_equal_to_reference(gpAs3, file = "gpAs3Comp")
  expect_equal_to_reference(gpAs4, file = "gpAs4Comp")
  expect_equal_to_reference(gpAs4Check, file = "gpAs4CheckComp")
  expect_equal_to_reference(gpAs5, file = "gpAs5Comp")
  expect_equal_to_reference(gpAs6, file = "gpAs6Comp")
  expect_equal_to_reference(gpAs6Check, file = "gpAs6CheckComp")
}
