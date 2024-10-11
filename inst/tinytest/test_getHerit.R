### Test getHerit.

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

expect_inherits(herit1, "data.frame")
expect_equal(dim(herit1), c(5, 3))
expect_equal(colnames(herit1), c("timeNumber", "timePoint", "h2"))

# Geno.decomp -> 1 column per level of geno.decomp.
expect_equal(dim(herit2), c(5, 4))
expect_equal(colnames(herit2), c("timeNumber", "timePoint", "1", "2"))

# Compare results.
expect_equal(herit1[["h2"]], c(0, 0.22, 0.02, 0.88, 0.43))
expect_equal(herit2[["1"]], c(0, 0.69, 0.02, 0.42, 0))
expect_equal(herit2[["2"]], c(0.3, 0.78, 0.59, 0.77, 0.73))

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getHerit(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
heritOut <- getHerit(testFitMod1, outFile = tmpFile)
heritIn <- read.csv(tmpFile)

# Ignore timePoint
# It is imported as factor, but is a date in the original data.
expect_equal(heritOut[, -2], heritIn[, -2])

## Test for models fitted using asreml.
## Limited set of tests, focused on diferences between SpATS and asreml.

if (at_home()) {
  ## Fit models.
  testFitModAs1 <- fitModels(testTP, trait = "t1", engine = "asreml",
                             quiet = TRUE)
  testFitModAs2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                             engine = "asreml", quiet = TRUE)
  testFitModAs3 <- fitModels(testTP, trait = "t1", what = "fixed",
                             engine = "asreml", quiet = TRUE)

  expect_error(getHerit(testFitModAs3),
               "Heritability can only be calculated when genotype is random")

  heritAs1 <- getHerit(testFitModAs1)
  heritAs2 <- getHerit(testFitModAs2)

  expect_equal(heritAs1[["h2"]],
               c(-2.02929559081078e-05, 0.167450198761112,
                 -4.37458765323306e-06, 0.650245799271745, 0.999982638701502))
  expect_equal(heritAs2[["1"]],
               c(-1.8990975917399e-05, 0.937541193808098, 0.886355776919209,
                 0.286419906650437, 8.85836689556996e-08))
  expect_equal(heritAs2[["2"]],
               c(-1.9124297023998e-05, 0.951117118858154, 0.941710391410272,
                 0.767914363565008, 0.646403039223904))
}

