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

## Test for models fitted using asreml.
## Limited set of tests, focused on diferences between SpATS and asreml.

if (at_home()) {
  ## Fit models.
  testFitModAs1 <- fitModels(testTP, trait = "t1", engine = "asreml",
                             quiet = TRUE)
  testFitModAs2 <- fitModels(testTP, trait = "t1", engine = "asreml",
                             geno.decomp = "repId", quiet = TRUE)
  testFitModAs3 <- fitModels(testTP, trait = "t1", useCheck = TRUE,
                             engine = "asreml", quiet = TRUE)
  ## This gives a lot of warnings about oscillating parameters.
  expect_warning(testFitModAs4 <-
                   fitModels(testTP, trait = "t1", useRepId = TRUE,
                             engine = "asreml", spatial = TRUE, quiet = TRUE),
                 "Oscillating parameter")

  varAs1 <- getVar(testFitModAs1)
  varAs2 <- getVar(testFitModAs2)
  varAs3 <- getVar(testFitModAs3)
  varAs4 <- getVar(testFitModAs4)

  expect_equal(varAs1[["varGen"]],
               c(1.65268879494901e-09, 5.64519368545327e-05, 3.15261383070164e-10,
                 0.000612456995455758, 2.03348262780288e-05))
  expect_equal(varAs1[["varRes"]],
               c(0.00103293042234413, 0.000719166385240099, 0.00095329228503577,
                 0.000342990030581744, 0.000639896746661389))
  expect_equal(varAs1[["varCol"]], rep(NA_real_, times = 5))
  expect_equal(varAs1[["varRow"]], rep(NA_real_, times = 5))

  expect_equal(varAs3[["varGen"]],
               c(1.27155338625288e-09, 0.000319265101756853, 1.14364124067873e-09,
                 0.000600635211857471, 0.000489856464945024))
  expect_equal(varAs3[["varRes"]],
               c(0.000794720809089666, 0.000516941433348486, 0.000714775723871777,
                 0.000348842294942067, 0.000290678907005537))

  expect_equal(varAs4[["varGen"]],
               c(1.38646720134564e-09, 1.00887137199698e-09, 0.000147384723762021,
                 0.000730185616933471, 0.000373421169180005))
  expect_equal(varAs4[["varRes"]],
               c(0.000866541938342623, 0.000630544562020765, 0.000456828756445013,
                 0.000175763890272973, 0.000389710248716934))
}





