### Test getVar

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

expect_inherits(var1, "data.frame")
expect_equal(dim(var1), c(5, 6))
expect_equal(colnames(var1), c("timeNumber", "timePoint", "varGen", "varRes",
                               "varCol", "varRow"))

# Extra column for 2 levels in geno.decomp.
expect_equal(dim(var2), c(5, 7))
expect_equal(colnames(var2), c("timeNumber", "timePoint", "var_geno.decomp_1",
                               "var_geno.decomp_2", "varRes", "varCol",
                               "varRow"))

# Compare results.

expect_equal_to_reference(var1, file = "var1Comp")
expect_equal_to_reference(var2, file = "var2Comp")
expect_equal_to_reference(var3, file = "var3Comp")
expect_equal_to_reference(var4, file = "var4Comp")

## Check that results can be written to a file.
tmpFile <- tempfile(fileext = ".csv")

expect_error(getVar(testFitMod1, outFile = "outfile"),
             "a single character string ending in .csv")
varOut <- getVar(testFitMod1, outFile = tmpFile)
varIn <- read.csv(tmpFile, stringsAsFactors = FALSE)

# Ignore timePoint
# It is imported as character, but is a date in the original data.
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
  ## This sometimes gives warnings about oscillating parameters.
  testFitModAs4 <- fitModels(testTP, trait = "t1", useRepId = TRUE,
                             engine = "asreml", spatial = TRUE, quiet = TRUE)

  varAs1 <- getVar(testFitModAs1)
  varAs2 <- getVar(testFitModAs2)
  varAs3 <- getVar(testFitModAs3)
  varAs4 <- getVar(testFitModAs4)

  expect_equal(varAs1[["varGen"]],
               c(1.65268879494901e-09, 5.64519368545327e-05,
                 3.15261383070164e-10, 0.000612456995455758,
                 2.03348262780288e-05))
  expect_equal(varAs1[["varRes"]],
               c(0.00103293042234413, 0.000719166385240099, 0.00095329228503577,
                 0.000342990030581744, 0.000639896746661389))
  expect_equal(varAs1[["varCol"]], rep(NA_real_, times = 5))
  expect_equal(varAs1[["varRow"]], rep(NA_real_, times = 5))

  expect_equal(varAs2[["var_geno.decomp_1"]],
               c(1.72671949906827e-09, 0.000593631973665045,
                 0.000566468562072854, 0.000130418179182102,
                 5.90967913706669e-10))
  expect_equal(varAs2[["var_geno.decomp_2"]],
               c(1.72671949906827e-09, 0.000807184697098741,
                 0.00123225562394439, 0.00114757984666827,
                 0.000675211544349767))
  expect_equal(varAs2[["varRes"]],
               c(0.00107919960908156, 4.15070259538522e-05,
                 7.62725046737966e-05, 0.000346567980990889,
                 0.000369354919427342))

  expect_equal(varAs3[["varGen"]],
               c(1.65397675029422e-09, 0.000312114772570697,
                 3.1433085674295e-05, 0.000626495864056879,
                 0.000479038160766828))
  expect_equal(varAs3[["varRes"]],
               c(0.00103373539437683, 0.000514334135629192,
                 0.000902330012987955, 0.000336341338876986,
                 0.000284339025219647))

  expect_equal(varAs4[["varGen"]],
               c(1.38646720134525e-09, 1.00887137199852e-09,
                 0.000147384723749649, 0.000729521481678914,
                 0.000373421169164799))
  expect_equal(varAs4[["varRes"]],
               c(0.00086654193834238, 0.00063054456202173,
                 0.000456828756440842, 0.000173447996381509,
                 0.000389710248725436))
}
