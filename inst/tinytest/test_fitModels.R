### Test fitModels.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x", addCheck = TRUE,
                           checkGenotypes = "check1")

## Check that input testing works properly.
## With correct input an object with class fitMod should be returned.

expect_error(fitModels(), "TP should be an object of class TP")
expect_error(fitModels(testTP, trait = 1), "trait should be a character string")
expect_error(fitModels(testTP, trait = "trait"),
             "trait should be a column in TP for all timePoints")
expect_error(fitModels(testTP, trait = "t1", extraFixedFactors = 1),
             "extraFixedFactors should be a character vector")
expect_error(fitModels(testTP, trait = "t1", extraFixedFactors = "c1"),
             "c1 should be a column in TP for all timePoints")
expect_error(fitModels(testTP, trait = "t1", geno.decomp = 1),
             "geno.decomp should be a character vector")
expect_error(fitModels(testTP, trait = "t1", geno.decomp = "gd"),
             "gd should be a column for all timePoints")
expect_error(fitModels(testTP, trait = "t1", geno.decomp = "Basin",
                       what = "fixed"),
             "genotype as fixed effect and geno.decomp is not possible")
expect_error(fitModels(testTP, trait = "t1", useCheck = TRUE,
                       what = "fixed"),
             "genotype as fixed effect and useCheck = TRUE is not possible")

## Create TP object without row, col info, repId and check.
testTP2 <- createTimePoints(dat = testDat, experimentName = "testExp",
                            genotype = "Genotype", timePoint = "timepoints",
                            plotId = "pos")

expect_error(fitModels(testTP2, trait = "t1", useCheck = TRUE),
             "check should be a column in TP for all timePoints")
expect_error(fitModels(testTP2, trait = "t1", useRepId = TRUE),
             "repId should be a column in TP for all timePoints")
expect_error(fitModels(testTP2, trait = "t1"),
             "columns in TP for all timePoints when fitting spatial")

## Create TP object with all NA for t1.
## Time points with only NA should be skipped.
## Error when there is nothing left.
testTP3 <- testTP
testTP3[[1]][["t1"]] <- NA
expect_warning(fitModels(testTP3, trait = "t1", quiet = TRUE),
               "t1 has only NA values for the following time points")
expect_error(suppressWarnings(fitModels(testTP3, trait = "t1", timePoints = 1)),
             "No time points left for fitting models")

### Fit correct models and check that output is as expected.

testFitMod <- fitModels(testTP, trait = "t1", quiet = TRUE)
expect_inherits(testFitMod, "fitMod")

## Check that fitMod structure is correct.

expect_equal(length(testFitMod), 5)
expect_equal(names(testFitMod), names(testTP))
expect_inherits(testFitMod[[1]], "SpATS")

## Check that attributes are added correctly to fitMod object.

expect_equal(attr(testFitMod, "experimentName"), "testExp")
expect_equal(attr(testFitMod, "timePoints"), attr(testTP, "timePoints"))
expect_equal(attr(testFitMod, "what"), "random")
expect_false(attr(testFitMod, "useRepId"))
expect_false(attr(testFitMod, "spatial"))

## Check that correct models are fitted.
fixPart <- function(fitMod) {fitMod[[1]]$terms$fixed}
randPart <- function(fitMod) {fitMod[[1]]$terms$random}
## genVar gives a vector of length 2 or 3:
## genotype, geno.decomp and as.random -> all converted to character.
geno <- function(fitMod) {unlist(fitMod[[1]]$model$geno, use.names = FALSE)}

expect_equal(fixPart(testFitMod), formula("~1"))
expect_equal(randPart(testFitMod), formula("~rowId + colId"))
expect_equal(geno(testFitMod), c("genotype", "TRUE"))

## Repeat model check for different combinations of input variables.

## Add extra fixed factors.
testFitMod1 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                         quiet = TRUE)

expect_equal(fixPart(testFitMod1), formula("~Basin"))
expect_equal(randPart(testFitMod1), formula("~rowId + colId"))
expect_equal(geno(testFitMod1), c("genotype", "TRUE"))

## Add geno.decomp.
testFitMod2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)

expect_equal(fixPart(testFitMod2), formula("~geno.decomp"))
expect_equal(randPart(testFitMod2), formula("~rowId + colId"))
expect_equal(geno(testFitMod2), c("genotype", "geno.decomp", "TRUE"))

## Add extra fixed factors and geno.decomp
testFitMod3 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         extraFixedFactors = "Basin", quiet = TRUE)

expect_equal(fixPart(testFitMod3),
             formula("~Basin + geno.decomp + Basin:geno.decomp"))
expect_equal(randPart(testFitMod3), formula("~rowId + colId"))
expect_equal(geno(testFitMod3), c("genotype", "geno.decomp", "TRUE"))

## Add repId.
testfitMod4 <- fitModels(testTP, trait = "t1", useRepId = TRUE, quiet = TRUE)

expect_equal(fixPart(testfitMod4), formula("~repId"))
expect_equal(randPart(testfitMod4), formula("~repId:rowId + repId:colId"))
expect_equal(geno(testfitMod4), c("genotype", "TRUE"))

## use check.
testfitMod5 <- fitModels(testTP, trait = "t1", useCheck = TRUE, quiet = TRUE)

expect_equal(fixPart(testfitMod5), formula("~check"))
expect_equal(randPart(testfitMod5), formula("~rowId + colId"))
expect_equal(geno(testfitMod5), c("genoCheck", "TRUE"))

## Use check + geno.decomp.
testfitMod6 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         useCheck = TRUE, quiet = TRUE)

expect_equal(fixPart(testfitMod6),
             formula("~check + geno.decomp + check:geno.decomp"))
expect_equal(randPart(testfitMod6), formula("~rowId + colId"))
expect_equal(geno(testfitMod5), c("genoCheck", "TRUE"))


## Repeat model checks for models fitted using asreml.

if (at_home()) {

  fixPartAs <- function(fitMod) {fitMod[[1]]$formulae$fixed}
  randPartAs <- function(fitMod) {fitMod[[1]]$formulae$random}

  testFitModAs <- fitModels(testTP, trait = "t1", engine = "asreml",
                            quiet = TRUE)

  expect_equal(fixPartAs(testFitModAs), formula("t1~1"))
  expect_equal(randPartAs(testFitModAs), formula("~genotype"))

  ## Repeat model check for different combinations of input variables.

  ## Add extra fixed factors.
  testFitModAs1 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                             engine = "asreml", quiet = TRUE)

  expect_equal(fixPartAs(testFitModAs1), formula("t1~Basin"))
  expect_equal(randPartAs(testFitModAs1), formula("~genotype"))

  ## Add geno.decomp.
  testFitModAs2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                             engine = "asreml", quiet = TRUE)

  expect_equal(fixPartAs(testFitModAs2), formula("t1~geno.decomp"))
  expect_equal(randPartAs(testFitModAs2), formula("~at(geno.decomp):genotype"))

  ## Add repId.
  testFitModAs4 <- fitModels(testTP, trait = "t1", useRepId = TRUE,
                             engine = "asreml", quiet = TRUE)

  expect_equal(fixPartAs(testFitModAs4), formula("t1~repId"))
  expect_equal(randPartAs(testFitModAs4), formula("~genotype"))

  ## use check.
  testFitModAs5 <- fitModels(testTP, trait = "t1", useCheck = TRUE,
                             engine = "asreml", quiet = TRUE)

  expect_equal(fixPartAs(testFitModAs5), formula("t1~check"))
  expect_equal(randPartAs(testFitModAs5), formula("~genoCheck"))

  # Use check + geno.decomp.
  testFitModAs6 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                             useCheck = TRUE, engine = "asreml", quiet = TRUE)

  expect_equal(fixPartAs(testFitModAs6),
               formula("t1~check + geno.decomp + check:geno.decomp"))
  expect_equal(randPartAs(testFitModAs6), formula("~at(geno.decomp):genoCheck"))

  ## Check that missing plots are added when fitting spatial models.

  # Remove observation from testTP.
  testTP4 <- testTP
  testTP4[[1]] <- testTP4[[1]][1:24, ]
  testFitModAs7 <- fitModels(testTP4, trait = "t1", engine = "asreml",
                             spatial = TRUE, quiet = TRUE)
  # Check that row is added again to the data.
  expect_equal(nrow(testFitModAs7[[1]]$call$data), 25)
}
