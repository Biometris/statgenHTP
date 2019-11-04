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

testFitMod <- fitModels(testTP, trait = "t1", quiet = TRUE)
expect_true(inherits(testFitMod, "fitMod"))

## Check that fitMod structure is correct.

expect_equal(length(testFitMod), 5)
expect_equal(names(testFitMod), names(testTP))
expect_true(inherits(testFitMod[[1]], "SpATS"))

## Check that attributes are added correctly to fitMod object.

expect_equal(attr(testFitMod, "experimentName"), "testExp")
expect_equal(attr(testFitMod, "timePoints"), attr(testTP, "timePoints"))
expect_equal(attr(testFitMod, "what"), "random")
expect_false(attr(testFitMod, "useRepId"))
expect_false(attr(testFitMod, "spatial"))

## Check that correct models are fitted.
fix <- function(fitMod) {fitMod[[1]]$terms$fixed}
rand <- function(fitMod) {fitMod[[1]]$terms$random}
## genVar gives a vector of length 2 or 3:
## genotype, geno.decomp and as.random -> all converted to character.
geno <- function(fitMod) {unlist(fitMod[[1]]$model$geno, use.names = FALSE)}

expect_equal(fix(testFitMod), formula("~1"))
expect_equal(rand(testFitMod), formula("~rowId + colId"))
expect_equal(geno(testFitMod), c("genotype", "TRUE"))

## Repeat model check for different combinations of input variables.

## Add extra fixed factors.
testFitMod1 <- fitModels(testTP, trait = "t1", extraFixedFactors = "Basin",
                         quiet = TRUE)

expect_equal(fix(testFitMod1), formula("~Basin"))
expect_equal(rand(testFitMod1), formula("~rowId + colId"))
expect_equal(geno(testFitMod1), c("genotype", "TRUE"))

## Add geno.decomp.
testFitMod2 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         quiet = TRUE)

expect_equal(fix(testFitMod2), formula("~geno.decomp"))
expect_equal(rand(testFitMod2), formula("~rowId + colId"))
expect_equal(geno(testFitMod2), c("genotype", "geno.decomp", "TRUE"))

## Add extra fixed factors and geno.decomp
testFitMod3 <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                         extraFixedFactors = "Basin", quiet = TRUE)

expect_equal(fix(testFitMod3), formula("~Basin"))
expect_equal(rand(testFitMod3), formula("~rowId + colId"))
expect_equal(geno(testFitMod3), c("genotype", "geno.decomp", "TRUE"))

## Add repId.
testfitMod4 <- fitModels(testTP, trait = "t1", useRepId = TRUE, quiet = TRUE)

expect_equal(fix(testfitMod4), formula("~repId"))
expect_equal(rand(testfitMod4), formula("~repId:rowId + repId:colId"))
expect_equal(geno(testfitMod4), c("genotype", "TRUE"))
