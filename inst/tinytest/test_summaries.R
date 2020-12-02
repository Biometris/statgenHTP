### Test summary functions.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP objects.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x", addCheck = TRUE,
                           checkGenotypes = "check1")
testTP2 <- createTimePoints(dat = testDat, experimentName = "testExp",
                            genotype = "Genotype", timePoint = "timepoints",
                            plotId = "pos", repId = "Replicate", rowNum = "y",
                            colNum = "x")

# Summary for TP with check genotypes differs slightly from
# TP without check genotypes.
sumTP1 <- capture.output(summary(testTP))
sumTP2 <- capture.output(summary(testTP2))

expect_true(any(grepl(pattern = "data for experiment testExp", x = sumTP1)))
expect_true(any(grepl(pattern = "contains 5 time points", x = sumTP1)))
expect_true(any(grepl(pattern = "First time point", x = sumTP1)))
expect_true(any(grepl(pattern = "Last time point", x = sumTP1)))
expect_true(any(grepl(pattern = "check genotypes: check1", x = sumTP1)))

expect_true(any(grepl(pattern = "No check genotypes", x = sumTP2)))

## Fit model for SpATS.
testFitMod <- fitModels(testTP, trait = "t1", quiet = TRUE)

## Summary for fitmod.
sumFitMod1 <- capture.output(summary(testFitMod))

expect_true(any(grepl(pattern = "fitted for experiment testExp",
                      x = sumFitMod1)))
expect_true(any(grepl(pattern = "contains 5 time points", x = sumFitMod1)))
expect_true(any(grepl(pattern = "fitted using SpATS", x = sumFitMod1)))

if (at_home() && FALSE) {
  ## Similar for asreml spatial.
  testFitMod2 <- fitModels(testTP, trait = "t1", engine = "asreml",
                           spatial = TRUE, quiet = TRUE)
  sumFitMod2 <- capture.output(summary(testFitMod2))

  expect_true(any(grepl(pattern = "fitted using asreml", x = sumFitMod2)))
  expect_true(any(grepl(pattern = "selected spatial model is AR1(x)id",
                        x = sumFitMod2, fixed = TRUE)))
  expect_true(any(grepl(pattern = "5 time points were used to select",
                        x = sumFitMod2)))
}
