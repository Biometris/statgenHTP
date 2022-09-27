### Test predict.psHDM

Sys.setlocale("LC_COLLATE", "C")

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)

## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           repId = "Replicate", plotId = "pos",
                           rowNum = "y", colNum = "x")

## Fit model.
testFitMod <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                        extraFixedFactors = "Basin", quiet = TRUE)

## Get corrected values.
corr <- getCorrected(testFitMod)
## Move all check1s to Basin 1 for a consistent population structure.
corr[corr[["genotype"]] == "check1", "Basin"] <- 1

## Fit a HDM spline on the corrected values.
splineRes <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          trace = FALSE)

## General input checks.
expect_error(predict(splineRes, newtimes = "a"),
             "newtimes should be a numerical vector")


## Make predictions using default settings.
pred <- predict(splineRes, trace = FALSE)

expect_inherits(pred, "psHDM")
expect_equal(length(pred), 4)
expect_equal(names(pred),
             c("newtimes", "popLevel", "genoLevel", "plotLevel"))

## Check structure of components.
expect_inherits(pred$newtimes, "numeric")
expect_inherits(pred$popLevel, "data.frame")
expect_inherits(pred$genoLevel, "data.frame")
expect_inherits(pred$plotLevel, "data.frame")

## Check results
expect_equal_to_reference(pred, "predHDM")


## Check that option newtimes functions correctly.
pred2 <- predict(splineRes, newtimes = 0)

expect_equal(pred2$newtimes, 0)
expect_equivalent(pred2$popLevel,
                  pred$popLevel[pred$popLevel$timeNumber == 0, ])
expect_equivalent(pred2$genoLevel,
                  pred$genoLevel[pred$genoLevel$timeNumber == 0, ])
expect_equivalent(pred2$plotLevel,
                  pred$plotLevel[pred$plotLevel$timeNumber == 0, 1:11])

