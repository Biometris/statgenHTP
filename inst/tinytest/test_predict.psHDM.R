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
expect_error(predict(splineRes, pred = "a"),
             "pred should be a named list of length 3")
expect_error(predict(splineRes,
                     pred = list(pop = FALSE, geno = TRUE, plot = TRUE)),
             paste("Predictions at plot level can only be made if predictions",
                   "are also made at geno and pop level"))
expect_error(predict(splineRes,
                     pred = list(pop = FALSE, geno = TRUE, plot = FALSE)),
             paste("Predictions at geno level can only be made if predictions",
                   "are also made at pop level"))

expect_error(predict(splineRes, se = "a"),
             "se should be a named list of length 3")
expect_error(predict(splineRes,
                     se = list(pop = FALSE, geno = TRUE, plot = TRUE)),
             paste("Standard errors at plot level can only be computed if",
                   "standard errors are also computed at geno and pop level"))
expect_error(predict(splineRes,
                     se = list(pop = FALSE, geno = TRUE, plot = FALSE)),
             paste("Standard errors at geno level can only be computed if",
                   "standard errors are also computed at pop level"))

expect_error(predict(splineRes,
                     pred = list(pop = FALSE, geno = FALSE, plot = FALSE),
                     se = list(pop = TRUE, geno = TRUE, plot = TRUE)),
             paste("Standard errors at population level can only be computed",
                   "if predictions are also made at population level"))
expect_error(predict(splineRes,
                     pred = list(pop = TRUE, geno = FALSE, plot = FALSE),
                     se = list(pop = TRUE, geno = TRUE, plot = TRUE)),
             paste("Standard errors at genotype level can only be computed",
                   "if predictions are also made at genotype level"))
expect_error(predict(splineRes,
                     pred = list(pop = TRUE, geno = TRUE, plot = FALSE),
                     se = list(pop = TRUE, geno = TRUE, plot = TRUE)),
             paste("Standard errors at plot level can only be computed if",
                   "predictions are also made at plot level"))


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
pred2 <- predict(splineRes, newtimes = 0, trace = FALSE)

expect_equal(pred2$newtimes, 0)
expect_equivalent(pred2$popLevel,
                  pred$popLevel[pred$popLevel$timeNumber == 0, ])
expect_equivalent(pred2$genoLevel,
                  pred$genoLevel[pred$genoLevel$timeNumber == 0, ])
expect_equivalent(pred2$plotLevel,
                  pred$plotLevel[pred$plotLevel$timeNumber == 0, 1:11])

## Check that computation of standard errors functions correctly.

pred3 <- predict(splineRes,
                 pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                 se = list(pop = TRUE, geno = TRUE, plot = TRUE),
                 trace = FALSE)

expect_equal_to_reference(pred3, "predHDMse")

## Check that option trace functions correctly.
expect_silent(predict(splineRes, trace = FALSE))

expect_stdout(predict(splineRes, trace = TRUE),
              "Population-specific growth curves OK")


