### Test createTP.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Add extra row for testing uniqueness criterion for plotId.
testDat2 <- rbind(testDat, testDat[1, ])

## Check that input testing works properly.
## With correct input an object with class TP should be returned.
expect_error(createTimePoints(), "dat has to be a data.frame")
expect_error(createTimePoints(dat = testDat), "experimentName should be")
expect_error(createTimePoints(dat = testDat, experimentName = "testExp"),
             "argument \"genotype\" is missing")
expect_error(createTimePoints(dat = testDat, experimentName = "testExp",
                              genotype = "Genotype"),
             "argument \"timePoint\" is missing")
expect_error(createTimePoints(dat = testDat, experimentName = "testExp",
                              genotype = "Genotype", timePoint = "timepoints"),
             "argument \"plotId\" is missing")
expect_error(createTimePoints(dat = testDat, experimentName = "testExp",
                              genotype = "genotype", timePoint = "timepoints",
                              plotId = "pos"),
             "genotype has to be NULL or a column in dat")
expect_error(createTimePoints(dat = testDat, experimentName = "testExp",
                              genotype = "Genotype", timePoint = "timepoints",
                              plotId = "pos", addCheck = TRUE),
             "checkGenotypes should be a character vector")
expect_error(createTimePoints(dat = testDat, experimentName = "testExp",
                              genotype = "Genotype", timePoint = "timepoints",
                              plotId = "pos", addCheck = TRUE,
                              checkGenotypes = "c1"),
             "All checkGenotypes should be genotypes in dat")
expect_error(createTimePoints(dat = testDat2, experimentName = "testExp",
                              genotype = "Genotype", timePoint = "timepoints",
                              plotId = "pos"),
             "plotId has to be unique within each time point.")
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos")
expect_true(inherits(testTP, "TP"))

## Check that option timeFormat is used correctly.
## Supplying incompatible timeFormat should give an error.
expect_error(createTimePoints(dat = testDat, experimentName = "testExp",
                              genotype = "Genotype", timePoint = "timepoints",
                              timeFormat = "%Y/%m/%d", plotId = "pos"),
             "Error when converting timePoints to Date format")
## Conversion of numeric columns is possible, but might give unexpected results.
testDatTF <- testDat
testDatTF[["timepoints"]] <- as.numeric(as.factor(testDatTF[["timepoints"]]))
expect_warning(createTimePoints(dat = testDatTF, experimentName = "testExp",
                                genotype = "Genotype", timePoint = "timepoints",
                                plotId = "pos"),
               "timePoint is a numeric column. It will be converted")
## Compatible timeFormat should be ok.
expect_silent(createTimePoints(dat = testDat, experimentName = "testExp",
                               genotype = "Genotype", timePoint = "timepoints",
                               timeFormat = "%Y-%m-%d %H:%M:%S", plotId = "pos"))

## Check that row and column numbers are added correctly.
## Should be added both as factor and as numeric column.
testTPrc <- createTimePoints(dat = testDat, experimentName = "testExp",
                             genotype = "Genotype", timePoint = "timepoints",
                             plotId = "pos", rowNum = "y", colNum = "x")
expect_true(hasName(testTPrc[[1]], "colId"))
expect_true(hasName(testTPrc[[1]], "rowId"))
expect_true(hasName(testTPrc[[1]], "colNum"))
expect_true(hasName(testTPrc[[1]], "rowNum"))
expect_true(inherits(testTPrc[[1]][["colId"]], "factor"))
expect_true(inherits(testTPrc[[1]][["rowId"]], "factor"))
expect_true(inherits(testTPrc[[1]][["colNum"]], "numeric"))
expect_true(inherits(testTPrc[[1]][["rowNum"]], "numeric"))

## Check that check genotypes are added correctly.
## When present two columns should be added: check and genoCheck.
testTPck <- createTimePoints(dat = testDat, experimentName = "testExp",
                             genotype = "Genotype", timePoint = "timepoints",
                             plotId = "pos", addCheck = TRUE,
                             checkGenotypes = "check1")
expect_true(hasName(testTPck[[1]], "check"))
expect_true(hasName(testTPck[[1]], "genoCheck"))
expect_true(inherits(testTPck[[1]][["check"]], "factor"))
expect_true(inherits(testTPck[[1]][["genoCheck"]], "factor"))
expect_equal(levels(testTPck[[1]][["check"]]), c("check1", "noCheck"))
expect_equal(sum(is.na(testTPck[[1]][["genoCheck"]])), 2)

## Check that attributes are added correctly.
expect_equal(names(testTP), unique(testDat[["timepoints"]]))
expect_equal(attr(testTP, "experimentName"), "testExp")

# Extract timePoints attribute.
timePoints <- attr(testTP, "timePoints")
expect_true(inherits(timePoints, "data.frame"))
expect_equal(names(timePoints), c("timeNumber", "timePoint"))
expect_equal(timePoints[["timeNumber"]], 1:5)
expect_equal(timePoints[["timePoint"]], names(testTP))

### Test subsetting TP.

## Check that input checking works properly.
expect_error(testTP[6], "All timePoints should be in x")
expect_error(testTP[c(1, 6)], "All timePoints should be in x")
expect_error(testTP["2018-05-01 16:37:00"], "All timePoints should be in x")

## Check that subsetting works with integer selection.
testTP23 <- testTP[2:3]
expect_true(inherits(testTP23, "TP"))
expect_equal(names(testTP23), names(testTP)[2:3])
expect_equal(attr(testTP23, "timePoints"), timePoints[2:3, ])

## Check that subsetting works with character selection.
timePoints15 <- unique(testDat[["timepoints"]])[c(1, 5)]
testTP15 <- testTP[timePoints15]
expect_true(inherits(testTP15, "TP"))
expect_equal(names(testTP15), names(testTP)[c(1, 5)])
expect_equal(attr(testTP15, "timePoints"), timePoints[c(1, 5), ])
