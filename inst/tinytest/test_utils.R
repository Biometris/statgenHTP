### Test calcPlotBorders.
tpDat <- data.frame(rowNum = rep(1:4, each = 4),
                    colNum = rep(1:4, times = 4),
                    repId = rep(1:2, each = 8))

## No missing values. Borders between replicates and around edge.
bord <- statgenHTP:::calcPlotBorders(tpDat = tpDat, bordVar = "repId")
expect_inherits(bord, "list")
expect_equal(names(bord), c("horW", "vertW"))
expect_inherits(bord$horW, "data.frame")
expect_inherits(bord$vertW, "data.frame")
expect_equal(bord$horW[["x"]], rep(1:4, each = 3))
expect_equal(bord$horW[["y"]], rep(c(1, 3, 5), times = 4))
expect_equal(bord$vertW[["x"]], rep(c(1, 5), times = 4))
expect_equal(bord$vertW[["y"]], rep(1:4, each = 2))

## Missing value. Only differs from above by borders moved inwards.
tpDat <- tpDat[-1, ]
bord <- statgenHTP:::calcPlotBorders(tpDat = tpDat, bordVar = "repId")
expect_equal(bord$horW[["x"]], c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4))
expect_equal(bord$horW[["y"]], c(3, 5, 1, 3, 5, 1, 3, 5, 1, 3, 5))
expect_equal(bord$vertW[["x"]], c(5, 1, 5, 1, 5, 1, 5))
expect_equal(bord$vertW[["y"]], c(1, 2, 2, 3, 3, 4, 4))

### Test dfBind

## Columns should be copied properly.
## Check that naming of the output is correct.
df1 <- data.frame(a = 1:2, b = 1:2)
df2 <- data.frame(a = 1:2, c = 1:2)
df3 <- data.frame(c = 1:2, d = 1:2)
expect_equal(colnames(statgenHTP:::dfBind(list(df1, df1))),
             c("a", "b"))
expect_equal(colnames(statgenHTP:::dfBind(list(df1, df2))),
             c("a", "b", "c"))
expect_equal(colnames(statgenHTP:::dfBind(list(df1, df3))),
             c("a", "b", "c", "d"))
expect_equal(colnames(statgenHTP:::dfBind(list(df1, df2, df3))),
             c("a", "b", "c", "d"))

## Check that content of the output is correct.
## NA should be inserted for missing columns.
expect_equivalent(unlist(statgenHTP:::dfBind(list(df1, df2))),
                  c(1, 2, 1, 2, 1, 2, NA, NA, NA, NA, 1, 2))
expect_equivalent(unlist(statgenHTP:::dfBind(list(df1, df2, df1))),
                  c(1, 2, 1, 2, 1,2, 1, 2, NA, NA, 1, 2, NA, NA, 1, 2, NA, NA))

## Check that empty data.frames are removed before binding.
expect_equal(statgenHTP:::dfBind(list(data.frame(), df1)), df1)
expect_equal(statgenHTP:::dfBind(list(df1, data.frame())), df1)
expect_equal(statgenHTP:::dfBind(list(data.frame())), data.frame())

### Test addMissVals

times <- strptime(c("1sep2019", "2sep2019"), "%d%b%Y")
timeNums <- c(1, 1, 2)
df1 <- data.frame(timePoint = times[timeNums], timeNumber = timeNums,
                  id = c("p1", "p2", "p1"), trait = 1:3,
                  stringsAsFactors = FALSE)

## Check that single missing value is added when no extra columns present.
dfOut1 <- statgenHTP:::addMissVals(df1, "trait")
expect_inherits(dfOut1, "data.frame")
expect_inherits(dfOut1[["timePoint"]], "POSIXct")
expect_equal(nrow(dfOut1), 4)
expect_equivalent(unlist(dfOut1[4, ]), c("p2", "1567382400", NA))

## Check that multiple missing values are added when no extra columns present.
df2 <- df1[-1, ]
dfOut2 <- statgenHTP:::addMissVals(df2, "trait")
expect_equal(nrow(dfOut2), 4)
expect_equivalent(unlist(dfOut2[1, ]), c("p1", "1567296000", NA))
expect_equivalent(unlist(dfOut2[4, ]), c("p2", "1567382400", NA))

## Check that single missing value is added when extra columns present.
df3 <- df1
df3[["treat"]] <- c("W", "D", "W")
dfOut3 <- statgenHTP:::addMissVals(df3, "trait")
expect_equal(dfOut3[colnames(dfOut3) != "treat"], dfOut1)
expect_equal(dfOut3[["treat"]], c("W", "D", "W" , "D"))

## Check that multiple missing values are added when extra columns present.
df4 <- df2
df4[["treat"]] <- c("D", "W")
dfOut4 <- statgenHTP:::addMissVals(df4, "trait")
expect_equal(dfOut4[colnames(dfOut4) != "treat"], dfOut2)
expect_equal(dfOut4[["treat"]], c("W", "D", "W" , "D"))

### Test prettier

## Different splits for n = 1, 2 and 3.
## Use 40 day time difference for ease of checking.
times <- strptime(c("1sep2019", "11oct2019"), "%d%b%Y")
expect_equal(statgenHTP:::prettier(n = 1)(times),
             strptime("21sep2019", "%d%b%Y"))
expect_equal(statgenHTP:::prettier(n = 2)(times),
             strptime(c("9sep2019", "3oct2019"), "%d%b%Y"))
expect_equal(statgenHTP:::prettier(n = 3)(times),
             strptime(c("6sep2019", "21sep2019", "6oct2019"), "%d%b%Y"))

### Test chkFile

tmpFile <- tempfile(fileext = ".csv")

expect_error(statgenHTP:::chkFile(1),
             "outFile should be a single character string ending in .csv")
expect_error(statgenHTP:::chkFile(tmpFile, fileType = "pdf"),
             "outFile should be a single character string ending in .pdf")
expect_silent(statgenHTP:::chkFile(tmpFile))

### Test chkTimePoints

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", repId = "Replicate", rowNum = "y",
                           colNum = "x")

## Check that input checks function correctly.
expect_error(statgenHTP:::chkTimePoints("testTP", 1),
             "should be an object of class TP or fitMod")
expect_error(statgenHTP:::chkTimePoints(testTP, strptime("1sep2019", "%d%b%Y")),
             "timePoints should be a character or numeric vector")
expect_error(statgenHTP:::chkTimePoints(testTP, 6),
             "All timePoints should be in testTP")
expect_error(statgenHTP:::chkTimePoints(testTP, c(1, 6)),
             "All timePoints should be in testTP")
expect_error(statgenHTP:::chkTimePoints(testTP, "2018-05-01 16:37:00"),
             "All timePoints should be in testTP")

## Check that integer subset is converted to character.
expect_inherits(statgenHTP:::chkTimePoints(testTP, 1), "character")
expect_equal(statgenHTP:::chkTimePoints(testTP, 1),
             statgenHTP:::chkTimePoints(testTP, "2018-06-01 16:37:00"))

### Test getTimePoints

expect_error(getTimePoints(x = 1:3),
             "x should be an object of class TP or fitMod")
expect_inherits(getTimePoints(testTP), "data.frame")
expect_equal(getTimePoints(testTP), attr(x = testTP, which = "timePoints"))

## Test countValid
expect_error(countValid(1),
             "TP should be an object of class TP")
expect_error(countValid(testTP, trait = 1),
             "trait should be a character string of length one")
valCount <- countValid(testTP, trait = "t1")
expect_equivalent(valCount, rep(24, times = 5))
expect_equal(names(valCount),
             c("2018-06-01 16:37:00", "2018-06-02 09:07:00",
               "2018-06-02 11:37:00", "2018-06-02 14:37:00",
               "2018-06-02 16:37:00"))


## Test countValidPlot
expect_error(countValidPlot(1),
             "TP should be an object of class TP")
expect_error(countValidPlot(testTP, trait = 1),
             "trait should be a character string of length one")
expect_error(countValidPlot(testTP, trait = "a"),
             "a should be a column in TP")
expect_error(countValidPlot(testTP, trait = "t1", plotIds = 1),
             "plotIds should be NULL or a character vector")
expect_error(countValidPlot(testTP, trait = "t1", plotIds = "a"),
             "All plotIds should be in TP")
valPlotCount <- countValidPlot(testTP, trait = "t1")
expect_equivalent(valPlotCount, c(rep(5, times = 7), 0, rep(5, times = 17)))
expect_equal(names(valPlotCount), unique(testDat$pos))
