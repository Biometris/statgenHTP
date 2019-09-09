### Test calcPlotBorders.
tpDat <- data.frame(rowNum = rep(1:4, each = 4),
                    colNum = rep(1:4, times = 4),
                    repId = rep(1:2, each = 8))

## No missing values. Borders between replicates and around edge.
bord <- statgenHTP:::calcPlotBorders(tpDat = tpDat, bordVar = "repId")
expect_true(inherits(bord, "list"))
expect_equal(names(bord), c("horW", "vertW"))
expect_true(inherits(bord$horW, "data.frame"))
expect_true(inherits(bord$vertW, "data.frame"))
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
                  id = c("p1", "p2", "p1"), trait = 1:3)

## Check that single missing value is added when no extra columns present.
dfOut1 <- statgenHTP:::addMissVals(df1, "trait")
expect_true(inherits(dfOut1, "data.frame"))
expect_true(inherits(dfOut1[["timePoint"]], "POSIXct"))
expect_equal(nrow(dfOut1), 4)
expect_equivalent(unlist(dfOut1[4, ]), c(2, 1567382400, NA))

## Check that multiple missing values are added when no extra columns present.
df2 <- df1[-1, ]
dfOut2 <- statgenHTP:::addMissVals(df2, "trait")
expect_equal(nrow(dfOut2), 4)
expect_equivalent(unlist(dfOut2[1, ]), c(1, 1567296000, NA))
expect_equivalent(unlist(dfOut2[4, ]), c(2, 1567382400, NA))

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


