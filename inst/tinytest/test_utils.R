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


