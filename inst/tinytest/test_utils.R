
tpDat <- data.frame(rowNum = rep(1:4, each = 4),
                    colNum = rep(1:4, times = 4),
                    repId = rep(1:2, each = 8))
bord <- statgenHTP:::calcPlotBorders(tpDat = tpDat, bordVar = "repId")

expect_true(inherits(bord, "list"))
expect_equal(names(bord), c("horW", "vertW"))
expect_true(inherits(bord$horW, "data.frame"))
expect_true(inherits(bord$vertW, "data.frame"))
expect_equal(bord$horW[["x"]], rep(1:4, each = 3))
expect_equal(bord$horW[["y"]], rep(c(1, 3, 5), times = 4))
expect_equal(bord$vertW[["x"]], rep(c(1, 5), times = 4))
expect_equal(bord$vertW[["y"]], rep(1:4, each = 2))

tpDat <- tpDat[-1, ]
bord <- statgenHTP:::calcPlotBorders(tpDat = tpDat, bordVar = "repId")
expect_equal(bord$horW[["x"]], c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4))
expect_equal(bord$horW[["y"]], c(3, 5, 1, 3, 5, 1, 3, 5, 1, 3, 5))
expect_equal(bord$vertW[["x"]], c(5, 1, 5, 1, 5, 1, 5))
expect_equal(bord$vertW[["y"]], c(1, 2, 2, 3, 3, 4, 4))


