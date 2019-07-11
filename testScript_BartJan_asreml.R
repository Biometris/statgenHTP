timePoint = 10

inDat <- data.table::fread("./data-raw/Original_PAM_reshape.csv",
                           data.table = FALSE)
# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
inDat$ID <- interaction(inDat[["x"]], inDat[["y"]], sep = "_")
inDat <- inDat[!is.na(inDat$pheno), ]
# Create an indicator for each plot (according to the row and column position)
inDat$pos <- paste0("c", inDat[["x"]], "r", inDat[["y"]])
# I removed a plant that has very few measurements
inDat <- inDat[inDat$pos != "c1r54",]
inTP <- createTimePoints(dat = inDat, genotype = "Genotype",
                         timePoint = "timepoints",
                         repId = "Sowing_Block",
                         plotId = "pos",
                         rowNum = "y", colNum = "x",
                         addCheck = TRUE,
                         checkGenotypes = c("col", "ely", "evo1", "ler"))

fitMods1a <- fitModels(TP = inTP[timePoint], trait = "pheno",
                       covariates = c("repId", "Image_pos"))
a1 <- predict(fitMods1a[[1]], which = c("repId", "Image_pos"))

fitMods1b <- fitModels(TP = inTP[timePoint], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml")
fitMod <- fitMods1b[[1]]
a2 <- predict(fitMod, classify = "repId + Image_pos")$pvals

aTot <- merge(a1, a2, by = c("repId", "Image_pos"))

max(abs(aTot$predicted.value - aTot$predicted.values))
max(abs(aTot$standard.errors - aTot$std.error))

cor(aTot$predicted.value, aTot$predicted.values)
plot(aTot$predicted.value, aTot$predicted.values)
abline(a = 0, b = 1, col = "red")


fitMods1c <- fitModels(TP = inTP[timePoint], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml",
                       spatial = TRUE)
fitMod2 <- fitMods1c[[1]]
a3 <- predict(fitMod2, classify = "repId + Image_pos + rowId + colId",
              present = c("rowId", "colId", "repId", "Image_pos"))$pvals

a0 <- predict(fitMod2, classify = "(Intercept)",
              present = c("rowId", "colId", "repId", "Image_pos"))$pvals

a4 <- predict(fitMods1a[[1]], which = c("rowId", "colId", "rowNum", "colNum",
                                        "repId", "Image_pos"))

aTot2 <- merge(a3, a4, by =  c("rowId", "colId", "repId", "Image_pos"))

max(abs(aTot2$predicted.value - aTot2$predicted.values), na.rm = TRUE)
max(abs(aTot2$standard.errors - aTot2$std.error), na.rm = TRUE)

cor(aTot2$predicted.value, aTot2$predicted.values, use = "pairwise.complete.obs")
plot(aTot2$predicted.value, aTot2$predicted.values)
abline(a = 0, b = 1, col = "red")

test <- merge(fitMod2$call$data, a3, by = c("rowId", "colId", "repId", "Image_pos"))
test$newTrait <- test$pheno - test$predicted.value +
  fitMod2$coefficients$fixed["(Intercept)",]
# fitMods1a[[1]]$coeff[["Intercept"]]
# a0$predicted.value

corSpat <- statgenHTP:::correctSpatial(fitMods1a[[1]])

comp <- merge(test, corSpat, by = c("rowId", "colId"))
mean(abs(comp$newTrait.x - comp$newTrait.y))
cor(comp$newTrait.x, comp$newTrait.y)
plot(comp$newTrait.x, comp$newTrait.y)
abline(a = 0, b = 1, col = "red")
