# rm(list=ls())
setwd(paste0("~/Documents/PostDoc/Pipeline_Biometris/statgenHTP/"))
# devtools::load_all(".")

#_________________ EXAMPLE 2 _____________________________________________________________

#### DATASET-SPECIFIC FORMATING #######

inDat2 <- data.table::fread("./data-raw/Data_modif_ZA17_anonymous.csv",
                            data.table = FALSE)
# Create an indicator for each plot (according to the row and column position)
inDat2$pos <- paste0("c", inDat2$Line, "r", inDat2$Position)

#### BEGINING FUNCTION USE #######

inTP2 <- createTimePoints(dat = inDat2,
                          genotype = "geno",
                          timePoint = "Date",
                          plotId = "pos",
                          rowNum = "Position",
                          colNum = "Line") #,
                          # repId = "Rep")

fitMods2c <- fitModels(TP = inTP2[20],
                       trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml",
                       spatial = TRUE)

fitMods2c[[1]]$converge


#### TRY TO GET CORRECTED DATA #######

trait = "LA_Estimated"
geno.decomp = c("Scenario", "population")
don <- inTP2[20][[1]]
#
don[["genotype.original"]] <- don[["genotype"]]
don[["geno.decomp"]] <- interaction(don[geno.decomp], sep = "_")
don[["genotype"]] <- interaction(don[["geno.decomp"]], don[["genotype"]], sep = "_")

# Test with SpATS
library(SpATS)

nseg = c(nlevels(don[["colId"]]), nlevels(don[["rowId"]]))

modTestSpats <- SpATS(response = trait,
                          fixed = ~  geno.decomp ,
                          random = ~ colId + rowId ,
                          spatial = ~ SpATS::PSANOVA(colNum, rowNum,
                                                     nseg = nseg,
                                                     nest.div=c(2,2)),
                          genotype = "genotype",
                          genotype.as.random = TRUE,
                          geno.decomp = "geno.decomp",
                          data = don,
                          control = list(maxit = 50,
                                         tolerance = 1e-03,
                                         monitoring = 0))

modTestSpats$coeff["Intercept"]

predSpats <- predict(modTestSpats,
                     which = c("colId","rowId","colNum","rowNum")) #
hist(predSpats$predicted.values)

# Test with ASReml-R

library(asreml)

TPTab <- as.data.frame(table(don[["colId"]], don[["rowId"]]))
TPTab <- TPTab[TPTab$Freq == 0, , drop = FALSE]

if (nrow(TPTab) > 0) {
  extObs <- setNames(as.data.frame(matrix(nrow = nrow(TPTab),
                                          ncol = ncol(don))),
                     colnames(don))
  extObs[["timePoint"]] <- don[["timePoint"]][1]
  extObs[, c("colId", "rowId")] <- TPTab[, c("Var1", "Var2")]
  extObs[, c("colNum", "rowNum")] <-
    c(as.numeric(levels(TPTab[, "Var1"]))[TPTab[, "Var1"]],
      as.numeric(levels(TPTab[, "Var2"]))[TPTab[, "Var2"]])
  don <- rbind(don, extObs)
}
don <- don[order(don[["rowId"]], don[["colId"]]), ]

modTestAsreml <- asreml(fixed = LA_Estimated ~ geno.decomp,
                  random = ~ rowId + colId + at(geno.decomp):genotype, # + units
                  resid = ~ ar1(rowId):ar1(colId),
                  data = don,
                  trace = FALSE, maxiter = 200,
                  na.action = na.method(x = "include"))

modTestAsreml$mf
ar1(rowId):ar1(colId)
predAsreml <- data.frame(predict(modTestAsreml, classify = "rowId + colId")$pvals)
hist(predAsreml$predicted.value)

predSpats$key <- interaction(predSpats[["rowId"]], predSpats[["colId"]], sep = "_")
predAsreml$key <- interaction(predAsreml[["rowId"]], predAsreml[["colId"]], sep = "_")

predSpats$pred.spats <- predSpats$predicted.values
predSpats$pred.asreml <- predAsreml$predicted.value[match(predSpats$key,predAsreml$key)]

plot(predSpats$pred.spats~predSpats$pred.asreml,pch=16,cex=1.2)
abline(0,1,col="red")




