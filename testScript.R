folder <- "C:/Projects/R packages/statgenhtp/"
setwd(paste0(folder, "output/"))

inDat <- data.table::fread("../rawdata/Original_PAM_reshape.csv",
                               data.table = FALSE)

# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
inDat$ID <- interaction(inDat[["x"]], inDat[["y"]], sep = "_")
inDat <- inDat[!is.na(inDat$pheno), ]
# Create an indicator for each plot (according to the row and column position)
inDat$pos <- paste0("c", inDat$x, "r", inDat$y)
# Create factors
inDat$Image_pos = as.factor(inDat$Image_pos)
inDat$Sowing_Block = as.factor(inDat$Sowing_Block)

# I removed a plant that has very few measurements
inDat <- inDat[inDat$pos != "c1r54",]
inDat <- droplevels(inDat)

inTD <- createTD(dat = inDat, genotype = "Genotype",
                 timePoint = "timepoints", plotId = "pos", rowNum = "y",
                 colNum = "x", addCheck = TRUE,
                 checkGenotypes = c("col", "ely", "evo1", "ler"))

plot(inTD, plotType = "layout",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     highlight = c("col", "ely", "evo1", "ler"))
plot(inTD, plotType = "cor", traits = "pheno")
plot(inTD, plotType = "box", traits = "pheno",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     colorBy = "Sowing_Block")

pdf("Phenovator_Rene_raw_data_na.pdf", height = 8, width = 12)
plot(inTD, plotType = "raw", traits = "pheno")
dev.off()

fitMods <- fitModels(TD = inTD, trait = "pheno",
                     covariates = c("Sowing_Block", "Image_pos"))

BLUPs <- getBLUPs(fitMods, outFile = "BLUPs_PAM_modRep.csv")
spatCorr <- getCorrected(fitMods, outFile = "Corrected_PAM_modRep.csv")
variance <- getVar(fitMods)
h2 <- getHerit(fitMods)

plot(fitMods, plotType = "corrPred")
plot(fitMods, plotType = "rawPred", traits = "pheno")
plot(fitMods, plotType = "herit")
plot(fitMods, plotType = "variance")
plot(fitMods, plotType = "effDim")


## Second example
inDat2 <- data.table::fread("../rawdata/Data_modif_ZA17_anonymous.csv",
                           data.table = FALSE)

# Create an indicator for each plot (according to the row and column position)
inDat2$pos <- paste0("c", inDat2$Line, "r", inDat2$Position)
# Create factors
inDat2$Treatment = as.factor(inDat2$Scenario)
inDat2$Population = as.factor(inDat2$population)
### This part is the columns that should be created to run SpATS with
#  a factor that split the genotypic variance
# here there is the combination of genotypic panel and water treatment
inDat2$TrtGeno <- interaction(inDat2$Treatment, inDat2$geno, sep = "_")
inDat2$TrtPop = interaction(inDat2$Treatment, inDat2$Population, sep = "_")

inTD2 <- createTD(dat = inDat2, genotype = "TrtGeno", timePoint = "Date",
                  plotId = "pos", rowNum = "Position", colNum = "Line")

plot(inTD2, plotType = "layout", timePoints = "2017-04-13")
plot(inTD2, plotType = "cor", traits = "LA_Estimated")
plot(inTD2, plotType = "box", traits = "LA_Estimated",
     timePoints = names(inTD2[1:3]))
pdf("Phenovator_ZA17_raw_data_na.pdf", height = 8, width = 12)
plot(inTD2, plotType = "raw", traits = "LA_Estimated")
dev.off()

basefunction(inTD2[1:2], trait = "LA_Estimated", covariates = "TrtPop",
             geno.decomp = "TrtPop",
             out1 = "BLUPs_ZA17_LeafArea.csv",
             out2 = "Corrected_ZA17_LeafArea.csv")
