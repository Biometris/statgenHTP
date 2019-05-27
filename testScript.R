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

basefunction(inTD[1:2], trait = "pheno",
             covariates = c("Sowing_Block", "Image_pos"),
             out1 = "BLUPs_PAM_modRep.csv",
             out2 = "Corrected_PAM_modRep.csv")

basefunction(inTD[1:2], trait = "pheno",
             covariates = c("Sowing_Block", "Image_pos"),
             useCheck = TRUE,
             out1 = "BLUPs_PAM_modRep_Check.csv",
             out2 = "Corrected_PAM_modRep_Check.csv")


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

inTD2 <- createTD(dat = inDat2, genotype = "TrtGeno",
                  timePoint = "time1", time = "Date", plotId = "pos",
                  rowNum = "Position", colNum = "Line")
basefunction(inTD2[1:2], trait = "LA_Estimated", covariates = "TrtPop",
             geno.decomp = "TrtPop",
             out1 = "BLUPs_ZA17_LeafArea.csv",
             out2 = "Corrected_ZA17_LeafArea.csv")
