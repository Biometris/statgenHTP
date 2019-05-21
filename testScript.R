folder <- "C:/WUR/statgenHTP/"
setwd(paste0(folder, "output/"))

inDat <- data.table::fread("../rawdata/Original_PAM_reshape.csv",
                               data.table = FALSE)

# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
inDat$ID <- interaction(inDat[["x"]], inDat[["y"]], sep = "_")
inDat <- inDat[!is.na(inDat$pheno), ]
# modification to be able to include the check in the model
# I forgot to mention that: quite often in the platform experiments,
# there is a highly replicated genotype (or 2, 3 ,4) that we can use in the model
# but we need to format the dataset
inDat$Check <- ifelse(inDat$Genotype %in% c("col", "ely", "evo1", "ler"),
                      inDat$Genotype, "_pop")
inDat$Genobis <- inDat$Genotype
inDat$Genobis[inDat$Check != "_pop"] <- NA
# Create an indicator for each plot (according to the row and column position)
inDat$pos <- paste0("c", inDat$x, "r", inDat$y)
# Create factors
inDat$Image_pos = as.factor(inDat$Image_pos)
inDat$Sowing_Block = as.factor(inDat$Sowing_Block)

# I removed a plant that has very few measurements
inDat <- inDat[inDat$pos != "c1r54",]
inDat <- droplevels(inDat)

inTD <- createTD(dat = inDat, genotype = "Genotype",
                 timePoint = "timepoints",  plotId = "ID", rowNum = "y",
                 colNum = "x")


inTD2 <- inTD[1:2]
basefunction(inTD2)
