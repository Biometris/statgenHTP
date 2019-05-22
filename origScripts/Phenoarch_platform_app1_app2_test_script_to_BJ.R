# This file contains the code for the two approaches (Coté Rodriguez)
# for the analyses of the ZA17 platform (Phenoarch)
# EM, 01/04/2019


# Clean the workspace
rm(list = ls())

# Load the needed packages
library(SpATS)
library(ggplot2)
# library(fields)
library(dplyr)
library(plyr)
library(lattice)
# library(latticeExtra)
# library(Matrix)
library(data.table)


folder <- "~/Documents/PostDoc/EPPN2020/Platform/M3P/Phenoarch/ZA17/"
setwd(folder)

dat.modif <- read.table("../rawdata/Data_modif_ZA17_anonymous.csv",sep=',',h=T)
dat.modif$ID <- dat.modif$XY
dat.modif$timepoint <- lubridate::ymd_hms(dat.modif$time1)


###################################################################################
######## Settings and output formatting

# In these analyses, the genotype is always included as random
	genotype.as.random = TRUE

# Trait to analyse
	trait = "LA_Estimated" # "Biomass_Estimated"

	# Create factors
	dat.modif$Col = as.factor(dat.modif$Line)
	dat.modif$Row = as.factor(dat.modif$Position)
	dat.modif$Genotype = as.factor(dat.modif$geno)
	dat.modif$Treatment = as.factor(dat.modif$Scenario)
	dat.modif$Population = as.factor(dat.modif$population)

	### This part is the columns that should be created to run SpATS with
	#  a factor that split the genotypic variance
	# here there is the combination of genotypic panel and water treatment
	dat.modif$TrtGeno = as.factor(paste0(dat.modif$Treatment,"_",dat.modif$Genotype))
	dat.modif$TrtPop = as.factor(paste0(dat.modif$Treatment,"_",dat.modif$Population))
	dat.modif$RepTrtPop = as.factor(paste0(dat.modif$Rep,"_",dat.modif$Treatment,"_",dat.modif$Population))
	# Create num
	dat.modif$Colnum = as.numeric(dat.modif$Line)
	dat.modif$Rownum = as.numeric(dat.modif$Position)

	# Transform the timepont to a numeric variable
	dat.modif$Time <- as.numeric(plyr::mapvalues(dat.modif$Date ,
	                                       from = levels(factor(dat.modif$Date )),
	                                       to = 1:length(levels(factor(dat.modif$Date )))))


	# Time points (factor codification)
	times.fact <- levels(factor(dat.modif$Date))

	# Time points (numeric codification)
	times.num <- sort(unique(dat.modif$Time))

# Create an indicator for each plot (according to the row and column position)
	dat.modif$pos <- paste0("c",dat.modif$Col,"r",dat.modif$Row, sep = "")

# Objects to save results for Approach 1
	Geno_pred <- Col_pred <- Row_pred <- NULL			# Predictions (may contain the intercept as well as other factors)
	Geno_BLUPs <- Col_BLUPs <- Row_BLUPs <- NULL		# BLUPs (only the estimated coefficients - they do not contain the intercept or other factors)

	ed_row <- ed_col <- ed_surface <- sigma2 <- sigma2_row <- sigma2_col  <- vector(length = length(times.fact))

	# Since we are using the geno.decomp argument, and we have 4 "populations"
	# (combination of treatment and population), we have four different heritabilities.
	h2 <- sigma2_geno <- matrix(0, ncol = 4, nrow = length(times.fact))
	colnames(h2) <- paste0("TrtPop",levels(dat.modif$TrtPop))

# Object to save results for Approach 2
	data_corr_spat <- NULL

################################################################################
######## Loop on timepoint to run SpATS

#for (ti in 1:length(times.fact)) {
for (ti in 1:2) {
	# ti = 18
	cat(times.fact[ti],'\n')

	### subset of dataset:
		dat.ti <- dat.modif[dat.modif$Date==times.fact[ti],]
		dat.ti <- droplevels(dat.ti)

	### number of segments for SpATS:
		nseg = c(nlevels(dat.ti$Col),
           	nlevels(dat.ti$Row))

	### Col and row (needed to obtain the BLUPs)
	Col_levels <- paste("Col", levels(dat.ti$Col), sep = "")
	Row_levels <- paste("Row", levels(dat.ti$Row), sep = "")

 # Fit the model
	fit.SpATS <- SpATS::SpATS(response = trait,
                              fixed = ~  TrtPop,
                              random = ~ Col + Row ,
                              spatial = ~ SpATS::PSANOVA(Colnum, Rownum,
                                                         nseg = nseg,
                                                         nest.div=c(2,2)),
                              genotype = "TrtGeno",
                              genotype.as.random = TRUE,
                              geno.decomp = "TrtPop",
                              data = dat.ti,
                              control = list(maxit = 50,
                                             tolerance = 1e-03,
                                             monitoring = 0))

	# plot(fit.SpATS)
	##########################################################################################
	# Approach 1: obtain the genotypic predictions
	##########################################################################################
		##############################
		# Genotype predictions
		##############################
		# Genotype prediction (including the effect of TrtPop, as well as the intercept)
		p_geno <- predict(fit.SpATS, which = c("TrtGeno","TrtPop"))

		# Include time point
		pred_g  <-  mutate(p_geno, Time = times.num[ti])

		# Select the needed variables for subsequent analyses
		pred_g  <- select(pred_g, Time, TrtGeno, predicted.values, standard.errors)

		# Add results
		Geno_pred = rbind(Geno_pred, pred_g)

		##############################
		# Col predictions
		##############################
 		# Col prediction (including intercept)
	  	p_Col <- predict(fit.SpATS, which = "Col")

	  	# Include time point
  		pred_c <- dplyr:::mutate(p_Col, Time = times.num[ti])

  		# Select the needed variables for subsequent analyses
  		pred_c <- dplyr:::select(pred_c, Time, Col, predicted.values, standard.errors)

  		# Add results
  		Col_pred = rbind(Col_pred, pred_c)

  		##############################
		# Row predictions
		##############################
		# Row prediction (including intercept)
	  	p_Row <- predict(fit.SpATS, which = "Row")

	  	# Include time point
  		pred_r <- dplyr:::mutate(p_Row, Time = times.num[ti])

  		# Select the needed variables for subsequent analyses
  		pred_r <- dplyr:::select(pred_r, Time, Row, predicted.values, standard.errors)

  		# Add results
  		Row_pred = rbind(Row_pred, pred_r)

  		##############################
		# Genotype BLUPs
		##############################
		# BLUPs
		BLUPs_geno <- fit.SpATS$coeff[fit.SpATS$terms$geno$geno_names]

		# Standard error
		se_BLUPs_geno <- sqrt(diag(fit.SpATS$vcov$C11_inv[fit.SpATS$terms$geno$geno_names,fit.SpATS$terms$geno$geno_names]))

		# Create data.frame the needed variables for subsequent analyses
		df_geno <- data.frame(TrtGeno = names(BLUPs_geno), predicted.values = BLUPs_geno, standard.errors = se_BLUPs_geno, Time = times.num[ti])

		# Add results
		Geno_BLUPs = rbind(Geno_BLUPs, df_geno)

		##############################
		# Col BLUPs
		##############################
		# BLUPs
  		BLUPs_Col <- fit.SpATS$coeff[Col_levels]

  		# Standard error
  		se_BLUPs_Col <- sqrt(diag(fit.SpATS$vcov$C22_inv[Col_levels,Col_levels]))

  		# Create data.frame the needed variables for subsequent analyses
  		df_Col <- data.frame(Col = names(BLUPs_Col), predicted.values = BLUPs_Col, standard.errors = se_BLUPs_Col, Time = times.num[ti])

  		# Add results
  		Col_BLUPs <- rbind(Col_BLUPs, df_Col)


  		##############################
		# Row BLUPs
		##############################
		# BLUPs
		BLUPs_Row <- fit.SpATS$coeff[Row_levels]
		se_BLUPs_Row <- sqrt(diag(fit.SpATS$vcov$C22_inv[Row_levels,Row_levels]))
		df_Row <- data.frame(Row = names(BLUPs_Row), predicted.values = BLUPs_Row, standard.errors = se_BLUPs_Row, Time = times.num[ti])
		Row_BLUPs <- rbind(Row_BLUPs, df_Row)


		##############################
		# Heritabilities
		##############################
		h2.aux <- getHeritability(fit.SpATS)

		### not optimised loop, sorry...
		if(length(h2.aux) == 4) {
		  h2[ti,] <- h2.aux
		} else {
		  aux <- rep(NA, 4)
		  aux[match(names(h2.aux), colnames(h2))] <- h2.aux
		  h2[ti,] <- aux
		}

 		##############################
		# Other information
		##############################
		sigma2[ti] <- fit.SpATS$psi[1]
		ed_col[ti] <- fit.SpATS$eff.dim["Col"]
		ed_row[ti] <- fit.SpATS$eff.dim["Row"]
		ed_surface[ti] <- sum(fit.SpATS$eff.dim[c("f(Colnum)", "f(Rownum)", "f(Colnum):Rownum", "Colnum:f(Rownum)","f(Colnum):f(Rownum)")])
		sigma2_col[ti] <- fit.SpATS$var.comp["Col"]
		sigma2_row[ti] <- fit.SpATS$var.comp["Row"]

	##########################################################################################
	# Approach 2: correct for spatial effects and other unnecesary factors
	##########################################################################################

	# Include in the prediction the factors (variables) whose effect we are interested in removing
		p_a2 <- predict(fit.SpATS, which = c("Colnum","Rownum","Col","Row"))

	# Order the data and prediction according to the previous covariates
		data.ord <-  dat.ti[order(dat.ti$Colnum, dat.ti$Rownum),]
		p_a2 <- p_a2[order(p_a2$TrtPop, p_a2$Colnum, p_a2$Rownum),]

	# Add to the predictions some needed information
		p_a2$TrtGeno <- data.ord[,"TrtGeno"]
		p_a2$TrtPop <- data.ord[,"TrtPop"]
		p_a2$Treatment <- data.ord[,"Treatment"]
		p_a2$pos <- data.ord[,"pos"]

	# Obtain the corrected trait
		p_a2$new_trait <- data.ord[,trait] - p_a2$predicted.values + fit.SpATS$coeff["Intercept"]

	# Include time point
		pred_a2 <- dplyr:::mutate(p_a2, Time = times.num[ti])

	# Select the needed variables for subsequent analyses
		pred_a2 <- dplyr:::select(pred_a2, new_trait, TrtGeno, TrtPop, Treatment, Time, Colnum, Rownum, Col, Row, pos)

	# Add results
		data_corr_spat <- rbind(data_corr_spat, pred_a2)

}

write.table(data_corr_spat,"../outputOrigScript/Corrected_ZA17_LeafArea.csv",sep=",",row.names=F)
write.table(Geno_pred,"../outputOrigScript/BULPs_ZA17_LeafArea.csv",sep=",",row.names=F)


#############################################################
# Some results (all genotypes)
#############################################################
	# Raw data
		# NOTE: in this case xyplot connect the lines even when there are missing data. Thus, this information is lost.
		pdf("New_Cote_et_Pspline/Phenoarch_ZA17_LeafArea_raw_data.pdf")
		xyplot(get(trait) ~ Time|TrtGeno,
		       data = dat.modif,
		       groups = dat.modif$pos,
		       type = c("l"),
		       xlab="Time", ylab = "Trait",
		       main = "Phenoarch ZA17 - Leaf Area \n Raw data",
		       layout = c(5,5))
		dev.off()

		# NOTE: With the following code, we add explicitly the rows of the data sets that are missing, thus they will appear when plotting
		# Very time consuming!!!
			# There are missing data. Fill in the gaps
			xy.coord <- data.table(expand.grid(Time = sort(unique(dat.modif$Time)), pos = unique(dat.modif$pos)))
			setkeyv(xy.coord, c("Time", "pos"))
			dat.modif <- data.table(dat.modif)
			setkeyv(dat.modif, c("Time", "pos"))
			dat.modif_na <- dat.modif[xy.coord]

			# Add the information we have for the missing values
			for(i in unique(dat.modif_na$pos)) {
				dat.modif_na[dat.modif_na$pos == i, "TrtGeno"] <- na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "TrtGeno"]))
				dat.modif_na[dat.modif_na$pos == i, "TrtPop"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "TrtPop"]))
				dat.modif_na[dat.modif_na$pos == i, "Treatment"] <- na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Treatment"]))
				dat.modif_na[dat.modif_na$pos == i, "Colnum"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Colnum"]))
				dat.modif_na[dat.modif_na$pos == i, "Rownum"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Rownum"]))
				dat.modif_na[dat.modif_na$pos == i, "Col"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Col"]))
				dat.modif_na[dat.modif_na$pos == i, "Row"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Row"]))
			}

			write.table(dat.modif_na,"New_Cote_et_Pspline/dat_modif_na_LeafArea.csv",sep=",",row.names=F)

			pdf("New_Cote_et_Pspline/Phenoarch_ZA17_LeafArea_raw_data_na.pdf")
			xyplot(get(trait) ~ Time|TrtGeno,
			       data = dat.modif_na,
			       groups = dat.modif_na$pos,
			       type = c("l"),
			       xlab="Time", ylab = "Trait",
			       main = "Phenoarch ZA17 - Leaf Area \n Raw data",
			       layout = c(5,5))
			dev.off()

	###
	### Raw data + genotypic predictions (Approach 1)
	###

		pdf("New_Cote_et_Pspline/Phenoarch_ZA17_LeafArea_raw_data_geno_pred_app1.pdf")
		aa <- xyplot(get(trait) ~ Time|TrtGeno,
		             data = dat.modif,
		             groups = dat.modif$pos,
								type = c("l"),
								xlab="Time", ylab = "Trait",
								main = "Phenoarch ZA17 - Leaf Area \n Raw data + Genotypic predictions (App 1)",
								layout = c(5,5))
		bb <- xyplot(predicted.values ~ Time|TrtGeno,
		             data = Geno_pred,
									xlab="", ylab = "", main = "",
									type = c("l"), col = "black", lwd = 2,
									layout = c(5,5))
		aa + as.layer(bb)
	 dev.off()

	 ###
	 ### Corrected trait (Approach 2) + + genotypic predictions (Approach 1)
	 ###

		pdf("New_Cote_et_Pspline/Phenoarch_platform_corrected_trait_app2_geno_pred_app1.pdf")
		aa <- xyplot(new_trait ~ Time|TrtGeno,
		             data = data_corr_spat,
		             groups = data_corr_spat$pos,
							type = c("l"),
							xlab="Time", ylab = "Corrected trait",
							main = "Phenoarch platform \n Corrected trait (App 2) + Genotypic predictions (App 1)",
							layout = c(5,5))
		bb <- xyplot(predicted.values ~ Time|TrtGeno,
		             data = Geno_pred,
							xlab="", ylab = "", main = "",
							type = c("l"), col = "black", lwd = 2,
							layout = c(5,5))
		aa + as.layer(bb)
		dev.off()

	# Heritabilities (Approach 1)
		pdf("New_Cote_et_Pspline/Phenoarch_platform_heritability.pdf")
		op <- par(mfrow = c(2,2))
		plot(h2[,1], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[1])
		lines(h2[,1])
		plot(h2[,2], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[2])
		lines(h2[,2])
		plot(h2[,3], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[3])
		lines(h2[,3])
		plot(h2[,4], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[4])
		lines(h2[,4])
		par(op)
		dev.off()

		# Effective dimensions (Approach 1)
		pdf("New_Cote_et_Pspline/Phenoarch_ZA17_ed_app1.pdf")
		par(mfrow = c(2,2))
		plot(times.num, ed_surface, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED spatial surface", cex.lab=1.5, type = "l")
		points(times.num, ed_surface)

		plot(times.num, ed_row, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED row RE", cex.lab=1.5, type = "l")
		points(times.num, ed_row)

		plot(times.num, ed_col, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED column RE", cex.lab=1.5, type = "l")
		points(times.num, ed_col)
		dev.off()

		# Variances (Approach 1)
		pdf("New_Cote_et_Pspline/Phenoarch_ZA17_LeafArea_var_app1.pdf")
		par(mfrow = c(2,2))
		plot(times.num, sigma2, cex.axis=1.5, xlab="Time", ylab = expression(sigma^2), main = "Residual variance", cex.lab=1.5, type = "l")
		points(times.num, sigma2)

		plot(times.num, sigma2_col, cex.axis=1.5,  xlab="Time", ylab = expression(sigma^2), main = "Column RE variance", cex.lab=1.5, type = "l")
		points(times.num, sigma2_col)

		plot(times.num, sigma2_row, cex.axis=1.5, xlab="Time", ylab = expression(sigma^2), main = "Row RE variance", cex.lab=1.5, type = "l")
		points(times.num, sigma2_row)
		dev.off()

		# Row and Column (BLUPs and predictions)
		pdf("New_Cote_et_Pspline/Phenoarch_ZA17_LeafArea_RowCol_pred_BLUPs.pdf",6,6)
		# Column
		xyplot(predicted.values ~ Time,  data = Col_pred, groups = Col_pred$Col,
		       xlab="Time", ylab = "Column random factor prediction",
		       main = "Phenoarch ZA17 - LeafArea \n Column predictions", type = c("g", "p", "o"))

		xyplot(predicted.values ~ Time,  data = Col_BLUPs, groups = Col_BLUPs$Col,
		       xlab="Time", ylab = "Column random factor BLUPs",
		       main = "Phenoarch ZA17 - LeafArea \n Column BLUPs", type = c("g", "p", "o"))

		# Row
		xyplot(predicted.values ~ Time,  data = Row_pred, groups = Row_pred$Row,
		       xlab="Time", ylab = "Row random factor prediction",
		       main = "Phenoarch ZA17 - LeafArea \n Row predictions", type = c("g", "p", "o"))

		xyplot(predicted.values ~ Time,  data = Row_BLUPs, groups = Row_BLUPs$Row,
		       xlab="Time", ylab = "Row random factor BLUPs",
		       main = "Phenoarch ZA17 - LeafArea \n Row BLUPs", type = c("g", "p", "o"))
		dev.off()


#############################################################
# Some results (some genotypes)
#############################################################
	sel.geno <- c("WD.GenoA30","WD.GenoA21","WW.GenoA30","WW.GenoA21","WW.GenoA24","WD.GenoA24", "WW.GenoB10", "WD.GenoB10")
	# Raw data, selected genotypes
		dat.subset <- dat.modif[dat.modif$TrtGeno %in% sel.geno,]
		dat.subset$TrtGeno <- droplevels(dat.subset$TrtGeno)
		dat.subset$TrtPop <- droplevels(dat.subset$TrtPop)

	# Genotype predictions (Approach 1), selected genotypes
		Geno_pred.subset <- Geno_pred[Geno_pred$TrtGeno %in% sel.geno,]
		Geno_pred.subset$TrtGeno <- droplevels(Geno_pred.subset$TrtGeno)

	# Corrected trait, selected genotypes
		data_corr_spat.subset <- data_corr_spat[data_corr_spat$TrtGeno %in% sel.geno,]
		data_corr_spat.subset$TrtGeno <- droplevels(data_corr_spat.subset$TrtGeno)

	# Raw data
		#pdf("Phenoarch_platform_Raw_data_selected_geno.pdf", height = 8, width = 12)
		xyplot(Trait ~ Time|TrtGeno,  data = dat.subset, groups = dat.subset$pos,
								type = c("l"),
								xlab="Time", ylab = "Trait", main = "Phenoarch platform \n Raw data", layout = c(4,2))
		#dev.off()

	# Raw data + genotypic predictions (Approach 1)
		#pdf("Phenoarch_platform_raw_data_geno_pred_app1_selected_geno.pdf", height = 8, width = 12)

		aa <- xyplot(Trait ~ Time|TrtGeno,  data = dat.subset, groups = dat.subset$pos,
									type = c("l"),
									xlab="Time", ylab = "Trait", main = "Phenoarch platform \n Raw data + Genotypic predictions (App 1)", layout = c(4,2))
		bb <- xyplot(predicted.values ~ Time|TrtGeno,  data = Geno_pred.subset,
									xlab="", ylab = "", main = "", type = c("l"), col = "black", lwd = 2)
		aa + as.layer(bb)

		#dev.off()

	# Corrected trait (Approach 2) + genotypic predictions (Approach 1)

		# NOTE: in this case xyplot connect the lines even when there are missing data. Thus, this information is lost.

		#pdf("Phenoarch_platform_corrected_trait_app2_geno_pred_app1_selected_geno", height = 8, width = 12)
		aa <- xyplot(new_trait ~ Time|TrtGeno,  data = data_corr_spat.subset, groups = data_corr_spat.subset$pos,
							type = c("l"),
							xlab="Time", ylab = "Corrected trait", main = "Phenoarch platform data \n Corrected trait (App 2) + Genotypic predictions (App 1)", layout = c(4,2))
		bb <- xyplot(predicted.values ~ Time|TrtGeno,  data = Geno_pred.subset,
							xlab="", ylab = "", main = "", type = c("l"), col = "black", lwd = 2)
		aa + as.layer(bb)
		#dev.off()

		# NOTE: With the following code, we add explicitly the rows of the data sets that are missing, thus they will appear when plotting
		# Very time consuming!!!
			# There are missing data. Fill in the gaps
			xy.coord <- data.table(expand.grid(Time = sort(unique(data_corr_spat$Time)), pos = unique(data_corr_spat$pos)))
			setkeyv(xy.coord, c("Time", "pos"))
			data_corr_spat <- data.table(data_corr_spat)
			setkeyv(data_corr_spat, c("Time", "pos"))
			data_corr_spat_na <- data_corr_spat[xy.coord]

			# Add the information we have for the missing values
			for(i in unique(data_corr_spat_na$pos)) {
				data_corr_spat_na[data_corr_spat_na$pos == i, "TrtGeno"] <- na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "TrtGeno"]))
				data_corr_spat_na[data_corr_spat_na$pos == i, "TrtPop"] <- 	na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "TrtPop"]))
				data_corr_spat_na[data_corr_spat_na$pos == i, "Treatment"] <- na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "Treatment"]))
				data_corr_spat_na[data_corr_spat_na$pos == i, "Colnum"] <- 	na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "Colnum"]))
				data_corr_spat_na[data_corr_spat_na$pos == i, "Rownum"] <- 	na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "Rownum"]))
				data_corr_spat_na[data_corr_spat_na$pos == i, "Col"] <- 	na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "Col"]))
				data_corr_spat_na[data_corr_spat_na$pos == i, "Row"] <- 	na.omit(unique(data_corr_spat_na[data_corr_spat_na$pos == i, "Row"]))
			}

			data_corr_spat_na.subset <- data_corr_spat_na[data_corr_spat_na$TrtGeno %in% sel.geno,]
			data_corr_spat_na.subset$TrtGeno <- droplevels(data_corr_spat_na.subset$TrtGeno)

			#pdf("Phenoarch_platform_corrected_trait_app2_geno_pred_app1_selected_geno", height = 8, width = 12)
			aa <- xyplot(new_trait ~ Time|TrtGeno,  data = data_corr_spat_na.subset, groups = data_corr_spat_na.subset$pos,
								type = c("l"),
								xlab="Time", ylab = "Corrected trait", main = "Phenoarch platform data \n Corrected trait (App 2) + Genotypic predictions (App 1)", layout = c(4,2))
			bb <- xyplot(predicted.values ~ Time|TrtGeno,  data = Geno_pred.subset,
								xlab="", ylab = "", main = "", type = c("l"), col = "black", lwd = 2)
			aa + as.layer(bb)
			#dev.off()

