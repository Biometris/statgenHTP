# This file contains the code for the two two-step approaches for the analyses of the phenoarch platform
# PENDING: Population and individual trajectories based on Approach 2


# Clean the workspace
rm(list = ls())

# Load the needed packages
library(SpATS)
library(ggplot2)
# library(fields)
library(dplyr)
library(plyr)
library(lattice)
library(latticeExtra)
# library(Matrix)
library(data.table)

folder <- "C:/WUR/statgenHTP/"
setwd(folder)

dat.modif <- read.table(paste0(folder,"rawdata/Original_PAM_reshape.csv"),sep=',',h=T)
# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
dat.modif$ID <- paste0(dat.modif$x,"_",dat.modif$y)
dat.modif <- dat.modif[!is.na(dat.modif$pheno),]
# recover the time as time, better visualisation of time course
dat.modif$timepoint <- lubridate::ymd_hms(dat.modif$timepoints)
dat.modif <- dat.modif[order(dat.modif$ID,dat.modif$timepoints),]

# modification to be able to include the check in the model
# I forgot to mention that: quite often in the platform experiments,
# there is a highly replicated genotype (or 2, 3 ,4) that we can use in the model
# but we need to format the dataset
dat.modif$Geno <- as.character(dat.modif$Genotype)
dat.modif$Check <- NA
dat.modif$Check[dat.modif$Geno %in% c("col","ely","evo1","ler")] <- dat.modif$Geno[dat.modif$Geno %in% c("col","ely","evo1","ler")]
dat.modif$Check[!dat.modif$Geno %in% c("col","ely","evo1","ler")] <- "_pop"
dat.modif$Genobis <- dat.modif$Geno
dat.modif$Genobis[dat.modif$Check!="_pop"] <- NA


###################################################################################
######## Settings and output formatting

# In these analyses, the genotype is always included as random
# I have to check with Coté if we can also include it as fixed
	genotype.as.random = TRUE

# Trait to analyse
	trait = "pheno"

# Create factors
	dat.modif$Col = as.factor(dat.modif$x)
	dat.modif$Row = as.factor(dat.modif$y)
	dat.modif$Basin = as.factor(dat.modif$Basin)
	dat.modif$Image_pos = as.factor(dat.modif$Image_pos)
	dat.modif$Sowing_Block = as.factor(dat.modif$Sowing_Block)
	# Create num
	dat.modif$Colnum = as.numeric(dat.modif$x)
	dat.modif$Rownum = as.numeric(dat.modif$y)

# Transform the timepont to a numeric variable
	dat.modif$Time <- as.numeric(mapvalues(dat.modif$timepoints,
	                                       from = levels(factor(dat.modif$timepoints)),
	                                       to = 1:length(levels(factor(dat.modif$timepoints)))))

	# Time points (factor codification)
	times.fact <- levels(factor(dat.modif$timepoint))

	# Time points (numeric codification)
	times.num <- sort(unique(dat.modif$Time))

# Create an indicator for each plot (according to the row and column position)
	dat.modif$pos <- paste0("c",dat.modif$Col,"r",dat.modif$Row, sep = "")
	# I removed a plant that has very few measurements
	dat.modif <- dat.modif[dat.modif$pos!="c1r54",]
	dat.modif <- droplevels(dat.modif)

# Objects to save results for Approach 1
	Geno_pred <- Col_pred <- Row_pred <- NULL			# Predictions (may contain the intercept as well as other factors)
	Geno_BLUPs <- Col_BLUPs <- Row_BLUPs <- NULL		# BLUPs (only the estimated coefficients - they do not contain the intercept or other factors)

	ed_row <- ed_col <- ed_surface <- sigma2 <- sigma2_row <- sigma2_col  <- vector(length = length(times.fact))

	# When we are using the geno.decomp argument,
	# we have several levels of "populations" (combination of treatment and population),
	# so we'll have four different heritabilities = several columns.
	h2 <- sigma2_geno <- matrix(0, ncol = 1, nrow = length(times.fact))
	# colnames(h2) <- paste0("TrtPop",levels(dat.modif$TrtPop))

	# Object to save results for Approach 2
	data_corr_spat <- NULL

################################################################################
######## Loop on timepoint to run SpATS

# for (ti in 1:length(times.fact)) {
  for (ti in 1:2) {
	# ti = 31 ti = 50 ti = 68
	cat(times.fact[ti],'\n')

	### subset of dataset:
		dat.ti <- dat.modif[dat.modif$timepoints==times.fact[ti],]
		dat.ti <- droplevels(dat.ti)

	### number of segments for SpATS:
		nseg = c(nlevels(dat.ti$Col),
           	nlevels(dat.ti$Row))

	### Col and row (needed to obtain the BLUPs)
	Col_levels <- paste("Col", levels(dat.ti$Col), sep = "")
	Row_levels <- paste("Row", levels(dat.ti$Row), sep = "")

 # Fit the model using check
	# fit.SpATS2 <- SpATS::SpATS(response = trait,
	#                           fixed = ~ Sowing_Block + Image_pos + Check ,
	#                           random = ~ Col + Row ,
	#                           spatial = ~ SpATS::PSANOVA(Colnum, Rownum,
	#                                                      nseg = nseg,
	#                                                      nest.div=c(2,2)),
	#                           genotype = "Genobis",
	#                           genotype.as.random = TRUE,
	#                           data = dat.ti,
	#                           control = list(maxit = 50,
	#                                          tolerance = 1e-03,
	#                                          monitoring = 0))

	# Fit the model without the check effect
	fit.SpATS <- SpATS::SpATS(response = trait,
	                          fixed = ~  Sowing_Block + Image_pos,
	                          random = ~ Col + Row ,
	                          spatial = ~ SpATS::PSANOVA(Colnum, Rownum,
	                                                     nseg = nseg,
	                                                     nest.div=c(2,2)),
	                          genotype = "Geno",
	                          genotype.as.random = TRUE,
	                          # geno.decomp = "TrtPop",
	                          data = dat.ti,
	                          control = list(maxit = 50,
	                                         tolerance = 1e-03,
	                                         monitoring = 0))


	##########################################################################################
	# Approach 1: obtain the genotypic predictions
	##########################################################################################

		##############################
		# Genotype predictions
		##############################
		# Genotype prediction (including the effect of TrtPop, as well as the intercept)
		p_geno <- predict(fit.SpATS, which = c("Geno")) # ,"Check"

			# Include time point
		pred_g  <-  dplyr:::mutate(p_geno, Time = times.num[ti])

		# Select the needed variables for subsequent analyses
		pred_g  <- dplyr:::select(pred_g, Time, Geno, predicted.values, standard.errors)

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
		se_BLUPs_geno <- sqrt(diag(fit.SpATS$vcov$C11_inv[fit.SpATS$terms$geno$geno_names,
		                                                  fit.SpATS$terms$geno$geno_names]))

		# Create data.frame the needed variables for subsequent analyses
		df_geno <- data.frame(TrtGeno = names(BLUPs_geno),
		                      predicted.values = BLUPs_geno,
		                      standard.errors = se_BLUPs_geno,
		                      Time = times.num[ti])

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
  		df_Col <- data.frame(Col = names(BLUPs_Col),
  		                     predicted.values = BLUPs_Col,
  		                     standard.errors = se_BLUPs_Col,
  		                     Time = times.num[ti])
  		# Add results
  		Col_BLUPs <- rbind(Col_BLUPs, df_Col)


  		##############################
		# Row BLUPs
		##############################
		# BLUPs
		BLUPs_Row <- fit.SpATS$coeff[Row_levels]
		se_BLUPs_Row <- sqrt(diag(fit.SpATS$vcov$C22_inv[Row_levels,Row_levels]))
		df_Row <- data.frame(Row = names(BLUPs_Row),
		                     predicted.values = BLUPs_Row,
		                     standard.errors = se_BLUPs_Row,
		                     Time = times.num[ti])
		Row_BLUPs <- rbind(Row_BLUPs, df_Row)


		##############################
		# Heritabilities
		##############################
		h2.aux <- getHeritability(fit.SpATS)
		h2[ti,] <- h2.aux

		# if the geno.decomp is used:
		  # if(length(h2.aux) == 4) {
		  #  	h2[ti,] <- h2.aux
		  # } else {
		  #   aux <- rep(NA, 4)
		  #   aux[match(names(h2.aux), colnames(h2))] <- h2.aux
		  # 	h2[ti,] <- aux
		  # }

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
		p_a2 <- predict(fit.SpATS, which = c("Colnum","Rownum","Col","Row","Sowing_Block","Image_pos"))

	# Order the data and prediction according to the previous covariates
		data.ord <-  dat.ti[order(dat.ti$Colnum, dat.ti$Rownum),]
		p_a2 <- p_a2[order( p_a2$Colnum, p_a2$Rownum),] #p_a2$Genotype,

	# Add to the predictions some needed information
		p_a2$Genotype <- data.ord[,"Geno"]
		p_a2$pos <- data.ord[,"pos"]

	# Obtain the corrected trait
		p_a2$new_trait <- data.ord[,trait] - p_a2$predicted.values + fit.SpATS$coeff["Intercept"]

	# Include time point
		pred_a2 <- dplyr:::mutate(p_a2, Time = times.num[ti])

	# Select the needed variables for subsequent analyses
		pred_a2 <- dplyr:::select(pred_a2, new_trait, Genotype, Time, Sowing_Block, Image_pos,
		                          Colnum, Rownum, Col, Row, pos) #Treatment,

	# Add results
		data_corr_spat <- rbind(data_corr_spat, pred_a2)

}
  # Maybe not to keep in the function but pretty convenient: :D
	beepr::beep(2)

#############################################################
# Some results (all genotypes)
#############################################################
	# Raw data
		# NOTE: in this case xyplot connect the lines even when there are missing data. Thus, this information is lost.
		pdf("Phenovator_Rene_raw_data.pdf")
		xyplot(pheno ~ timepoint|Geno,
		       data = dat.modif, #[dat.modif$Genotype %in% c("ler","890","col","828","88","764"),],
		       groups = dat.modif$pos,
		       type = c("l"),
		       xlab="Time",
		       ylab = "Trait",
		       main = "Phenovator platform - Rene \n Raw data",
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
				dat.modif_na[dat.modif_na$pos == i, "Genotype"] <- na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Genotype"]))
				# dat.modif_na[dat.modif_na$pos == i, "Treatment"] <- na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Treatment"]))
				dat.modif_na[dat.modif_na$pos == i, "Colnum"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Colnum"]))
				dat.modif_na[dat.modif_na$pos == i, "Rownum"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Rownum"]))
				dat.modif_na[dat.modif_na$pos == i, "Col"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Col"]))
				dat.modif_na[dat.modif_na$pos == i, "Row"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Row"]))
			}

			pdf("Phenovator_Rene_raw_data_na.pdf")
			xyplot(pheno ~ timepoint|Genotype,
			       data = dat.modif_na,
			       groups = dat.modif_na$pos,
			       type = c("l"),
			       xlab="Time",
			       ylab = "Trait",
			       main = "Phenovator platform - Rene \n Raw data",
			       layout = c(5,5))
			dev.off()

managetime <- data.frame(timepoint = unique(dat.modif$timepoint),
                         Time = 1:73)
Geno_pred$timepoint <- managetime$timepoint[match(Geno_pred$Time,managetime$Time)]


	# Raw data + genotypic predictions (Approach 1)
		pdf("Phenovator_Rene_raw_data_geno_pred_app1_modRep.pdf")
		aa <- xyplot(pheno ~ timepoint|Geno,
		             data = dat.modif,
		             groups = dat.modif$pos,
		             type = c("l"),
		             xlab="Time", ylab = "Trait",
		             main = "Phenovator platform - Rene \n Raw data + Genotypic predictions (App 1)",
		             layout = c(5,5))


		bb <- xyplot(predicted.values ~ timepoint|Geno,
		             data = Geno_pred,
		             xlab="", ylab = "", main = "",
		             type = c("l"), col = "black", lwd = 2, layout = c(5,5))
		aa + as.layer(bb)
		dev.off()

	# Corrected trait (Approach 2) + + genotypic predictions (Approach 1)
		data_corr_spat$timepoint <- managetime$timepoint[match(data_corr_spat$Time,managetime$Time)]


		pdf("Phenovator_Rene_corrected_trait_app2_geno_pred_app1_modRep.pdf")
		aa <- xyplot(new_trait ~ timepoint|Genotype,
		             data = data_corr_spat,
		             groups = data_corr_spat$pos,
		             type = c("l"),
		             xlab="Time", ylab = "Corrected trait",
		             main = "Phenovator platform - Rene \n Corrected trait (App 2) + Genotypic predictions (App 1)",
		             layout = c(5,5))
		bb <- xyplot(predicted.values ~ timepoint|Geno,
		             data = Geno_pred,
		             xlab="", ylab = "", main = "",
							type = c("l"), col = "black", lwd = 2)
		aa + as.layer(bb)
		dev.off()

	# Heritabilities (Approach 1)
		pdf("Phenovator_Rene_heritability_modRep.pdf")
		# op <- par(mfrow = c(2,2))
		plot(h2[,1], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[1])
		lines(h2[,1])
		# plot(h2[,2], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[2])
		# lines(h2[,2])
		# plot(h2[,3], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[3])
		# lines(h2[,3])
		# plot(h2[,4], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[4])
		# lines(h2[,4])
		# par(op)
		dev.off()

	# Effective dimensions (Approach 1)
		pdf("Phenovator_Rene_ed_app1_modRep.pdf")
		par(mfrow = c(2,2))
		plot(times.num, ed_surface, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED spatial surface", cex.lab=1.5, type = "l")
		points(times.num, ed_surface)

		plot(times.num, ed_row, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED row RE", cex.lab=1.5, type = "l")
		points(times.num, ed_row)

		plot(times.num, ed_col, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED column RE", cex.lab=1.5, type = "l")
		points(times.num, ed_col)
		dev.off()

	# Variances (Approach 1)
		pdf("Phenovator_Rene_var_app1_modRep.pdf")
		par(mfrow = c(2,2))
	  	plot(times.num, sigma2, cex.axis=1.5, xlab="Time", ylab = expression(sigma^2), main = "Residual variance", cex.lab=1.5, type = "l")
			points(times.num, sigma2)

			plot(times.num, sigma2_col, cex.axis=1.5,  xlab="Time", ylab = expression(sigma^2), main = "Column RE variance", cex.lab=1.5, type = "l")
			points(times.num, sigma2_col)

			plot(times.num, sigma2_row, cex.axis=1.5, xlab="Time", ylab = expression(sigma^2), main = "Row RE variance", cex.lab=1.5, type = "l")
			points(times.num, sigma2_row)
		dev.off()

	# Row and Column (BLUPs and predictions)
		pdf("Phenovator_Rene_RowCol_pred_BLUPs_modRep.pdf",6,6)
		# Column
		xyplot(predicted.values ~ Time,  data = Col_pred, groups = Col_pred$Col,
					xlab="Time", ylab = "Column random factor prediction", main = "Phenovator Rene\n Column predictions", type = c("g", "p", "o"))

			xyplot(predicted.values ~ Time,  data = Col_BLUPs, groups = Col_BLUPs$Col,
					xlab="Time", ylab = "Column random factor BLUPs", main = "Phenovator Rene\n Column BLUPs", type = c("g", "p", "o"))

		# Row
			xyplot(predicted.values ~ Time,  data = Row_pred, groups = Row_pred$Row,
					xlab="Time", ylab = "Row random factor prediction", main = "Phenovator Rene\n Row predictions", type = c("g", "p", "o"))

			xyplot(predicted.values ~ Time,  data = Row_BLUPs, groups = Row_BLUPs$Row,
					xlab="Time", ylab = "Row random factor BLUPs", main = "Phenovator Rene\n Row BLUPs", type = c("g", "p", "o"))
		dev.off()

#############################################################
# Some results (some genotypes)
#############################################################
	sel.geno <- c("col","ely","evo1","ler","1293","1070","1724","1845")

	# Raw data, selected genotypes
		dat.subset <- dat.modif[dat.modif$Geno %in% sel.geno,]
		dat.subset <- droplevels(dat.subset)

	# Genotype predictions (Approach 1), selected genotypes
		Geno_pred.subset <- Geno_pred[Geno_pred$Geno %in% sel.geno,]
		Geno_pred.subset <- droplevels(Geno_pred.subset)

	# Corrected trait, selected genotypes
		data_corr_spat.subset <- data_corr_spat[data_corr_spat$Genotype %in% sel.geno,]
		data_corr_spat.subset <- droplevels(data_corr_spat.subset)

	# Raw data
		pdf("Phenovator_Rene_raw_and_corr_data_selected_geno_TrainingBis.pdf", height = 8, width = 12)
		xyplot(pheno ~ timepoint|Geno,
		       data = dat.subset,
		       groups = dat.subset$pos,
		       type = c("l"),
		       xlab="Time", ylab = "Trait",
		       main = "Phenovator platform - Rene \n Raw data",
		       layout = c(4,2))
		xyplot(new_trait ~ timepoint|Genotype,
		       data = data_corr_spat.subset,
		       groups = data_corr_spat.subset$pos,
		       type = c("l"),
		       xlab="Time", ylab = "Corrected trait",
		       main = "Phenovator platform - Rene \n Corrected trait (App 2)",
		       layout = c(4,2))
		dev.off()

	# Raw data + genotypic predictions (Approach 1)
		pdf("Phenovator_Rene_raw_data_geno_pred_app1_selected_geno_modRep.pdf", height = 8, width = 12)
		aa <- xyplot(pheno ~ timepoint|Geno,
		             data = dat.subset,
		             groups = dat.subset$pos,
		             type = c("l"),
		             xlab="Time", ylab = "Trait",
		             main = "Phenovator platform - Rene \n Raw data + Genotypic predictions (App 1)",
		             layout = c(4,2))


		bb <- xyplot(predicted.values ~ timepoint|Geno,
		             data = Geno_pred.subset,
		             xlab="", ylab = "", main = "",
		             type = c("l"), col = "black", lwd = 2, layout = c(4,2))
		aa + as.layer(bb)
		dev.off()

	# Corrected trait (Approach 2) + genotypic predictions (Approach 1)

		# NOTE: in this case xyplot connect the lines even when there are missing data. Thus, this information is lost.

		pdf("Phenovator_Rene_corrected_trait_app2_geno_pred_app1_selected_geno_modRep.pdf", height = 8, width = 12)
		aa <- xyplot(new_trait ~ timepoint|Genotype,
		             data = data_corr_spat.subset,
		             groups = data_corr_spat.subset$pos,
		             type = c("l"),
		             xlab="Time", ylab = "Corrected trait",
		             main = "Phenovator platform - Rene \n Corrected trait (App 2) + Genotypic predictions (App 1)",
		             layout = c(4,2))
		bb <- xyplot(predicted.values ~ timepoint|Genotype,
		             data = Geno_pred.subset,
		             xlab="", ylab = "", main = "",
		             type = c("l"), col = "black", lwd = 2)
		aa + as.layer(bb)
		dev.off()

write.table(data_corr_spat,paste0(folder,"outputOrigScript/Corrected_PAM_modRep.csv"),sep=',',row.names=F)
write.table(Geno_pred,paste0(folder,"outputOrigScript/BULPs_PAM_modRep.csv"),sep=',',row.names=F)
