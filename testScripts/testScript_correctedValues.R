rm(list = ls())
# Load the needed packages

# install.packages(c("fields","dplyr","plyr","lattice","latticeExtra","Matrix","data.table","car"))
library(SpATS)
library(fields)
library(dplyr)
library(plyr)
library(lattice)
library(latticeExtra)
library(Matrix)
library(data.table)


######### Raw data
dat.modif <- read.table("C:/Users/rossu027/Downloads/Data_Phenoarchplatform_Cote.csv",sep=',',h=T, stringsAsFactors = TRUE)
dat.modif <- na.omit(dat.modif)

# Create factors
dat.modif$Col = as.factor(dat.modif$Col)
dat.modif$Row = as.factor(dat.modif$Row)

# Transform the timepont to a numeric variable
dat.modif$Time <- as.numeric(mapvalues(dat.modif$timepoint, from = levels(factor(dat.modif$timepoint)), to = 1:length(levels(factor(dat.modif$timepoint)))))

# Time points (factor codification)
times.fact <- levels(factor(dat.modif$Time))

# Time points (numeric codification)
times.num <- sort(unique(dat.modif$Time))

# Create an indicator for each plot (according to the row and column position)
dat.modif$pos <- paste0("c",dat.modif$Col,"r",dat.modif$Row, sep = "")

# Order the data
dat.modif <- dat.modif[order(dat.modif$TrtPop, dat.modif$TrtGeno, dat.modif$pos, dat.modif$Time),]

######## Initialize the first-stage approach
# In these analyses, the genotype is always included as random
genotype.as.random = TRUE

# Trait to analyse
trait = "Trait"

# Object to save results for Approach 2
data.fit   <- NULL # Predictions
Geno_BLUPs <- NULL # BLUPs (only the estimated coefficients - they do not contain the intercept or other factors)
ed_row <- ed_col <- ed_surface <- sigma2 <- sigma2_row <- sigma2_col  <- vector(length = length(times.fact))
h2 <- sigma2_geno <- matrix(0, ncol = 4, nrow = length(times.fact)) # Since we are using the geno.decomp argument, and we have 4 "populations" (combination of treatment and population), we have four different heritabilities.
# colnames(h2) <- paste0("TrtPop",levels(dat.modif_na$TrtPop))
colnames(h2) <- paste0("TrtPop",levels(dat.modif$TrtPop))
model <- list()

################################################################################
######## Loop on timepoint to run SpATS
start <- proc.time()[3]

for (ti in 1:length(times.fact)) {
  # for (ti in 1:1) {
  cat(times.fact[ti],'\n')

  ### subset of dataset:
  dat.ti <- dat.modif[dat.modif$Time==times.fact[ti],]
  # dat.ti <- dat.modif_na[dat.modif_na$Time==times.fact[ti],]
  dat.ti <- droplevels(dat.ti)

  ### number of segments for SpATS:
  nseg = c(nlevels(dat.ti$Col) / 2,
           nlevels(dat.ti$Row) / 2)

  # Fit the model
  # With Geno as random
  fit.SpATS <- SpATS::SpATS(response = trait,
                            fixed = ~ TrtPop ,
                            random = ~ Col + Row ,
                            spatial = ~ SpATS::PSANOVA(Colnum, Rownum, nseg = nseg, nest.div=c(2,2)),
                            genotype = "TrtGeno",
                            genotype.as.random = TRUE,
                            geno.decomp = "TrtPop",
                            data = dat.ti,
                            control = list(maxit = 50, tolerance = 1e-03, monitoring = 1))
  model[[ti]] <- fit.SpATS

  # # Include in the prediction the factors (variables) whose effect we are interested in removing
  # p_a2 <- predict(fit.SpATS, which = c("Colnum","Rownum","Col","Row"))
  # Obtain BLUPs: with return.vcov.matrix = TRUE we obtain the variace-covariance matrix for the predictions
  # Include in the prediction the factors (variables) whose effects we are interested in
  which.cov <- c("TrtPop","TrtGeno")
  p         <- predict(fit.SpATS, which = which.cov, return.vcov.matrix = TRUE, predFixed = "marginal")

  # Data used in the fit plus residuals
  data.fit.ti     <- fit.SpATS$data
  data.fit.ti$res <- fit.SpATS$residuals

  # Match observations with predictions, based on covariates in "which.cov"
  data.fit.ti <- left_join(data.fit.ti, select(p, TrtPop, TrtGeno, predicted.values, standard.errors), by = which.cov)

  # New response
  # Obtain the corrected trait
  # data.fit.ti$new_trait <- dat.ti[,trait] - data.fit.ti$predicted.values + fit.SpATS$coeff["Intercept"]
  data.fit.ti <- mutate(data.fit.ti, yield_corr = predicted.values + res) # BLUPs + residuals

  ###############################################
  # Weights based on standard errors
  ###############################################
  # We add the residual error fit.SpATS$psi[1] to the variance of the predictions (since we add the residuals)
  data.fit.ti <- mutate(data.fit.ti, weights_1 = 1/(sqrt(data.fit.ti$standard.errors^2 + fit.SpATS$psi[1])))

  ###################################################################################################
  # Weights based on the inverse of the var-cov (vcov) matrix
  # NOTE: not sure whether this is completely correct.
  # We obtain the vcov + residual variance for the predictions but not for the observations
  # The inverse of both matrices (predictions vs observations) is not the same
  ###################################################################################################
  # Get the vcov matrix from p
  vcov        <- attr(p, "vcov")

  # We add the residual error m0$psi[1] to the diagonal of the vcov matrix (since we add the residuals)
  vcov_comb   <- as.matrix(vcov) + fit.SpATS$psi[1]*diag(nrow(vcov))
  p$weights_2 <- sqrt(diag(solve(vcov_comb)))

  # Match observations with (new/corrected) predictions, based on covariate in which.cov
  data.fit.ti <- left_join(data.fit.ti, select(p, TrtPop, TrtGeno, weights_2), by = which.cov)

  data.fit    <- rbind(data.fit, data.fit.ti)
  # # Plot the weights
  # plot(data.fit.ti$weights_2, data.fit.ti$weights_1)
  # abline(lm(data.fit.ti$weights_1~ data.fit.ti$weights_2))

  ##############################
  # Genotype BLUPs
  ##############################
  # BLUPs
  BLUPs_geno <- fit.SpATS$coeff[fit.SpATS$terms$geno$geno_names[fit.SpATS$terms$geno$ndx]]

  # Standard error
  se_BLUPs_geno <- sqrt(diag(fit.SpATS$vcov$C11_inv[fit.SpATS$terms$geno$geno_names[fit.SpATS$terms$geno$ndx],fit.SpATS$terms$geno$geno_names[fit.SpATS$terms$geno$ndx]]))

  # Create data.frame the needed variables for subsequent analyses
  df_geno <- data.frame(TrtPop = p$TrtPop, TrtGeno = names(BLUPs_geno), predicted.values = BLUPs_geno, standard.errors = se_BLUPs_geno, Time = times.num[ti])

  # Add results
  Geno_BLUPs = rbind(Geno_BLUPs, df_geno)

  ##############################
  # Heritabilities
  ##############################
  h2.aux <- getHeritability(fit.SpATS)

  if(length(h2.aux) == 4) {
    h2[ti,] <- h2.aux
  } else {
    aux <- rep(NA, 4)
    aux[match(names(h2.aux), colnames(h2))] <- h2.aux
    h2[ti,] <- aux
  }

  # Object to be returned (for the momment only few information)
  res            <- list()
  res$data.fit   <- data.fit
  res$Geno_BLUPs <- Geno_BLUPs
  res$h2         <- h2
  res$model      <- model
}

end <- proc.time()[3]
spats.time <- end - start
spats.time
# Save results
save("res", file = "Data/Data_PhenoArch_First_Stage_GenoRan_SpATS_Weights_Paper.RData")

# Number of variance components: (error variances are missing)
# sum(unlist(lapply(res$model, function(x) length(x$var.comp))))
# Number of coefficients
# sum(unlist(lapply(res$model, function(x) length(x$coeff))))

# ##########################################################
# # Comparison of the results with statgenHTD (v0.003) and SpATS (v1.0-11) with Weights and marginal prediction
# ##########################################################
#
# # Load the results for both methods
# load("Data/Data_PhenoArch_First_Stage_GenoRan_SpATS_Weights_August_2020_onetime.RData")
# load("Data/Data_PhenoArch_First_Stage_GenoRan_statgenHTP_August_2020.RData")
#
# # Comparing the corrected phenotype
# statgen <- pheno.cor.ran[order(pheno.cor.ran$geno.decomp, pheno.cor.ran$genotype, pheno.cor.ran$plotId, pheno.cor.ran$timeNumber),]
# spats   <- res$data.fit[order(res$data.fit$TrtPop, res$data.fit$TrtGeno, res$data.fit$pos, res$data.fit$Time),]
# cbind(statgen[,c("timeNumber","Trait","Trait_corr")], spats[,c("Time","Trait","yield_corr")])
# plot(statgen[,"Trait_corr"], spats[,"yield_corr"], type = "l", xlab = "statgenHTP", ylab = "SpATS")
# range(statgen[,"Trait_corr"] - spats[,"yield_corr"])
# lines(statgen[,"Trait"], spats[,"Trait"], col = 2)
# range(statgen[,"Trait"] - spats[,"Trait"])
#
# # Comparing the heritabilities
# statgen.h2 <- Pheno.h2
# spats.h2 <- res$h2[1,]
#
# # Comparing the variance components
# # For SpATS
# res$model$var.comp
# # For statgenHTP
# Pheno.vc.ran
#
# # Comparing the effective dimensions
# # For SpATS
# res$model$eff.dim
# # For statgenHTP
# Pheno.ed.ran

library(statgenHTP)

phenoTParch <- createTimePoints(dat = dat.modif,
                                experimentName = "Phenoarch",
                                genotype = "Genotype",
                                timePoint = "timepoint",
                                plotId = "pos",
                                rowNum = "Rownum",
                                colNum = "Colnum")

modPhenoSpGD <- fitModels(TP = phenoTParch,
                          trait = "Trait",
                          geno.decomp = c("Treatment", "Population"))

phenoTPcorr <- getCorrected(modPhenoSpGD)


modPhenoSpAs <- fitModels(TP = phenoTParch,
                          trait = "Trait",
                          geno.decomp = c("Treatment", "Population"),
                          engine = "asreml", spatial = TRUE)

phenoTPcorrAs <- getCorrected(modPhenoSpAs)


comp <- merge(phenoTPcorr, phenoTPcorrAs,
              by = c("timeNumber", "timePoint", "genotype", "geno.decomp", "rowId", "colId", "plotId"))

comp1 <- comp[comp$timeNumber == 1,]
plot(comp1$Trait_corr.x, comp1$Trait_corr.y)
abline(0,1,col = "blue")

plot(comp1$wt.x, comp1$wt.y)
abline(0,1,col = "blue")

comp35 <- comp[comp$timeNumber == 35,]
plot(comp35$Trait_corr.x, comp35$Trait_corr.y)
abline(0,1,col = "blue")

plot(comp35$wt.x, comp35$wt.y)
abline(0,1,col = "blue")

