# This file contains the code for the two two-step approaches for the analyses of the phenoarch platform
# PENDING: Population and individual trajectories based on Approach 2

#' @import data.table
#' @export
basefunction <- function(TD,
                         trait,
                         covariates = NULL,
                         geno.decomp = NULL,
                         useCheck = FALSE,
                         out1,
                         out2) {

  fitMods <- fitModels(TD = TD, trait = trait, covariates = covariates,
                       geno.decomp = geno.decomp, useCheck = useCheck)
  ## APPROACH 1
  pred_a1 <- lapply(X = fitMods, FUN = method1)
  # Add results
  genoPred <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "predGeno"))
  colPred <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "predCol"))
  rowPred <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "predRow"))
  genoBLUPs <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "BLUPsGeno"))
  colBLUPs <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "BLUPsCol"))
  rowBLUPs <- Reduce(f = rbind, x = lapply(X = pred_a1, FUN = `[[`, "BLUPsRow"))

  ##############################
  # Heritabilities
  ##############################
  h2 <- sapply(X = fitMods, FUN = SpATS::getHeritability)

  # # When we are using the geno.decomp argument,
  # # we have several levels of "populations" (combination of treatment and population),
  # # so we'll have four different heritabilities = several columns.
  # h2 <- sigma2_geno <- matrix(0, ncol = 1, nrow = length(times.fact))
  # # colnames(h2) <- paste0("TrtPop",levels(dat.modif$TrtPop))



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
  sigma2 <- sapply(X = fitMods, FUN = function(x) x$psi[1])
  edCol <- sapply(X = fitMods, FUN = function(x) x$eff.dim["colId"])
  edRow <- sapply(X = fitMods, FUN = function(x) x$eff.dim["rowId"])
  edSurface <- sapply(X = fitMods, FUN = function(x) {
    sum(x$eff.dim[c("f(colNum)", "f(rowNum)", "f(colNum):rowNum",
                            "colNum:f(rowNum)","f(colNum):f(rowNum)")])
    })
  sigma2Col <- sapply(X = fitMods, FUN = function(x) x$var.comp["colId"])
  sigma2Row <- sapply(X = fitMods, FUN = function(x) x$var.comp["rowId"])

  ## APPROACH2
  pred_a2 <- lapply(X = fitMods, FUN = method2)
  # Add results
  dataCorrSpat <- Reduce(f = rbind, x = pred_a2)


  # Maybe not to keep in the function but pretty convenient: :D
  beepr::beep(2)

  #
  #
  # #############################################################
  # # Some results (all genotypes)
  # #############################################################
  # # Raw data
  # # NOTE: in this case xyplot connect the lines even when there are missing data. Thus, this information is lost.
  # pdf("Phenovator_Rene_raw_data.pdf")
  # lattice::xyplot(pheno ~ timepoint|Geno,
  #                 data = dat.modif, #[dat.modif$Genotype %in% c("ler","890","col","828","88","764"),],
  #                 groups = dat.modif$pos,
  #                 type = c("l"),
  #                 xlab="Time",
  #                 ylab = "Trait",
  #                 main = "Phenovator platform - Rene \n Raw data",
  #                 layout = c(5,5))
  # dev.off()
  #
  # # NOTE: With the following code, we add explicitly the rows of the data sets that are missing, thus they will appear when plotting
  # # Very time consuming!!!
  # # There are missing data. Fill in the gaps
  # xy.coord <- data.table::data.table(expand.grid(Time = sort(unique(dat.modif$Time)), pos = unique(dat.modif$pos)))
  # data.table::setkeyv(xy.coord, c("Time", "pos"))
  # dat.modif <- data.table::data.table(dat.modif)
  # data.table::setkeyv(dat.modif, c("Time", "pos"))
  # dat.modif_na <- dat.modif[xy.coord]
  #
  # # Add the information we have for the missing values
  # for(i in unique(dat.modif_na$pos)) {
  #   dat.modif_na[dat.modif_na$pos == i, "Genotype"] <- na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Genotype"]))
  #   # dat.modif_na[dat.modif_na$pos == i, "Treatment"] <- na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Treatment"]))
  #   dat.modif_na[dat.modif_na$pos == i, "colNum"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "colNum"]))
  #   dat.modif_na[dat.modif_na$pos == i, "rowNum"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "rowNum"]))
  #   dat.modif_na[dat.modif_na$pos == i, "Col"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Col"]))
  #   dat.modif_na[dat.modif_na$pos == i, "Row"] <- 	na.omit(unique(dat.modif_na[dat.modif_na$pos == i, "Row"]))
  # }
  #
  # pdf("Phenovator_Rene_raw_data_na.pdf")
  # lattice::xyplot(pheno ~ timepoint|Genotype,
  #                 data = dat.modif_na,
  #                 groups = dat.modif_na$pos,
  #                 type = c("l"),
  #                 xlab="Time",
  #                 ylab = "Trait",
  #                 main = "Phenovator platform - Rene \n Raw data",
  #                 layout = c(5,5))
  # dev.off()
  #
  # managetime <- data.frame(timepoint = unique(dat.modif$timepoint),
  #                          Time = 1:73)
  # Geno_pred$timepoint <- managetime$timepoint[match(Geno_pred$Time,managetime$Time)]
  #
  #
  # # Raw data + genotypic predictions (Approach 1)
  # pdf("Phenovator_Rene_raw_data_geno_pred_app1_modRep.pdf")
  # aa <- lattice::xyplot(pheno ~ timepoint|Geno,
  #                       data = dat.modif,
  #                       groups = dat.modif$pos,
  #                       type = c("l"),
  #                       xlab="Time", ylab = "Trait",
  #                       main = "Phenovator platform - Rene \n Raw data + Genotypic predictions (App 1)",
  #                       layout = c(5,5))
  #
  #
  # bb <- lattice::xyplot(predicted.values ~ timepoint|Geno,
  #                       data = Geno_pred,
  #                       xlab="", ylab = "", main = "",
  #                       type = c("l"), col = "black", lwd = 2, layout = c(5,5))
  # aa + latticeExtra::as.layer(bb)
  # dev.off()
  #
  # # Corrected trait (Approach 2) + + genotypic predictions (Approach 1)
  # data_corr_spat$timepoint <- managetime$timepoint[match(data_corr_spat$Time,managetime$Time)]
  #
  #
  # pdf("Phenovator_Rene_corrected_trait_app2_geno_pred_app1_modRep.pdf")
  # aa <- lattice::xyplot(newTrait ~ timepoint|Genotype,
  #                       data = data_corr_spat,
  #                       groups = data_corr_spat$pos,
  #                       type = c("l"),
  #                       xlab="Time", ylab = "Corrected trait",
  #                       main = "Phenovator platform - Rene \n Corrected trait (App 2) + Genotypic predictions (App 1)",
  #                       layout = c(5,5))
  # bb <- lattice::xyplot(predicted.values ~ timepoint|Geno,
  #                       data = Geno_pred,
  #                       xlab="", ylab = "", main = "",
  #                       type = c("l"), col = "black", lwd = 2)
  # aa + latticeExtra::as.layer(bb)
  # dev.off()
  #
  # # Heritabilities (Approach 1)
  # pdf("Phenovator_Rene_heritability_modRep.pdf")
  # # op <- par(mfrow = c(2,2))
  # plot(h2[,1], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[1])
  # lines(h2[,1])
  # # plot(h2[,2], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[2])
  # # lines(h2[,2])
  # # plot(h2[,3], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[3])
  # # lines(h2[,3])
  # # plot(h2[,4], ylim = c(0.5,1), xlab = "Time", ylab = "h2", main =  colnames(h2)[4])
  # # lines(h2[,4])
  # # par(op)
  # dev.off()
  #
  # # Effective dimensions (Approach 1)
  # pdf("Phenovator_Rene_ed_app1_modRep.pdf")
  # par(mfrow = c(2,2))
  # plot(times.num, ed_surface, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED spatial surface", cex.lab=1.5, type = "l")
  # points(times.num, ed_surface)
  #
  # plot(times.num, ed_row, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED row RE", cex.lab=1.5, type = "l")
  # points(times.num, ed_row)
  #
  # plot(times.num, ed_col, cex.axis=1.5, xlab="Time", ylab = "ED", main = "ED column RE", cex.lab=1.5, type = "l")
  # points(times.num, ed_col)
  # dev.off()
  #
  # # Variances (Approach 1)
  # pdf("Phenovator_Rene_var_app1_modRep.pdf")
  # par(mfrow = c(2,2))
  # plot(times.num, sigma2, cex.axis=1.5, xlab="Time", ylab = expression(sigma^2), main = "Residual variance", cex.lab=1.5, type = "l")
  # points(times.num, sigma2)
  #
  # plot(times.num, sigma2_col, cex.axis=1.5,  xlab="Time", ylab = expression(sigma^2), main = "Column RE variance", cex.lab=1.5, type = "l")
  # points(times.num, sigma2_col)
  #
  # plot(times.num, sigma2_row, cex.axis=1.5, xlab="Time", ylab = expression(sigma^2), main = "Row RE variance", cex.lab=1.5, type = "l")
  # points(times.num, sigma2_row)
  # dev.off()
  #
  # # Row and Column (BLUPs and predictions)
  # pdf("Phenovator_Rene_RowCol_pred_BLUPs_modRep.pdf",6,6)
  # # Column
  # lattice::xyplot(predicted.values ~ Time,  data = Col_pred, groups = Col_pred$Col,
  #                 xlab="Time", ylab = "Column random factor prediction", main = "Phenovator Rene\n Column predictions", type = c("g", "p", "o"))
  #
  # lattice::xyplot(predicted.values ~ Time,  data = Col_BLUPs, groups = Col_BLUPs$Col,
  #                 xlab="Time", ylab = "Column random factor BLUPs", main = "Phenovator Rene\n Column BLUPs", type = c("g", "p", "o"))
  #
  # # Row
  # lattice::xyplot(predicted.values ~ Time,  data = Row_pred, groups = Row_pred$Row,
  #                 xlab="Time", ylab = "Row random factor prediction", main = "Phenovator Rene\n Row predictions", type = c("g", "p", "o"))
  #
  # lattice::xyplot(predicted.values ~ Time,  data = Row_BLUPs, groups = Row_BLUPs$Row,
  #                 xlab="Time", ylab = "Row random factor BLUPs", main = "Phenovator Rene\n Row BLUPs", type = c("g", "p", "o"))
  # dev.off()
  #
  # #############################################################
  # # Some results (some genotypes)
  # #############################################################
  # sel.geno <- c("col","ely","evo1","ler","1293","1070","1724","1845")
  #
  # # Raw data, selected genotypes
  # dat.subset <- dat.modif[dat.modif$Geno %in% sel.geno,]
  # dat.subset <- droplevels(dat.subset)
  #
  # # Genotype predictions (Approach 1), selected genotypes
  # Geno_pred.subset <- Geno_pred[Geno_pred$Geno %in% sel.geno,]
  # Geno_pred.subset <- droplevels(Geno_pred.subset)
  #
  # # Corrected trait, selected genotypes
  # data_corr_spat.subset <- data_corr_spat[data_corr_spat$Genotype %in% sel.geno,]
  # data_corr_spat.subset <- droplevels(data_corr_spat.subset)
  #
  # # Raw data
  # pdf("Phenovator_Rene_raw_and_corr_data_selected_geno_TrainingBis.pdf", height = 8, width = 12)
  # lattice::xyplot(pheno ~ timepoint|Geno,
  #                 data = dat.subset,
  #                 groups = dat.subset$pos,
  #                 type = c("l"),
  #                 xlab="Time", ylab = "Trait",
  #                 main = "Phenovator platform - Rene \n Raw data",
  #                 layout = c(4,2))
  # lattice::xyplot(newTrait ~ timepoint|Genotype,
  #                 data = data_corr_spat.subset,
  #                 groups = data_corr_spat.subset$pos,
  #                 type = c("l"),
  #                 xlab="Time", ylab = "Corrected trait",
  #                 main = "Phenovator platform - Rene \n Corrected trait (App 2)",
  #                 layout = c(4,2))
  # dev.off()
  #
  # # Raw data + genotypic predictions (Approach 1)
  # pdf("Phenovator_Rene_raw_data_geno_pred_app1_selected_geno_modRep.pdf", height = 8, width = 12)
  # aa <- lattice::xyplot(pheno ~ timepoint|Geno,
  #                       data = dat.subset,
  #                       groups = dat.subset$pos,
  #                       type = c("l"),
  #                       xlab="Time", ylab = "Trait",
  #                       main = "Phenovator platform - Rene \n Raw data + Genotypic predictions (App 1)",
  #                       layout = c(4,2))
  #
  #
  # bb <- lattice::xyplot(predicted.values ~ timepoint|Geno,
  #                       data = Geno_pred.subset,
  #                       xlab="", ylab = "", main = "",
  #                       type = c("l"), col = "black", lwd = 2, layout = c(4,2))
  # aa + latticeExtra::as.layer(bb)
  # dev.off()
  #
  # # Corrected trait (Approach 2) + genotypic predictions (Approach 1)
  #
  # # NOTE: in this case xyplot connect the lines even when there are missing data. Thus, this information is lost.
  #
  # pdf("Phenovator_Rene_corrected_trait_app2_geno_pred_app1_selected_geno_modRep.pdf", height = 8, width = 12)
  # aa <- lattice::xyplot(newTrait ~ timepoint|Genotype,
  #                       data = data_corr_spat.subset,
  #                       groups = data_corr_spat.subset$pos,
  #                       type = c("l"),
  #                       xlab="Time", ylab = "Corrected trait",
  #                       main = "Phenovator platform - Rene \n Corrected trait (App 2) + Genotypic predictions (App 1)",
  #                       layout = c(4,2))
  # bb <- lattice::xyplot(predicted.values ~ timepoint|Geno,
  #                       data = Geno_pred.subset,
  #                       xlab="", ylab = "", main = "",
  #                       type = c("l"), col = "black", lwd = 2)
  # aa + latticeExtra::as.layer(bb)
  # dev.off()

  write.csv(genoPred, file = out1, row.names = FALSE)
  write.csv(dataCorrSpat, file = out2, row.names = FALSE)

  return(list(genoPred = genoPred, dataCorrSpat = dataCorrSpat))
}
