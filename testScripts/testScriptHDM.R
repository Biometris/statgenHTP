library(statgenHTP)
Pheno.cor.ran.sub3 <- readRDS("../htp_two_stage_approach/PhenoArchHDM")

testDat <- Pheno.cor.ran.sub3[Pheno.cor.ran.sub3$genotype %in%
                                c("WD_GenoA44","WD_GenoB20","WW_GenoA44","WW_GenoB20"),]

fit.psHDM  <- fitSplineHDM(response = "LeafArea_corr",
                           time = "timeNumber",
                           pop = "geno.decomp",
                           geno = "genotype",
                           plant = "plotId",
                           weights = "wt",
<<<<<<< HEAD
                           #data = testDat,
                           data = Pheno.cor.ran.sub3,
=======
                           data = testDat,
>>>>>>> 0b82344 (Added creation of full grid to fitSplineHDM)
                           dif.var = list(geno = FALSE, plant = FALSE),
                           smooth.pop = list(nseg = 10, bdeg = 3, pord = 2),
                           smooth.geno = list(nseg = 10, bdeg = 3, pord = 2),
                           smooth.plant = list(nseg = 10, bdeg = 3, pord = 2),
                           offset = NULL, family = gaussian(), maxit = 200,
                           trace = TRUE, thr = 1e-03)

<<<<<<< HEAD


=======
>>>>>>> 0b82344 (Added creation of full grid to fitSplineHDM)
pred.psHDM <- predict(object = fit.psHDM,
                      newtimes = fit.psHDM$time,
                      pred = list(pop = TRUE, geno = TRUE, plant = TRUE),
                      se = list(pop = TRUE, geno = TRUE, plant = FALSE))

plot(object = pred.psHDM,
     geno.sub = c("WD_GenoA44","WD_GenoB20","WW_GenoA44","WW_GenoB20"),
     # Order of geno.sub.names must match with the orden in geno.sub
     geno.sub.names = c("Geno 44 - Panel 1 - WD","Geno 20 - Panel 2 - WD","Geno 44 - Panel 1 - WW","Geno 20 - Panel 2 - WW"),
     # Specify the order you wish in the plot
     geno.sub.order = c(1,3,2,4),
     ylab = expression(paste("Spatially corrected leaf area (",m^2, " plan",t^{-1},")")),
     xlab = "Days since January 1st",
     my.theme = my.theme(my.size = 15),
     ask = TRUE,
     global.main = list(pop.tra = "(a) Population-specific growth curves",
                        pop.tra.deriv1 = "First-order derivative of the population-specific growth curves",
                        geno.tra = "Genotype-specific growth curves",
                        geno.tra.deriv1 = "First-order derivative of the genotype-specific growth curves",
                        geno.dev = "(b) Genotype-specific deviations (all genotypes)",
                        plant.tra = "(d) Plant- and genotype-specific growth curves"))
