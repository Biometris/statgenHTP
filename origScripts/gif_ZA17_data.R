

### save spatial trends as gif
### requires the pakcage «magick» : https://github.com/ropensci/magick
library(animation)

### also calls a function by Martin
source("~/Documents/PostDoc/EPPN2020/Platform/M3P/Phenoarch/ZA17/PSANOVA.plots3.R")

saveGIF(for (ti in 1:length(times.fact)) 
# saveGIF(for (ti in 1:5) 
  {
    # ti = 18

    ### subset of dataset:
    dat.ti <- dat.modif[dat.modif$Date==times.fact[ti],]
    dat.ti <- droplevels(dat.ti)
    
    day = times.fact[ti]
    time = ti
    title = paste0("Trait: ",trait,", day: ",day,", time: ",time)
    cat(ti,'\n')
  
  ### number of segments for SpATS:
  nseg = c(nlevels(dat.ti$Col),
           nlevels(dat.ti$Row))
  
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
  
  
  PSANOVA.plots(fit.SpATS, 
                which='global',
                cex.lab = 1.5, 
                cex.axis = 1.2, 
                main = title)
  
  
} , movie.name = 'spatialtrends.gif')
