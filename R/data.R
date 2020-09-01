#' Growth chamber data for an Arabidopsis experiment in the Phenovator platform.
#'
#' A dataset containing growth chamber data from an experiment with Arabidopsis
#' in the Phenovator platform (WUR, Netherlands, Flood et al. 2016). It consists of
#' one experiment with 1440 plants grown in a growth chamber. The number of tested
#' genotypes is 192 with 6-7 replicates per genotype. Four reference genotypes were
#' also tested with 15 or 30 replicates. The studied trait is the photosystem II
#' efficiency (EffpsII) extracted from the pictures over time (van Rooijen et al. 2017).
#'
#' @format A data.frame with 103,839 rows and 10 columns:
#' \describe{
#'   \item{Genotype}{Genotypes}
#'   \item{Basin}{Table of experiment}
#'   \item{Replicate}{Block define after sowing for post-blocking. They are not full-resolvable blocks.}
#'   \item{Image_pos}{Position of the camera}
#'   \item{x}{Row coordinate}
#'   \item{y}{Column coordinate}
#'   \item{Sowing_Position}{Unique pot ID}
#'   \item{timepoints}{time of picture}
#'   \item{EffpsII}{Efficiency of the photosystem II}
#'   \item{pos}{Unique pot ID using rowcol coordinates}
#' }
#'
#'
#' @references Rooijen, Roxanne van, Willem Kruijer, René Boesten,
#' Fred A. van Eeuwijk, Jeremy Harbinson, and Mark G. M. Aarts. 2017.
#' “Natural Variation of YELLOW SEEDLING1 Affects Photosynthetic Acclimation
#' of Arabidopsis Thaliana.” Nature Communications 8 (1).
#' doi:10.1038/s41467-017-01576-3.
#'
#' Flood, Pádraic J., Willem Kruijer, Sabine K. Schnabel, Rob van der Schoor,
#' Henk Jalink, Jan F. H. Snel, Jeremy Harbinson, and Mark G. M. Aarts. 2016.
#' “Phenomics for Photosynthesis, Growth and Reflectance in Arabidopsis Thaliana
#' Reveals Circadian and Long-Term Fluctuations in Heritability.”
#' Plant Methods 12 (1): 14.
#' doi:10.1186/s13007-016-0113-y.
"PhenovatorDat1"


#' Arabidopsis data corrected for spatial trends.
#'
#' This dataset contains the corrected data obtained by (1) removing point outliers
#' and (2) running a spatial model on the PhenovatorDat1 dataset.
#' See the vignettes for details.
#'
#' @format A data.frame with 103,801 rows and 11 columns:
#' \describe{
#'   \item{timeNumber}{Time number obtained after formatting the original dataset with the function createTP.}
#'   \item{timePoint}{Original time point.}
#'   \item{EffpsII_corr}{Efficiency of the photosystem II, corrected data}
#'   \item{EffpsII}{Efficiency of the photosystem II, raw data}
#'   \item{genotype}{Genotypes}
#'   \item{repId}{Block define after sowing for post-blocking.}
#'   \item{Image_pos}{Position of the camera}
#'   \item{check}{Status of the genotypes: check for the reference genotypes, noCheck for the others.}
#'   \item{colId}{Column coordinate}
#'   \item{rowId}{Row coordinate}
#'   \item{plotId}{Unique pot ID using rowcol coordinates}
#' }
"spatCorrVator"

#' Greenhouse data for a maize experiment in the Phenoarch platform.
#'
#' A dataset containing greenhouse data from an experiment with maize in the
#' Phenoarch platform (INRAE, France, Cabrera-Bosquet et al. 2016). It consists of
#' one experiment with 1671 plants grown in a greenhouse under two water scenarios,
#' well-watered (WW) and water deficit (WD). There are two populations of genotypes,
#' panel 1 and panel 2. Panel 1 contains 60 genotypes with 14 replicates: 7 in WW
#' and 7 in WD. Note that there are more plants per replicates than one (about 24
#' plants per genotypes). Panel 2 contains 30 genotypes with 8 replicates, 4 in WW
#' and 4 in WD. The studied trait is the leaf area extracted from the pictures over
#' time. Plants were pictured every day for 35 days.
#'
#' @format A data.frame with 42,536 rows and 14 columns:
#' \describe{
#'   \item{Date}{Date of measurement}
#'   \item{pos}{Unique pot using rowcol coordinate}
#'   \item{geno}{Genotype}
#'   \item{Scenario}{Water regime, WW or WD}
#'   \item{population}{Panel 1 or 2}
#'   \item{Row}{Pot position on the conveyor belt (i.e. row coordinate)}
#'   \item{Col}{Line of conveyor belt (i.e. column coordinate)}
#'   \item{Biomass_Estimated}{Biomass from the picture}
#'   \item{LA_Estimated}{Leaf area from the picture}
#'   \item{Height_Estimated}{Plant height from the picture}
#'   \item{count_leaf}{Number of leaves manually scored}
#'   \item{phyllocron}{Leaf emission rate}
#'   \item{TT}{Thermal time cumulated from emergence}
#' }
#'
#'
#' @references Cabrera-Bosquet, Llorenç, Fournier Christian, Brichet Nicolas,
#' Welcker Claude, Suard Benoît, and Tardieu François. 2016.
#' “High-throughput estimation of incident light, light interception and
#' radiation-use efficiency of thousands of plants in a phenotyping platform.”
#' New Phytologist 212 (1): 269-81.
#' doi:10.1111/nph.14027.
"PhenoarchDat1"

#' Maize data corrected for spatial trends.
#'
#' This dataset contains the corrected data obtained by (1) removing point outliers
#' and (2) running a spatial model on the PhenoarchDat1 dataset.
#' See the vignettes for details.
#'
#' @format A data.frame with 37,038 rows and 9 columns:
#' \describe{
#'   \item{timeNumber}{Time number obtained after formatting the original dataset with the function createTP.}
#'   \item{timePoint}{Original time point.}
#'   \item{Biomass_Estimated_corr}{Biomass, corrected data}
#'   \item{Biomass_Estimated}{Biomass from the picture, raw data}
#'   \item{genotype}{Genotypes}
#'   \item{geno.decomp}{Combination of treatment levels to decompose the genotypic variance (see vignettes)}
#'   \item{colId}{Column coordinate}
#'   \item{rowId}{Row coordinate}
#'   \item{plotId}{Unique pot ID using rowcol coordinates}
#' }
"spatCorrArch"

#' Greenhouse data for an experiment in the RootPhAir platform.
#'
#' A dataset containing greenhouse data from two experiments in the RootPhAir
#' platform (UCLouvain, France, ref). It consists of two experiments, each in two
#' aeroponic tanks with a capacity of 495 plants each tank, located in the same
#' greenhouse.
#' with 1671 plants grown in a greenhouse in two tanks
#' The studied traits are the root tip coordinates in y and x axis, extracted from the
#' pictures over time. Plants were pictured every 2 hours for 10 days.
#'
#' @format A data.frame with 16,275 rows and 10 columns:
#' \describe{
#'   \item{Exp}{Experiment number}
#'   \item{thermalTime}{Thermal time cumulated}
#'   \item{Genotype}{Genotype}
#'   \item{plantId}{Unique pot using tank and rowcol coordinate}
#'   \item{Tank}{Tank A or B}
#'   \item{Strip}{Number of strip of five plants (i.e. row coordinate)}
#'   \item{Pos}{Position within th strip (i.e. column coordinate)}
#'   \item{tipPos_x}{Position of the root tip in x axis}
#'   \item{tipPos_y}{Position of the root tip in y axis}
#'   \item{Time}{Time of measurement}
#' }
#'
#'
#' @references Author 1 author 2. date.
#' “A Title.”
#' Journal.
#' doi:xxxx.
"RootDat1"

#' Root data corrected for point outliers.
#'
#' This dataset contains the corrected data obtained by removing point outliers
#' on the RootDat1 dataset. See the vignettes for details.
#'
#' @format A data.frame with 15,934 rows and 8 columns:
#' \describe{
#'   \item{timePoint}{Original time points, date and time}
#'   \item{Date}{Date}
#'   \item{thermalTime}{Thermal time cumulated}
#'   \item{Exp}{Experiment number}
#'   \item{genotype}{Genotypes}
#'   \item{Tank}{Tank in the greenhouse}
#'   \item{plotId}{Unique pot ID using rowcol coordinates}
#'   \item{rowId}{Row coordinate}
#'   \item{colId}{Column coordinate}
#'   \item{tipPos_y}{Position of the root tip in y axis}
#' }
"noCorrRoot"

