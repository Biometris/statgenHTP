#' Growth chamber data for an Arabidopsis experiment in the Phenovator platform.
#'
#' A dataset containing growth chamber data from an experiment with Arabidopsis
#' in the Phenovator platform (WUR, Netherlands, Flood et al. 2016). It consists of
#' one experiment with 1440 plants grown in a growth chamber. The number of tested
#' genotypes is 192 with 6-7 replicates per genotype. Four reference genotypes were
#' also tested with 15 or 30 replicates. The studied trait is the photosystem II
#' efficiency (EffpsII) extracted from the pictures over time (van Rooijen et al. 2017).
#'
#' @format A data.frame with 103839 rows and 11 columns:
#' \describe{
#'   \item{Genotype}{Genotypes}
#'   \item{Basin}{Table of experiment}
#'   \item{Replicate}{Replicate}
#'   \item{Image_pos}{Position of the camera}
#'   \item{Sowing_Block}{Block defined for post-blocking}
#'   \item{x}{row coordinate}
#'   \item{y}{column coordinate}
#'   \item{Sowing_Position}{unique pot ID}
#'   \item{timepoints}{time of picture}
#'   \item{EffpsII}{trial, a combination of location and year}
#'   \item{pos}{unique pot using rowcol coordinate}
#' }
#'
#' @source \url{}
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

#' Greenhouse data for a maize experiment in the Phenoarch platform.
#'
#' A dataset containing greenhouse data from an experiment with maize in the
#' Phenoarch platform (INRA, France, Cabrera-Bosquet et al. 2016). It consists of
#' one experiment with 1671 plants grown in a greenhouse under two water scenarios,
#' well-watered (WW) and water deficit (WD). There are two populations of genotypes,
#' panel 1 and panel 2. Panel 1 contains 60 genotypes with 14 replicates: 7 in WW
#' and 7 in WD. Note that thereare more plants per replicates than one (about 24
#' plants per genotypes). Panel 2 contains 30 genotypes with 8 replicates, 4 in WW
#' and 4 in WD. The studied trait is the leaf area extracted from the pictures over
#' time. Plants were pictured every day for 35 days.
#'
#' @format A data.frame with 50405 rows and 10 columns:
#' \describe{
#'   \item{geno}{Genotype}
#'   \item{population}{Panel 1 or 2}
#'   \item{Scenario}{Water regime}
#'   \item{Rep}{Replicate}
#'   \item{Line}{line of conveyor belt (i.e column coordinate)}
#'   \item{Position}{pot position on the conveyor belt (i.e. row coordinate)}
#'   \item{Pot}{unique pot ID}
#'   \item{LA_Estimated}{Leaf area from the picture (m^2^)}
#'   \item{Date}{Date of measurement}
#'   \item{pos}{unique pot using rowcol coordinate}
#' }
#'
#' @source \url{}
#'
#' @references Cabrera-Bosquet, Llorenç, Fournier Christian, Brichet Nicolas,
#' Welcker Claude, Suard Benoît, and Tardieu François. 2016.
#' “High-throughput estimation of incident light, light interception and
#' radiation-use efficiency of thousands of plants in a phenotyping platform.”
#' New Phytologist 212 (1): 269-81.
#' doi:10.1111/nph.14027.
"PhenoarchDat1"
