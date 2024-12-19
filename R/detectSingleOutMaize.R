#' detectSingleOutMaize
#'
#' Function to detect plant outliers in a temporal lattice experiment on Maize
#' which can be extended to other experiment types.
#' The criteria needs three phenotypes (ex for maize: the estimated biomass,
#' plant height and phyllocron)
#' \describe{
#'   \item{plants are identified as "small outlier plant"}{if for biomass AND
#'         phyllocron \eqn{res_i < \mu_{res} - qnorm(threshold) * sd_{res}}}
#'   \item{plants are identified as "big outlier plant"}{ if for biomass AND
#'         plant height \eqn{res_i > \mu_{res} + qnorm(threshold) * sd_{res}}}
#' }
#'
#' @param TP An object of class TP.
#' @param timeBeforeTrt A character or numeric value indicating the date just
#' before treatment in the experiment. When using a character string to
#' reference a time point, the value has to be an exact match to one of the
#' existing timePoints. When using a number it will be matched by its number
#' ("timeNumber") in the timePoints attribute of the TP object.
#' @param trait1 A character vector indicating the first trait to model in TP.
#' @param trait2 A character vector indicating the second trait to model in TP.
#' @param trait3 A character vector indicating the third trait to model in TP.
#' @param thr A numeric value indicating the threshold.
#'
#' @returns A list with three data.frames, \code{modDat} containing the
#' data used for fitting the models, \code{smallPlants} containing the plants
#' identified as small plants and \code{bigPlants} containing the plants
#' identified as big plants.
#'
#' @examples
#' \donttest{
#' ## Create a TP object containing the data from the PhenoArch.
#' phenoTParch <- createTimePoints(dat = PhenoarchDat1,
#'                                 experimentName = "Phenoarch",
#'                                 genotype = "Genotype",
#'                                 timePoint = "Date",
#'                                 plotId = "pos",
#'                                 rowNum = "Row",
#'                                 colNum = "Col")
#' singleOutMaize <- detectSingleOutMaize(phenoTParch,
#'                                        timeBeforeTrt = "2017-04-27",
#'                                        trait1 = "Biomass",
#'                                        trait2 = "PlantHeight",
#'                                        trait3 = "phyllocron",
#'                                        thr = 0.95)
#' }
#'
#' @family functions for detecting outliers for single observations
#'
#' @export
detectSingleOutMaize <- function(TP,
                                 timeBeforeTrt,
                                 trait1 = "Biomass",
                                 trait2 = "PlantHeight",
                                 trait3 = "phyllocron",
                                 thr = 0.95) {
  ## Checks.
  if (!inherits(TP, "TP")) {
    stop("TP should be an object of class TP.\n")
  }
  timePoint <- chkTimePoints(TP, timeBeforeTrt)
  modDat <- TP[[timePoint]]
  ## Fit spatial models for the three traits.
  modTr1 <- fitModels(TP, trait = trait1, timePoints = timeBeforeTrt,
                      quiet = TRUE)
  modTr2 <- fitModels(TP, trait = trait2, timePoints = timeBeforeTrt,
                      quiet = TRUE)
  modTr3 <- fitModels(TP, trait = trait3, timePoints = timeBeforeTrt,
                      quiet = TRUE)
  ## Add residuals and fitted values to modDat.
  modDat[["fittedTr1"]] <- fitted(modTr1[[1]])
  modDat[["devTr1"]] <- resid(modTr1[[1]])
  modDat[["fittedTr2"]] <- fitted(modTr2[[1]])
  modDat[["devTr2"]] <- resid(modTr2[[1]])
  modDat[["fittedTr3"]] <- fitted(modTr3[[1]])
  modDat[["devTr3"]] <- resid(modTr3[[1]])
  ## Calculate outlier criteria.
  qVal <- qnorm(thr)
  modDat[["meanDevTr1"]] <- mean(modDat[["devTr1"]], na.rm = TRUE)
  modDat[["sdDevTr1"]] <- sd(modDat[["devTr1"]], na.rm = TRUE)
  modDat[["meanDevTr2"]] <- mean(modDat[["devTr2"]], na.rm = TRUE)
  modDat[["sdDevTr2"]] <- sd(modDat[["devTr2"]], na.rm = TRUE)
  modDat[["meanDevTr3"]] <- mean(modDat[["devTr3"]], na.rm = TRUE)
  modDat[["sdDevTr3"]] <- sd(modDat[["devTr3"]], na.rm = TRUE)
  modDat[["lowDevTr1"]] <- modDat[["meanDevTr1"]] - modDat[["sdDevTr1"]] * qVal
  modDat[["lowDevTr2"]] <- modDat[["meanDevTr2"]] - modDat[["sdDevTr2"]] * qVal
  modDat[["lowDevTr3"]] <- modDat[["meanDevTr3"]] - modDat[["sdDevTr3"]] * qVal
  modDat[["upDevTr1"]] <- modDat[["meanDevTr1"]] + modDat[["sdDevTr1"]] * qVal
  modDat[["upDevTr2"]] <- modDat[["meanDevTr2"]] + modDat[["sdDevTr2"]] * qVal
  modDat[["upDevTr3"]] <- modDat[["meanDevTr3"]] + modDat[["sdDevTr3"]] * qVal
  modDat[["lowSpTr1"]] <- ifelse(modDat[["devTr1"]] > modDat[["lowDevTr1"]], 1, 0)
  modDat[["lowSpTr2"]] <- ifelse(modDat[["devTr2"]] > modDat[["lowDevTr2"]], 1, 0)
  modDat[["lowSpTr3"]] <- ifelse(modDat[["devTr3"]] > modDat[["lowDevTr3"]], 1, 0)
  modDat[["upSpTr1"]] <- ifelse(modDat[["devTr1"]] < modDat[["upDevTr1"]], 1, 0)
  modDat[["upSpTr2"]] <- ifelse(modDat[["devTr2"]] < modDat[["upDevTr2"]], 1, 0)
  modDat[["upSpTr3"]] <- ifelse(modDat[["devTr3"]] < modDat[["upDevTr3"]], 1, 0)
  ## Remove missing values for response variables.
  ## Leaving them in gives NA only rows in small and big plants.
  modDat <- modDat[!is.na(modDat[[trait1]]) & !is.na(modDat[[trait2]]) &
                     !is.na(modDat[[trait3]]), ]
  ## Detect small plants.
  modDat[["smallPlant"]] <-
    ifelse(modDat[["lowSpTr1"]] + modDat[["lowSpTr3"]] == 0, 1, 0)
  ## Detect big plants.
  modDat[["bigPlant"]] <-
    ifelse(modDat[["upSpTr1"]] + modDat[["upSpTr2"]] == 0, 1, 0)
  ## Create data.frames with outliers.
  smallPlants <- modDat[modDat[["smallPlant"]] == 1, ]
  bigPlants <- modDat[modDat[["bigPlant"]] == 1, ]
  return(list(modDat = modDat, smallPlants = smallPlants, bigPlants = bigPlants))
}

