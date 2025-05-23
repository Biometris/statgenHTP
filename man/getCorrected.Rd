% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMethods.R
\name{getCorrected}
\alias{getCorrected}
\title{Extract corrected phenotypic values}
\usage{
getCorrected(fitMod, timePoints = names(fitMod), outFile = NULL)
}
\arguments{
\item{fitMod}{An object of class \code{fitMod}.}

\item{timePoints}{A character or numeric vector indicating the time point(s)
for which the corrected values should be extracted. When using a character
string to reference a time point, the value has to be an exact match to one
of the existing time points. When using a number it will be matched by its
number ("timeNumber") in the timePoints attribute of the TP object.}

\item{outFile}{A character string indicating the .csv file to which the
results should be written. If \code{NULL} no file is written.}
}
\value{
A data.frame with spatially corrected values per time point.
}
\description{
Extract corrected phenotypic values from an object of class fitMod. After
fitting a spatial model at each time point, the raw phenotypic data is
corrected by subtracting the (estimated) sources of variation (environmental,
design  effect) that are of no interest (nuisances). This allows keeping
the data resolution at the plot/plant level.
}
\examples{
\donttest{
## Using the first example dataset (PhenovatorDat1).
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            experimentName = "Phenovator",
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            repId = "Replicate",
                            plotId = "pos",
                            rowNum = "y", colNum = "x",
                            addCheck = TRUE,
                            checkGenotypes = c("check1", "check2",
                                               "check3", "check4"))

## Fit a SpATS model on few time points.
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        timePoints = c(1, 6, 20))

## Extract the corrected values for one time point:
spatCorrSp <- getCorrected(modPhenoSp,
                           timePoints = 6)
head(spatCorrSp)
}

}
\seealso{
Other functions for spatial modeling: 
\code{\link{fitModels}()},
\code{\link{getEffDims}()},
\code{\link{getGenoPred}()},
\code{\link{getHerit}()},
\code{\link{getVar}()},
\code{\link{plot.fitMod}()},
\code{\link{summary.fitMod}()}
}
\concept{functions for spatial modeling}
