% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectSingleOut.R
\name{detectSingleOut}
\alias{detectSingleOut}
\title{Detect outliers for single observations}
\usage{
detectSingleOut(
  TP,
  trait,
  plotIds = NULL,
  checkEdges = TRUE,
  confIntSize = 5,
  nnLocfit = 0.5
)
}
\arguments{
\item{TP}{An object of class \code{TP}.}

\item{trait}{A character vector indicating the trait to model in \code{TP}.}

\item{plotIds}{A character vector of plotIds for which the outliers should be
detected. If \code{NULL}, all plotIds in \code{TP} are used.}

\item{checkEdges}{Before fitting the local regression should a check be done
if the first and last time point for a plot are outlying observations?}

\item{confIntSize}{A numeric value defining the confidence interval (see
Details).}

\item{nnLocfit}{A numeric value defining the constant component of the
smoothing parameter nn (see Details).}
}
\value{
An object of class singleOut, a \code{data.frame} with the following
columns.
\describe{
\item{plotId}{plotId}
\item{timePoint}{time point}
\item{trait}{modeled trait}
\item{yPred}{prediction from the local regression}
\item{sd_yPred}{standard deviation of the prediction}
\item{lwr}{lower bound of the confidence interval}
\item{upr}{upper bound of the confidence interval}
\item{outlier}{flag for detected outlier (a value of 1 indicates the
observation is an outlier)}
}
}
\description{
Detect outlying observations in a time series by modeling each plotId using
a local regression.
}
\details{
See locfit() help function from the locfit R library. The user can act on:
\describe{
\item{nnLocfit}{the constant of the smoothing parameter. Increase nnLocfit
to have a very smooth curve}
\item{confIntSize}{the level to calculate the confidence interval. Increase
confIntSize to exclude less outliers}
}
}
\examples{
## Create a TP object containing the data from the Phenovator.
PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos \%in\%
                                 c("c24r41", "c7r18", "c7r49"), ]
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

## First select a subset of plants, for example here 9 plants
plantSel <- phenoTP[[1]]$plotId[1:9]
# Then run on the subset
resuVatorHTP <- detectSingleOut(TP = phenoTP,
                                trait = "EffpsII",
                                plotIds = plantSel,
                                confIntSize = 3,
                                nnLocfit = 0.1)

}
\seealso{
Other functions for detecting outliers for single observations: 
\code{\link{detectSingleOutMaize}()},
\code{\link{plot.singleOut}()},
\code{\link{removeSingleOut}()}
}
\concept{functions for detecting outliers for single observations}
