% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{RootDat1}
\alias{RootDat1}
\title{Greenhouse data for an experiment in the RootPhAir platform.}
\format{
A data.frame with 16,275 rows and 10 columns:
\describe{
\item{Exp}{Experiment number}
\item{thermalTime}{Thermal time cumulated}
\item{Genotype}{Genotype}
\item{plantId}{Unique pot using tank and rowcol coordinate}
\item{Tank}{Tank A or B}
\item{Strip}{Number of strip of five plants (i.e. row coordinate)}
\item{Pos}{Position within th strip (i.e. column coordinate)}
\item{tipPos_x}{Position of the root tip in x axis}
\item{tipPos_y}{Position of the root tip in y axis}
\item{Time}{Time of measurement}
}
}
\usage{
RootDat1
}
\description{
A dataset containing greenhouse data from the RootPhAir platform
(UCLouvain, Belgium). It consists of one experiment with one
aeroponic tanks with 340 maize plants. The studied traits are the root tip
coordinates in y and x axis, extracted from the pictures over time.
Plants were pictured every 2 hours for 10 days.
This dataset was kindly provided by Xavier Draye.
}
\keyword{datasets}
