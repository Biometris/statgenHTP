# statgenHTP 1.0.6.1

* Functions no longer rely on soft-deprecated ggplot2 functions.

# statgenHTP 1.0.6

* The package has been extended by functions for modelling evolution of the genetic signal. This extension consists of functions `fitSplineHDM()` for fitting a hierarchical data model, `predict.psHDM()` for making predictions based on this model, and `plot.psHDM()` for plotting the results. The methods used are described in a new vignette.
* A small bug in `detectSingleOutMaize()` is fixed. Observations with a missing value for one of the involved traits are no longer tagged as outliers.

# statgenHTP 1.0.5

* No user visible changes

# statgenHTP 1.0.4

* The `removeSerieOut()` function now has an extra argument reason allowing for restricting removal of outliers to one or more reason the outliers where tagged.
* The `plot` function for `serieOut` objects now has an extra argument reason allowing for restricting the plotting of outliers to one or more reason the outliers where tagged.
* The `detectSerieOut()` function is now able to handle plotIds with irregular naming, i.e. plotIds starting with a number.
* The results of `estimateSplineParameters()` can now be plotted in a box plot and histogram.
* In the `estimateSplineParameters()` function multiple parameters can now be estimated at once.
* A bug in `detectSerieOut()` that caused slope outliers to be never detected is fixed.
* A bug causing predictions to be made in `fitSpline()` for missing values at the beginning or end of a time course is fixed.
* The `detectSerieOut()` function now checks for number of plotIds per plant in the correct location leading to a more user friendly error message.

# statgenHTP 1.0.3

* No user visible changes

# statgenHTP 1.0.2

* Plotting second derivatives of fitted splines, plotted first derivatives. This is fixed now. 
* The `detectSerieOut()` function now uses an extra criterion for checking if time courses our outlying. See the function documentation and vignettes for a full explanation of this new criterion.
* The parameter `trait` in the `removeSerieOut()` function is renamed to `traits` and now accepts a vector of traits that for which outlier values can be replaced by `NA`.
* Plotting the output of `detectSerieOut()` now has an extra parameter `geno.decomp` to restrict the output to selected levels of the geno.decomp variable in the data.
* Computing the area under the curve in `estimateSplineParameters()` now allows the specification of the time unit used on the x-axis.
* Fitting models using SpATS is now done using the option `centered = TRUE` by default.

# statgenHTP 1.0.1

* No user visible changes

# statgenHTP 1.0.0

* Initial CRAN Release
