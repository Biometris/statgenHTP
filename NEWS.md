# statgenHTP 1.0.3

* The `removeSerieOut` function now has an extra argument reason allowing for restricting removal of outliers to one or more reason the outliers where tagged.

# statgenHTP 1.0.3

* No user visible changes

# statgenHTP 1.0.2

* Plotting second derivatives of fitted splines, plotted first derivatives. This is fixed now. 
* The `detectSerieOut` function now uses an extra criterion for checking if time courses our outlying. See the function documentation and vignettes for a full explanation of this new criterion.
* The parameter `trait` in the `removeSerieOut` function is renamed to `traits` and now accepts a vector of traits that for which outlier values can be replaced by `NA`.
* Plotting the output of `detectSerieOut` now has an extra parameter `geno.decomp` to restrict the output to selected levels of the geno.decomp variable in the data.
* Computing the area under the curve in `estimateSplineParameters` now allows the specification of the time unit used on the x-axis.
* Fitting models using SpATS is now done using the option `centered = TRUE` by default.

# statgenHTP 1.0.1

* No user visible changes

# statgenHTP 1.0.0

* Initial CRAN Release
