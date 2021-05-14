# statgenHTP 1.0.1.1

* Plotting second derivatives of fitted splines, plotted first derivatives. This is fixed now. 
* The ```detectSerieOut```` function now uses an extra criterion for checking if time courses our outlying. See the function documentation and vignettes for a full explanation of this new criterion.
* The parameter ```trait``` in the ```removeSerieOut``` function is renamed to ```traits``` and now accepts a vector of traits that for which outlier values can be replaced by ```NA```.
* Plotting the output of ```detectSerieOut``` now has an extra parameter ```geno.decomp``` to restrict the output to selected levels of the geno.decomp variable in the data.

# statgenHTP 1.0.1

* No user visible changes

# statgenHTP 1.0.0

* Initial CRAN Release
