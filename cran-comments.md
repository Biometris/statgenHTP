## Patch release fixing MKL check error

Unable to produce locally or using rhub, but highly likely fixed

----

## Test environments

* local Windows 11 install, R 4.6.1
* Ubuntu (on github actions, devel, release, and oldrelease)
* macOS (on github actions, release)
* mkl, intel, atlas (using rhub/github actions)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTES:

  * checking CRAN incoming feasibility ... NOTE    

   Suggests or Enhances not in mainstream repositories:
     asreml

    - asreml is a commercial R package that is used as one of the alternatives for modeling data.

  * checking package dependencies ... NOTE  
    Package suggested but not available for checking: 'asreml'
    
    - asreml is a commercial R package that is used as one of the alternatives
    for modeling data.

