## Update to fix additional issues in checks

Release adding functionality to existing functions by extending the values some parameters can take and by adding new parameters. In addition so minor bug fixes and additional tests.

----

## Test environments

* local Windows 10 install, R 4.1
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

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

