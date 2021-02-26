## Initial release

----

## Test environments

* local Windows 10 install, R 4.0.2
* local Debian
* R-hub (devel and release)
* winbuilder (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTES:

  * checking CRAN incoming feasibility ... NOTE    
  
    Possibly mis-spelled words in DESCRIPTION:
    EPPN (37:67)
    HTP (2:37)
    Phenotyping (2:24)
    phenotyping (39:41)
    statgenHTP (37:18)
    
    - These are all spelled correctly.

    Suggests or Enhances not in mainstream repositories: asreml

    - asreml is a commercial R package that is used as one of three alternatives for modeling data.

  * checking package dependencies ... NOTE  
    Package suggested but not available for checking: 'asreml'
    
    - asreml is a commercial R package that is used as one of three alternatives for      modeling data.

