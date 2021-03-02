## Initial release

After earlier submission, the 5 vignettes are now condensed to a single overview vignette. Still takes about a minute to check but I hope this is acceptable. This also reduced the package size from 8Mb to 3.5Mb.

----

## Test environments

* local Windows 10 install, R 4.0.4
* local Debian, R 3.6.3
* R-hub (devel and release)
* winbuilder (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTES:

  * checking CRAN incoming feasibility ... NOTE    
  
    New submission
  
    Possibly mis-spelled words in DESCRIPTION:
    EPPN (40:16)
    HTP (2:37, 38:18)
    Phenotyping (2:24)
    phenotyping (38:5)

    - These are all spelled correctly.

    Suggests or Enhances not in mainstream repositories: asreml

    - asreml is a commercial R package that is used as one of the alternatives for modeling data.

  * checking package dependencies ... NOTE  
    Package suggested but not available for checking: 'asreml'
    
    - asreml is a commercial R package that is used as one of the alternatives
    for modeling data.

