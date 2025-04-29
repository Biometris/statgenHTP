## Patch release

- Minor bugfix and documentation improvements.

----

## Test environments

* local Windows 11 install, R 4.5.0
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There wer 2 NOTES:

  * checking CRAN incoming feasibility ... NOTE    

   Suggests or Enhances not in mainstream repositories:
     asreml

    - asreml is a commercial R package that is used as one of the alternatives for modeling data.

  * checking package dependencies ... NOTE  
    Package suggested but not available for checking: 'asreml'
    
    - asreml is a commercial R package that is used as one of the alternatives
    for modeling data.

