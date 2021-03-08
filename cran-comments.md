## Update to fix additional issues in checks

Fixed test issues with noLD, OpenBLAS and MKL

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
  
       Days since last update: 3
       
    - In the initial version of the package one test didn't pass the noLD, OpenBLAS and MKL checks after being accepted. This is fixed now.
   
   Suggests or Enhances not in mainstream repositories:
     asreml

    - asreml is a commercial R package that is used as one of the alternatives for modeling data.

  * checking package dependencies ... NOTE  
    Package suggested but not available for checking: 'asreml'
    
    - asreml is a commercial R package that is used as one of the alternatives
    for modeling data.

