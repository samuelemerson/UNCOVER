# UNCOVER release notes

## V1.0.0

* Initial release to CRAN

## V1.1.0

* Fixed bug which caused an error if a deforestation criterion only finds one cluster before the deforestation stage 

* Fixed bug with memoisation and allowed a parameter to be specified which gives the user control of when to use RIBIS

* Changed name of balanced deforestation criterion to "diverse"

* Allowed UNCOVER to predict a single observation

* Accounted for duplicates in the data by stacking duplicates when constructing
the graph
