
## R CMD check results 

Duration: 2m 31.7s

0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

## Downstream dependencies

I have also run R CMD check on downstream dependencies of `aorsf`: 

- Rcpp
- data.table. 

Both packages passed.

## Revision 1, 2022-08-22

Updates: 

- ensured all links use https instead of http and have added trailing slashes where appropriate.

- added a reference to `DESCRIPTION` that discusses the methods used in `aorsf`.


## Revision 2, 2022-08-22

Updates: 

- the trailing slash on URL's to the arXiv paper was removed in the readme.

## Revision 3, 2022-08-22

Updates: 

- the trailing slash on URL's to the arXiv paper was removed in the aorsf vignette.

## Revision 4, 2022-08-23

- Extra spaces in `DESCRIPTION` removed.

- `LICENSE` file fixed and names of authors added.

## Revision 4, 2022-08-24

I reviewed the additional issues at https://cran.r-project.org/web/checks/check_results_aorsf.html and modified my testing scripts based on results appearing from ATLAS and MKL checks.

Updates: 

- added tolerance for the equivalence test in `test-orsf_predict.R` that failed on ATLAS.

- relegated time to train tests (`test-orsf.R`, line 617-655) to development version only. These tests can fail randomly (as it did on MKL) due to variations in computing environment.


