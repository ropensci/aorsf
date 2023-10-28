## Version 0.1.1

Update, October 25: Thank you for helping me with this. I have tidied up threads and avoided usage of the `data.table` functions that were creating possible memory leaks.  I have checked this submission locally with valgrind and on rhub, with both indicating 0 memory leaks. However, if this submission does not pass on your end, I would like to request an extension on the October 28th deadline.

Update, October 21: I have updated the submission to fix memory leaks. Many of the leaks were caused by my omission of a virtual de-constructor for derived classes or by omission of a delete statement for dynamically allocated memory. I apologize for these oversights. After reviewing, you may still see a possible memory leak from `orsf_ice` functions. From what I can tell, this possible leak could either be measurement error or could be attributed to `data.table`. I do not think it's from `aorsf`.

Initial submission: This version is being submitted to CRAN early due to a memory error that was identified in version 1.0.0. I apologize for the oversight. As `aorsf` would be removed from CRAN if the issue is not fixed before October 28, I would like to request an expedited submission. I have run the current submission's tests and examples with valgrind to ensure the memory error has been fixed. 

## R CMD check results

Duration: 4m 8.3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded

## valgrind

The error spotted by valgrind on initial submission of 0.1.0, "a conditional jump or move depends on uninitialized values has been fixed".

## Downstream dependencies

I have also run R CMD check on downstream dependencies of `aorsf`: 

- Rcpp
- data.table. 
- collapse

All packages passed.

## Version 0.1.0

## R CMD check results

Duration: 3m 53.1s

❯ checking C++ specification ... NOTE
    Specified C++14: please drop specification unless essential

0 errors ✔ | 0 warnings ✔ | 1 note ✖

I have specified C++14 for this release. C++14 is essential, as this release uses `std::make_unique`.

## Downstream dependencies

I have also run R CMD check on downstream dependencies of `aorsf`: 

- Rcpp
- data.table. 
- collapse

All packages passed.

## Version 0.0.7

I noticed the package was removed from CRAN on 1-10 so I made some additional changes and am now re-submitting as these changes should fix the problems that appear on the ATLAS check page. 

## Version 0.0.6

I reviewed the additional issues at https://cran.r-project.org/web/checks/check_results_aorsf.html and modified my testing scripts based on results appearing from ATLAS checks.

## Version 0.0.5, submitted 2022-12-13

## Version 0.0.4, submitted 2022-11-06

## Version 0.0.3, submitted 2022-10-07

R CMD check succeeded with no warnings or notes for `aorsf` and its dependencies.

## Version 0.0.2, submitted 2022-08-24

I reviewed the additional issues at https://cran.r-project.org/web/checks/check_results_aorsf.html and modified my testing scripts based on results appearing from ATLAS and MKL checks.

Updates: 

- added tolerance for the equivalence test in `test-orsf_predict.R` that failed on ATLAS.

- relegated time to train tests (`test-orsf.R`, line 617-655) to development version only. These tests can fail randomly (as it did on MKL) due to variations in computing environment.

## Version 0.0.1, submitted 2022-08-22

Updates: 

- ensured all links use https instead of http and have added trailing slashes where appropriate.

- added a reference to `DESCRIPTION` that discusses the methods used in `aorsf`.

- Extra spaces in `DESCRIPTION` removed.

- `LICENSE` file fixed and names of authors added.

- the trailing slash on URL's to the arXiv paper was removed in the readme.

- the trailing slash on URL's to the arXiv paper was removed in the aorsf vignette.
