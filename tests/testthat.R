library(testthat)
library(aorsf)

## force tests to be executed if in dev release which we define as
## having a sub-release, eg 0.9.15.5 is one whereas 0.9.16 is not
if (length(strsplit(packageDescription("aorsf")$Version, "\\.")[[1]]) > 3) {
 Sys.setenv("run_all_aorsf_tests" = "yes")
} else {

}

test_check("aorsf")
