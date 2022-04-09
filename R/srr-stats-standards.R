#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'





















#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#'
#' @srrstatsNA {G2.4c} *orsf() and other functions in aorsf generally throw errors when data types do not meet expectations. I believe this is the design I want to pursue because I think trusting users to learn how to supply inputs correctly will lead to more accurate and reproducible science than relying on my software to continuously fix mistakes in the inputs.*
#' @srrstatsNA {G2.4d} see above
#' @srrstatsNA {G2.4e} see above
#' @srrstatsNA {G2.12} data.frame or list columns are not accepted as valid types for orsf().
#'
#' @srrstatsNA {G2.14} *I have made orsf and its associated functions throw an error when there is a missing value in the relevant data. Here is why I made this decision: (1) imputation of missing data is an involved process that many other packages have been designed to engage with. I want aorsf to be good at one thing, which is oblique random survival forests. If I try to make routines to handle missing data, I am kind of re-inventing the wheel when I could be working on things that are more relevant to the oblique random survival forest. (2) ignoring missing data would be okay from a programmatic point of view, but I have chosen not to implement this because it would not be helpful if the user was unaware of their missing data. I want orsf to perform a hard stop when it detects missing data because in many cases, junior analysts are not familiar enough with their data to know it has missing values, and perpetuating the unawareness of missing data by handling it on the back-end of analysis functions just creates downstream issues when the analysis is written up.*
#' @srrstatsNA {G2.14a} see above
#' @srrstatsNA {G2.14b} see sbove
#' @srrstatsNA {G2.14c} see above
#'
#' @srrstatsNA {G3.1} *there is no user-facing covariance calculation.*
#'
#' @srrstatsNA {G3.1a} *Specific covariance methods are not applied, although there are explicitly documented control methods for fitting oblique survival trees.*
#'
#' @srrstatsNA {G4.0} *outputs are not written to local files.*
#'
#' @srrstatsNA {G5.1} *No data sets are created within the package. Random values are created from a set seed, and they can be reproduced this way if needed*
#'
#' @srrstatsNA {G5.4c} *There are no tests where published paper outputs are applicable and where code from original implementations is not available*
#'
#' @srrstatsNA {G5.6b} *Parameter recovery tests are not run with multiple random seeds as neither data simulation nor the algorithms we test contain a random component.*
#'
#' @srrstatsNA {G5.10} *Extended tests are not used in this package*
#' @srrstatsNA {G5.11} *Extended tests are not used in this package*
#' @srrstatsNA {G5.11a} *Extended tests are not used in this package*
#' @srrstatsNA {G5.12} *Extended tests are not used in this package*
#'
#' @srrstatsNA {ML1.0a} *Training and testing data are not eschewed.*
#' @srrstatsNA {ML1.1a, ML1.1b} *The labels 'train' and 'test' are not used by functions in aorsf, so are not explicitly confirmed via pre-processing steps.*
#'
#' @srrstatsNA {ML1.2} *orsf() only requires training data, and predict.aorsf() only requires testing data, so I do not not apply methods to distinguish training and testing data in these functions.*
#'
#' @srrstatsNA {ML1.4} *aorsf does not provide explicit training and test data sets, only a small dataset for illustrative purposes. If testing data are needed for illustration, they are created by sub-setting.*
#'
#' @srrstatsNA {ML1.5} *aorsf is a very lightweight package that focuses on fitting, applying, and interpreting oblique random survival forests. There are numerous other R packages that summarize the contents of data sets, and their extensively developed functions are much better at summarizing training and testing data than something I would write*
#'
#' @srrstatsNA {ML1.7, ML1.7a, ML1.7b, ML1.8} *aorsf does not admit or impute missing values*
#'
#' @srrstatsNA {ML2.1} *aorsf does not use broadcasting to reconcile dimensionally incommensurate input data.*
#' @srrstatsNA {ML2.2, ML2.2a, ML2.2b} *aorsf does not require numeric transformation of input data and therefore does not have a dedicated input data specification stage. When internal transformations are performed, they are always reversed for compatibility with the original input. The strategy used by aorsf is the same as Terry Therneau's routine in the survival::coxph function. I do not intend to allow specification of target values for this transformation. I want to keep my Newton-Raphson scoring procedure as close to identical to Terry Therneau's as possible.*
#'
#' @srrstatsNA {ML3.0b} *Since aorsf doesn't have a dedicated input data specification state, the output of orsf(no_fit=TRUE) is passed directly to orsf_train().*
#'
#' @srrstatsNA {ML3.1} *aorsf allows users to print and fit untrained models. Since aorsf doesn't include a data pre-processing stage, I am not sure if there are any other helpful functions to include for an untrained random forest. For pre-trained models, users can always run saveRDS and readRDS on aorsf objects, which are just lists with attributes. Thus, additional functions for these purposes have not been added*
#'
#' @srrstatsNA {ML3.2} *aorsf does not have a dedicated input data specification step.*
#'
#' @srrstatsNA {ML3.4, ML3.4a, ML3.4b} *aorsf does not use training rates.*
#'
#' @srrstatsNA {ML3.7} *This software uses C++, facilitated through Rcpp, which does not currently allow user-controlled use of either CPUs or GPUs.*
#'
#' @srrstatsNA {ML4.1c} *The random forest trees do not depend on each other, so there is no information used to advance from one tree to the next.*
#'
#' @srrstatsNA {ML4.3} *aorsf does not use batch processing.*
#'
#' @srrstatsNA {ML4.4} *aorsf does not use batch processing.*
#'
#' @srrstatsNA {ML4.6} *aorsf does not use batch jobs.*
#'
#' @srrstatsNA {ML4.7, ML4.8, ML4.8a} *aorsf does not currently include functions for re-sampling as random forests generally do not need a lot of tuning and because there are other R packages that are dedicated to providing robust resampling routines (e.g., rsample).*

#' @srrstatsNA {ML5.1} *The properties and behaviours of ORSF models were explicitly compared with objects produced by other ML software in Jaeger et al, 2019 (DOI: 10.1214/19-AOAS1261). These comparisons focused on comparing model performance. I am not including comparisons such as this in the aorsf package because I want aorsf to include or suggest including as few other R packages as possible.*
#'
#' @srrstatsNA {ML5.2c} *General functions for saving or serializing objects, such as [`saveRDS`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) are  appropriate for storing local copies of trained aorsf models.*
#'
#' @srrstatsNA {ML7.0} *aorsf does not have text inputs with labels of "test", "train", or "validation" data. However, aorsf does implement a function, check_arg_is_valid(), which assesses validity of text inputs based on a case-sensitive set of valid options.*
#'
#' @srrstatsNA {ML7.2} *aorsf does not impute missing data.*
#'
#' @srrstatsNA {ML7.3, ML7.3a, ML7.3b} *I am not including comparisons such as this in the aorsf package because I want aorsf to include or suggest including as few other R packages as possible. I don't want to overload the imported or suggested packages for aorsf because it becomes exponentially harder to get a package onto CRAN the more it depends on other packages. Jaeger et al, 2019 (DOI: 10.1214/19-AOAS1261) made comparisons like these formally using several ML software packages, and I plan on writing a similar paper for aorsf that will make meaningful comparisons similar to the ones I made in Jaeger et al, 2019. *
#'
#' @srrstatsNA {ML7.4} *aorsf does not use training rates*
#'
#' @srrstatsNA {ML7.5} *aorsf does not use training rates*
#'
#' @srrstatsNA {ML7.6} *aorsf does not use training epochs.*
#'
#' @srrstatsNA {ML7.11a} *aorsf does not implement multiple metrics and therefore cannot demonstrate relative advantages and disadvantages of different metrics. However, when verifying the accuracy of aorsf's scripts to compute certain metrics (e.g., the likelihood ratio test, and cox PH regression), aorsf tests to make sure the metrics are correctly computed over a wide range of inputs*
#'
#' @noRd
NULL
