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













#' @srrstatsTODO {ML5.1} *As for the untrained model objects produced according to the above standards, and in particular as a direct extension of **ML3.3**, the properties and behaviours of trained models produced by ML software should be explicitly compared with equivalent objects produced by other ML software. (Such comparison will generally be done in terms of comparing model performance, as described in the following standard **ML5.3**--**ML5.4**).*
#' @srrstatsTODO {ML5.2} *The structure and functionality of objects representing trained ML models should be thoroughly documented. In particular,*
#' @srrstatsTODO {ML5.2a} *Either all functionality extending from the class of model object should be explicitly documented, or a method for listing or otherwise accessing all associated functionality explicitly documented and demonstrated in example code.*
#' @srrstatsTODO {ML5.2b} *Documentation should include examples of how to save and re-load trained model objects for their re-use in accordance with **ML3.1**, above.*
#' @srrstatsTODO {ML5.2c} *Where general functions for saving or serializing objects, such as [`saveRDS`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html) are not appropriate for storing local copies of trained models, an explicit function should be provided for that purpose, and should be demonstrated with example code.*
#' @srrstatsTODO {ML5.3} *Assessment of model performance should be implemented as one or more functions distinct from model training.*
#' @srrstatsTODO {ML5.4} *Model performance should be able to be assessed according to a variety of metrics.*
#' @srrstatsTODO {ML5.4a} *All model performance metrics represented by functions internal to a package must be clearly and distinctly documented.*
#' @srrstatsTODO {ML5.4b} *It should be possible to submit custom metrics to a model assessment function, and the ability to do so should be clearly documented including through example code.*
#' @srrstatsTODO {ML6.0} *Descriptions of ML software should make explicit reference to a workflow which separates training and testing stages, and which clearly indicates a need for distinct training and test data sets.*
#' @srrstatsTODO {ML6.1} *ML software intentionally designed to address only a restricted subset of the workflow described here should clearly document how it can be embedded within a typical full ML workflow in the sense considered here.*
#' @srrstatsTODO {ML6.1a} *Such demonstrations should include and contrast embedding within a full workflow using at least two other packages to implement that workflow.*
#' @srrstatsTODO {ML7.0} *Test should explicitly confirm partial and case-insensitive matching of "test", "train", and, where applicable, "validation" data.*
#' @srrstatsTODO {ML7.1} *Tests should demonstrate effects of different numeric scaling of input data (see **ML2.2**).*
#' @srrstatsTODO {ML7.2} *For software which imputes missing data, tests should compare internal imputation with explicit code which directly implements imputation steps (even where such imputation is a single-step implemented via some external package). These tests serve as an explicit reference for how imputation is performed.*
#' @srrstatsTODO {ML7.3} *Where model objects are implemented as distinct classes, tests should explicitly compare the functionality of these classes with functionality of equivalent classes for ML model objects from other packages.*
#' @srrstatsTODO {ML7.3a} *These tests should explicitly identify restrictions on the functionality of model objects in comparison with those of other packages.*
#' @srrstatsTODO {ML7.3b} *These tests should explicitly identify functional advantages and unique abilities of the model objects in comparison with those of other packages.*
#' @srrstatsTODO {ML7.4} *ML software should explicit document the effects of different training rates, and in particular should demonstrate divergence from optima with inappropriate training rates.*
#' @srrstatsTODO {ML7.5} *ML software which implements routines to determine optimal training rates (see **ML3.4**, above) should implement tests to confirm the optimality of resultant values.*
#' @srrstatsTODO {ML7.6} *ML software which implement independent training "epochs" should demonstrate in tests the effects of lesser versus greater numbers of epochs.*
#' @srrstatsTODO {ML7.7} *ML software should explicitly test different optimization algorithms, even where software is intended to implement one specific algorithm.*
#' @srrstatsTODO {ML7.8} *ML software should explicitly test different loss functions, even where software is intended to implement one specific measure of loss.*
#' @srrstatsTODO {ML7.9} *Tests should explicitly compare all possible combinations in categorical differences in model architecture, such as different model architectures with same optimization algorithms, same model architectures with different optimization algorithms, and differences in both.*
#' @srrstatsTODO {ML7.9a} *Such combinations will generally be formed from multiple categorical factors, for which explicit use of functions such as [`expand.grid()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/expand.grid.html) is recommended.*
#' @srrstatsTODO {ML7.10} *The successful extraction of information on paths taken by optimizers (see **ML5.1**, above), should be tested, including testing the general properties, but not necessarily actual values of, such data.*
#' @srrstatsTODO {ML7.11} *All performance metrics available for a given class of trained model should be thoroughly tested and compared.*
#' @srrstatsTODO {ML7.11a} *Tests which compare metrics should do so over a range of inputs (generally implying differently trained models) to demonstrate relative advantages and disadvantages of different metrics.*

#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#'
#' @srrstatsNA {G2.14}, @srrstatsNA {G2.14a}, @srrstatsNA {G2.14b}, @srrstatsNA {G2.14c} *I have made orsf and its associated functions throw an error when there is a missing value in the relevant data. Here is why I made this decision: (1) imputation of missing data is an involved process that many other packages have been designed to engage with. I want aorsf to be good at one thing, which is oblique random survival forests. If I try to make routines to handle missing data, I am kind of re-inventing the wheel when I could be working on things that are more relevant to the oblique random survival forest. (2) ignoring missing data would be okay from a programmatic point of view, but I have chosen not to implement this because it would not be helpful if the user was unaware of their missing data. I want orsf to perform a hard stop when it detects missing data because in many cases, junior analysts are not familiar enough with their data to know it has missing values, and perpetuating the unawareness of missing data by handling it on the back-end of analysis functions just creates downstream issues when the analysis is written up.*
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

#'
#' @noRd
NULL

