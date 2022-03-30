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


#' @srrstatsTODO {ML1.0} *Documentation should make a clear conceptual distinction between training and test data (even where such may ultimately be confounded as described above.)*
#' @srrstatsTODO {ML1.0a} *Where these terms are ultimately eschewed, these should nevertheless be used in initial documentation, along with clear explanation of, and justification for, alternative terminology.*
#' @srrstatsTODO {ML1.1} *Absent clear justification for alternative design decisions, input data should be expected to be labelled "test", "training", and, where applicable, "validation" data.*
#' @srrstatsTODO {ML1.1a} *The presence and use of these labels should be explicitly confirmed via pre-processing steps (and tested in accordance with **ML7.0**, below).*
#' @srrstatsTODO {ML1.1b} *Matches to expected labels should be case-insensitive and based on partial matching such that, for example, "Test", "test", or "testing" should all suffice.*
#' @srrstatsTODO {ML1.2} *Training and test data sets for ML software should be able to be input as a single, generally tabular, data object, with the training and test data distinguished either by* - *A specified variable containing, for example, `TRUE`/`FALSE` or `0`/`1` values, or which uses some other system such as missing (`NA`) values to denote test data); and/or* - *An additional parameter designating case or row numbers, or labels of test data.*
#' @srrstatsTODO {ML1.3} *Input data should be clearly partitioned between training and test data (for example, through having each passed as a distinct `list` item), or should enable an additional means of categorically distinguishing training from test data (such as via an additional parameter which provides explicit labels). Where applicable, distinction of validation and any other data should also accord with this standard.*
#' @srrstatsTODO {ML1.4} *Training and test data sets, along with other necessary components such as validation data sets, should be stored in their own distinctly labelled sub-directories (for distinct files), or according to an explicit and distinct labelling scheme (for example, for database connections). Labelling should in all cases adhere to **ML1.1**, above.*
#' @srrstatsTODO {ML1.5} *ML software should implement a single function which summarises the contents of test and training (and other) data sets, minimally including counts of numbers of cases, records, or files, and potentially extending to tables or summaries of file or data types, sizes, and other information (such as unique hashes for each component).*
#' @srrstatsTODO {ML1.6} *ML software which does not admit missing values, and which expects no missing values, should implement explicit pre-processing routines to identify whether data has any missing values, and should generally error appropriately and informatively when passed data with missing values. In addition, ML software which does not admit missing values should:*
#' @srrstatsTODO {ML1.6a} *Explain why missing values are not admitted.*
#' @srrstatsTODO {ML1.6b} *Provide explicit examples (in function documentation, vignettes, or both) for how missing values may be imputed, rather than simply discarded.*
#' @srrstatsTODO {ML1.7} *ML software which admits missing values should clearly document how such values are processed.*
#' @srrstatsTODO {ML1.7a} *Where missing values are imputed, software should offer multiple user-defined ways to impute missing data.*
#' @srrstatsTODO {ML1.7b} *Where missing values are imputed, the precise imputation steps should also be explicitly documented, either in tests (see **ML7.2** below), function documentation, or vignettes.*
#' @srrstatsTODO {ML1.8} *ML software should enable equal treatment of missing values for both training and test data, with optional user ability to control application to either one or both.*
#' @srrstatsTODO {ML2.0} *A dedicated function should enable pre-processing steps to be defined and parametrized.*
#' @srrstatsTODO {ML2.0a} *That function should return an object which can be directly submitted to a specified model (see section 3, below).*
#' @srrstatsTODO {ML2.0b} *Absent explicit justification otherwise, that return object should have a defined class minimally intended to implement a default `print` method which summarizes the input data set (as per **ML1.5** above) and associated transformations (see the following standard).*
#' @srrstatsTODO {ML2.1} *ML software which uses broadcasting to reconcile dimensionally incommensurate input data should offer an ability to at least optionally record transformations applied to each input file.*
#' @srrstatsTODO {ML2.2} *ML software which requires or relies upon numeric transformations of input data (such as change in mean values or variances) should allow optimal explicit specification of target values, rather than restricting transformations to default generic values only (such as transformations to z-scores).*
#' @srrstatsTODO {ML2.2a} *Where the parameters have default values, reasons for those particular defaults should be explicitly described.*
#' @srrstatsTODO {ML2.2b} *Any extended documentation (such as vignettes) which demonstrates the use of explicit values for numeric transformations should explicitly describe why particular values are used.*
#' @srrstatsTODO {ML2.3} *The values associated with all transformations should be recorded in the object returned by the function described in the preceding standard (**ML2.0**).*
#' @srrstatsTODO {ML2.4} *Default values of all transformations should be explicitly documented, both in documentation of parameters where appropriate (such as for numeric transformations), and in extended documentation such as vignettes.*
#' @srrstatsTODO {ML2.5} *ML software should provide options to bypass or otherwise switch off all default transformations.*
#' @srrstatsTODO {ML2.6} *Where transformations are implemented via distinct functions, these should be exported to a package's namespace so they can be applied in other contexts.*
#' @srrstatsTODO {ML2.7} *Where possible, documentation should be provided for how transformations may be reversed. For example, documentation may demonstrate how the values retained via **ML2.3**, above, can be used along with transformations either exported via **ML2.6** or otherwise exemplified in demonstration code to independently transform data, and then to reverse those transformations.*
#' @srrstatsTODO {ML3.0} *Model specification should be implemented as a distinct stage subsequent to specification of pre-processing routines (see Section 2, above) and prior to actual model fitting or training (see Section 4, below). In particular,*
#' @srrstatsTODO {ML3.0a} *A dedicated function should enable models to be specified without actually fitting or training them, or if this (**ML3**) and the following (**ML4**) stages are controlled by a single function, that function should have a parameter enabling models to be specified yet not fitted (for example, `nofit = FALSE`).*
#' @srrstatsTODO {ML3.0b} *That function should accept as input the objects produced by the previous Input Data Specification stage, and defined according to **ML2.0**, above.*
#' @srrstatsTODO {ML3.0c} *The function described above (**ML3.0a**) should return an object which can be directly trained as described in the following sub-section (**ML4**).*
#' @srrstatsTODO {ML3.0d} *That return object should have a defined class minimally intended to implement a default `print` method which summarises the model specification, including values of all relevant parameters.*
#' @srrstatsTODO {ML3.1} *ML software should allow the use of both untrained models, specified through model parameters only, as well as pre-trained models. Use of the latter commonly entails an ability to submit a previously-trained model object to the function defined according to **ML3.0a**, above.*
#' @srrstatsTODO {ML3.2} *ML software should enable different models to be applied to the object specifying data inputs and transformations (see sub-sections 1--2, above) without needing to re-define those preceding steps.*
#' @srrstatsTODO {ML3.3} *Where ML software implements its own distinct classes of model objects, the properties and behaviours of those specific classes of objects should be explicitly compared with objects produced by other ML software. In particular, where possible, ML software should provide extended documentation (as vignettes or equivalent) comparing model objects with those from other ML software, noting both unique abilities and restrictions of any implemented classes.*
#' @srrstatsTODO {ML3.4} *Where training rates are used, ML software should provide explicit documentation both in all functions which use training rates, and in extended form such as vignettes, of the importance of, and/or sensitivity to, different values of training rates. In particular,*
#' @srrstatsTODO {ML3.4a} *Unless explicitly justified otherwise, ML software should offer abilities to automatically determine appropriate or optimal training rates, either as distinct pre-processing stages, or as implicit stages of model training.*
#' @srrstatsTODO {ML3.4b} *ML software which provides default values for training rates should clearly document anticipated restrictions of validity of those default values; for example through clear suggestions that user-determined and -specified values may generally be necessary or preferable.*
#' @srrstatsTODO {ML3.5} *Parameters controlling optimization algorithms should minimally include:*
#' @srrstatsTODO {ML3.5a} *Specification of the type of algorithm used to explore the search space (commonly, for example, some kind of gradient descent algorithm)*
#' @srrstatsTODO {ML3.5b} *The kind of loss function used to assess distance between model estimates and desired output.*
#' @srrstatsTODO {ML3.6} *Unless explicitly justified otherwise (for example because ML software under consideration is an implementation of one specific algorithm), ML software should:*
#' @srrstatsTODO {ML3.6a} *Implement or otherwise permit usage of multiple ways of exploring search space*
#' @srrstatsTODO {ML3.6b} *Implement or otherwise permit usage of multiple loss functions.*
#' @srrstatsTODO {ML3.7} *For ML software in which algorithms are coded in C++, user-controlled use of either CPUs or GPUs (on NVIDIA processors at least) should be implemented through direct use of [`libcudacxx`](https://github.com/NVIDIA/libcudacxx).*
#' @srrstatsTODO {ML4.0} *ML software should generally implement a unified single-function interface to model training, able to receive as input a model specified according to all preceding standards. In particular, models with categorically different specifications, such as different model architectures or optimization algorithms, should be able to be submitted to the same model training function.*
#' @srrstatsTODO {ML4.1} *ML software should at least optionally retain explicit information on paths taken as an optimizer advances towards minimal loss. Such information should minimally include:*
#' @srrstatsTODO {ML4.1a} *Specification of all model-internal parameters, or equivalent hashed representation.*
#' @srrstatsTODO {ML4.1b} *The value of the loss function at each point*
#' @srrstatsTODO {ML4.1c} *Information used to advance to next point, for example quantification of local gradient.*
#' @srrstatsTODO {ML4.2} *The subsequent extraction of information retained according to the preceding standard should be explicitly documented, including through example code.*
#' @srrstatsTODO {ML4.3} *All parameters controlling batch processing and associated terminology should be explicitly documented, and it should not, for example, be presumed that users will understand the definition of "epoch" as implemented in any particular ML software.*
#' @srrstatsTODO {ML4.4} *Explicit guidance should be provided on selection of appropriate values for parameter controlling batch processing, for example, on trade-offs between batch sizes and numbers of epochs (with both terms provided as Control Parameters in accordance with the preceding standard, **ML3**).*
#' @srrstatsTODO {ML4.5} *ML software may optionally include a function to estimate likely time to train a specified model, through estimating initial timings from a small sample of the full batch.*
#' @srrstatsTODO {ML4.6} *ML software should by default provide explicit information on the progress of batch jobs (even where those jobs may be implemented in parallel on GPUs). That information may be optionally suppressed through additional parameters.*
#' @srrstatsTODO {ML4.7} *ML software should provide an ability to combine results from multiple re-sampling iterations using a single parameter specifying numbers of iterations.*
#' @srrstatsTODO {ML4.8} *Absent any additional specification, re-sampling algorithms should by default partition data according to proportions of original test and training data.*
#' @srrstatsTODO {ML4.8a} *Re-sampling routines of ML software should nevertheless offer an ability to explicitly control or override such default proportions of test and training data.*
#' @srrstatsTODO {ML5.0} *The result of applying the training processes described above should be contained within a single model object returned by the function defined according to **ML4.0**, above. Even where the output reflects application to a test data set, the resultant object need not include any information on model performance (see **ML5.3**--**ML5.4**, below).*
#' @srrstatsTODO {ML5.0a} *That object should either have its own class, or extend some previously-defined class.*
#' @srrstatsTODO {ML5.0b} *That class should have a defined `print` method which summarises important aspects of the model object, including but not limited to summaries of input data and algorithmic control parameters.*
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
#' @noRd
NULL




