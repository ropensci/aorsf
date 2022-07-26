#' Mayo Clinic Primary Biliary Cholangitis Data
#'
#' These data are a light modification of the [survival::pbc] data. The
#'  modifications are:
#'
#'  1. removed rows with missing data
#'
#'  2. converted `status` into 0 for censor or transplant, 1 for dead
#'
#'  3. converted `stage` into an ordered factor.
#'
#'  4. converted `trt`, `ascites`, `hepato`, `spiders`, and `edema` into factors.
#'
#'
#'
#' @format A data frame with 276 rows and 20 variables:
#' \describe{
#'   \item{id}{case number}
#'   \item{time}{number of days between registration and the earlier of death, transplantion, or study analysis in July, 1986}
#'   \item{status}{ status at endpoint, 0 for censored or transplant, 1 for dead}
#'   \item{trt}{randomized treatment group: D-penicillmain or placebo}
#'   \item{age}{in years}
#'   \item{sex}{m/f}
#'   \item{ascites}{presence of ascites}
#'   \item{hepato}{presence of hepatomegaly or enlarged liver}
#'   \item{spiders}{blood vessel malformations in the skin}
#'   \item{edema}{0 no edema, 0.5 untreated or successfully treated, 1 edema despite diuretic therapy}
#'   \item{bili}{serum bilirubin (mg/dl)}
#'   \item{chol}{serum cholesterol (mg/dl)}
#'   \item{albumin}{serum albumin (g/dl)}
#'   \item{copper}{urine copper (ug/day)}
#'   \item{alk.phos}{alkaline phosphotase (U/liter)}
#'   \item{ast}{aspartate aminotransferase, once called SGOT (U/ml)}
#'   \item{trig}{triglycerides (mg/dl)}
#'   \item{platelet}{platelet count}
#'   \item{protime}{standardized blood clotting time}
#'   \item{stage}{histologic stage of disease (needs biopsy)}
#' }
#' @source T Therneau and P Grambsch (2000), Modeling Survival Data: Extending the Cox Model, Springer-Verlag, New York. ISBN: 0-387-98784-3.
"pbc_orsf"



