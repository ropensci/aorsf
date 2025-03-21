

#' Get predictors from a formula
#'
#' @param formula a `formula` object
#' @param outcomes a character vector containing names of outcome
#'  variable(s) in `formula`
#' @param data the data where `formula` will be evaluated.
#'
#' @return a character vector of predictors in `formula`
#' @noRd
#'
#' @examples
#'
#' f <- a ~ b + c
#' get_predictors(f, outcomes = 'a')
#'
#' data_example <- matrix(0, ncol = length(letters))
#' colnames(data_example) <- letters
#' data_example <- as.data.frame(data_example)
#' f <- a + b ~ . - z
#'
#' get_predictors(f, outcomes = c("a", "b"), data = data_example)
#'
get_predictors <- function(formula,
                           outcomes = NULL,
                           data = NULL) {

 stopifnot(inherits(formula, "formula"))

 if(!is.null(data)){
  stopifnot(is.data.frame(data))
 }

 if(!is.null(outcomes)){
  stopifnot(is.character(outcomes))
 } else {
  outcomes <- all.vars(formula[[2]])
 }

 # parse to a character valued vector
 formula_chr <- as.character(formula)
 # Right-hand side
 rhs <- formula_chr[length(formula_chr)]

 # Split RHS by '+' and '-' to get variable names and operations
 # Add spaces around operators to help splitting
 rhs_spaced <- gsub("([+-])", " \\1 ", rhs)
 tokens <- unlist(strsplit(rhs_spaced, "\\s+"))
 tokens <- tokens[tokens != ""]  # Remove empty elements

 # Handle '.' if present
 if ("." %in% tokens) {

  if (is.null(data))
   stop("Formula contains '.' but 'data' was not provided.",
        call. = FALSE)

  # Exclude response variable
  dot_vars <- setdiff(names(data), outcomes)

  # Replace '.' with actual variable names
  dot_index <- which(tokens == ".")
  tokens <- append(tokens, dot_vars, after = dot_index)

  # Remove the '.'
  tokens <- tokens[-dot_index]

 }

 drop_vars <- tokens[which(tokens == '-') + 1]

 unique(setdiff(tokens, c("+", "-", drop_vars)))

}
