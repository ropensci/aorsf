
# Forest class ----

Forest <- R6::R6Class(
 "Forest",
 cloneable = FALSE,


 public = list( # public ----

  # user-facing fields
  data = NULL,
  trained = FALSE,
  forest = list(),
  importance = vector(mode = 'double', length = 0),
  pred_oobag = matrix(ncol = 0, nrow = 0),
  eval_oobag = list(),

  # fitting fields
  formula = NULL,
  control = NULL,
  weights = NULL,

  # count fields
  n_tree = NULL,
  n_split = NULL,
  n_retry = NULL,
  n_thread = NULL,
  n_obs = NULL,
  mtry = NULL,

  # sampling fields
  sample_with_replacement = NULL,
  sample_fraction = NULL,

  # tree fields
  tree_seeds = NULL,
  tree_type = NULL,
  leaf_min_events = NULL,
  leaf_min_obs = NULL,
  split_rule = NULL,
  split_min_events = NULL,
  split_min_obs = NULL,
  split_min_stat = NULL,

  # out-of-bootstrap aggregate (oobag) fields
  oobag_pred_mode = NULL,
  oobag_pred_type = NULL,
  pred_horizon = NULL,
  oobag_eval_every = NULL,
  oobag_eval_function = NULL,

  # variable importance fields
  importance_type = NULL,
  importance_max_pvalue = NULL,
  importance_group_factors = NULL,

  # missing data fields
  na_action = NULL,

  # miscellaneous
  verbose_progress = NULL,

  # initialization
  initialize = function(data,
                        formula,
                        control,
                        weights = NULL,
                        n_tree,
                        n_split,
                        n_retry,
                        n_thread,
                        mtry = NULL,
                        sample_with_replacement,
                        sample_fraction,
                        leaf_min_events,
                        leaf_min_obs,
                        split_rule,
                        split_min_events,
                        split_min_obs,
                        split_min_stat,
                        oobag_pred_type,
                        oobag_pred_horizon,
                        oobag_eval_every,
                        oobag_fun = NULL,
                        importance_type,
                        importance_max_pvalue,
                        importance_group_factors,
                        tree_seeds,
                        na_action,
                        verbose_progress) {

   self$data     <- data
   self$formula  <- formula
   self$control  <- control
   self$weights  <- weights
   self$n_tree   <- n_tree
   self$n_split  <- n_split
   self$n_retry  <- n_retry
   self$n_thread <- n_thread
   self$mtry     <- mtry
   self$sample_with_replacement  <- sample_with_replacement
   self$sample_fraction          <- sample_fraction
   self$leaf_min_events          <- leaf_min_events
   self$leaf_min_obs             <- leaf_min_obs
   self$split_rule               <- split_rule
   self$split_min_events         <- split_min_events
   self$split_min_obs            <- split_min_obs
   self$split_min_stat           <- split_min_stat
   self$oobag_pred_type          <- oobag_pred_type
   self$pred_horizon             <- oobag_pred_horizon
   self$oobag_eval_every         <- oobag_eval_every
   self$oobag_eval_function      <- oobag_fun
   self$importance_type          <- importance_type
   self$importance_max_pvalue    <- importance_max_pvalue
   self$importance_group_factors <- importance_group_factors
   self$tree_seeds               <- tree_seeds
   self$na_action                <- na_action
   self$verbose_progress         <- verbose_progress

   private$init()

  },

  print = function(){

   info_model_type <- switch(self$tree_type,
                             'survival' = "Cox regression",
                             'regression' = "Linear regression",
                             'classification' = "Logistic regression")

   info_lincomb_type <- switch(self$control$lincomb_type,
                               'glm'    = NULL,
                               'net'    = "Penalized ",
                               'custom' = "Custom user function")


   if(self$control$lincomb_type != 'custom'){

    if(self$control$lincomb_iter_max == 1){
     info_lincomb_type <- "Accelerated "
    }

    info_lincomb_type <- paste0(info_lincomb_type, info_model_type)

   }

   if(!is_empty(self$eval_oobag)){

    oobag_type <- self$eval_oobag$stat_type

    oobag_stat <- "Not estimated"

    if(self$trained){
     oobag_stat <- last_value(self$eval_oobag$stat_values)
     oobag_stat <- round_magnitude(oobag_stat)
    }

   } else {

    oobag_stat <- "TBD"
    oobag_type <- "TBD"

   }


   forest_text <- paste("Oblique random", self$tree_type, "forest")

   header <- paste0('---------- ', forest_text, '\n')

   if(!self$trained){
    header <- paste0('Untrained ', tolower(forest_text), '\n')
   }


   n_predictors <- length(private$data_names$x_original)

   cat(header,
       paste0('     Linear combinations: ', info_lincomb_type             ),
       paste0('          N observations: ', self$n_obs                    ),

       if(self$tree_type == 'survival')
        paste0('                N events: ', self$n_events                ),

       paste0('                 N trees: ', self$n_tree                   ),
       paste0('      N predictors total: ', n_predictors                  ),
       paste0('   N predictors per node: ', self$mtry                     ),
       paste0(' Average leaves per tree: ', self$n_leaves_mean %||% "TBD" ),
       paste0('Min observations in leaf: ', self$leaf_min_obs             ),
       paste0('      Min events in leaf: ', self$leaf_min_events          ),
       paste0('          OOB stat value: ', oobag_stat                    ),
       paste0('           OOB stat type: ', oobag_type                    ),
       paste0('     Variable importance: ', self$importance_type          ),
       '\n-----------------------------------------',
       sep = '\n')

   invisible(self)

  },

  # getters

  get_var_type = function(.name){
   class(self$data[[.name]])[1]
  }


 ),

 private = list( # private ----

  data_rows_complete = NULL,
  data_types = NULL,
  data_names = NULL,
  data_fctrs = NULL,
  data_units = NULL,
  data_means = NULL,
  data_stdev = NULL,
  data_modes = NULL,
  data_bounds = NULL,

  # runs checks and sets defaults where needed
  init = function() {

   private$check_data()
   private$check_formula()

   private$init_data()
   private$init_mtry()
   private$init_weights()

   private$check_control()
   private$check_n_tree()
   private$check_n_split()
   private$check_n_retry()
   private$check_n_thread()
   private$check_mtry()
   private$check_sample_with_replacement()
   private$check_sample_fraction()
   private$check_leaf_min_obs()
   private$check_split_rule()
   private$check_split_min_obs()
   private$check_split_min_stat()
   private$check_oobag_pred_type()
   private$check_oobag_eval_every()
   private$check_importance_type()
   private$check_importance_max_pvalue()
   private$check_importance_group_factors()
   private$check_tree_seeds()
   private$check_na_action()


   private$init_oobag_eval_function()

   private$init_oobag_pred_mode()

   private$init_tree_seeds()

   private$init_internal()


  },
  init_internal = function(){
   NULL
  },
  init_data = function(){

   formula_terms <- stats::terms(x = self$formula, data = self$data)

   if(attr(formula_terms, 'response') == 0)
    stop("formula should have a response", call. = FALSE)

   names_x_data <- attr(formula_terms, 'term.labels')
   names_y_data <- all.vars(self$formula[[2]])

   # check factors in x data
   fctr_check(self$data, names_x_data)
   fctr_id_check(self$data, names_x_data)

   private$data_rows_complete <-
    stats::complete.cases(select_cols(self$data, names_x_data))

   self$n_obs <- ifelse(test = self$na_action == 'omit',
                        yes  = length(private$data_rows_complete),
                        no   = nrow(self$data))



   private$check_var_names(c(names_x_data, names_y_data))

   private$data_names <- list(y = names_y_data,
                              x_original = names_x_data)

   private$check_var_missing()

   types_y_data <- vapply(names_y_data, self$get_var_type, character(1))
   types_x_data <- vapply(names_x_data, self$get_var_type, character(1))

   valid_y_types <- c('numeric', 'integer', 'units', 'factor')
   valid_x_types <- c(valid_y_types, 'ordered')


   private$check_var_types(types_y_data, valid_types = valid_y_types)
   private$check_var_types(types_x_data, valid_types = valid_x_types)

   private$data_types <- list(y = types_y_data, x = types_x_data)
   private$data_fctrs <- fctr_info(self$data, private$data_names$x_original)

   private$init_ref_code_names()

   unit_names <- c(names_y_data[types_y_data == 'units'],
                   names_x_data[types_x_data == 'units'])

   private$data_units <- unit_info(data = self$data, .names = unit_names)

  },
  init_tree_seeds = function(){
   if(is.null(self$tree_seeds)){
    self$tree_seeds <- sample(1e6, size = 1)
   }

   if(length(self$tree_seeds) == 1){
    if(self$n_tree > 1) set.seed(self$tree_seeds)
    self$tree_seeds <- sample(self$n_tree*10, size = self$n_tree)
   }
  },
  init_ref_code_names = function(){

   xref <- xnames <- private$data_names$x_original
   fi <- private$data_fctrs

   for (i in seq_along(fi$cols)){

    if(fi$cols[i] %in% xnames){

     if(!fi$ordr[i]){

      xref <- insert_vals(
       vec = xref,
       where = xref %==% fi$cols[i],
       what = fi$keys[[i]][-1]
      )

     }
    }

   }

   private$data_names$x_ref_code = xref

  },
  init_oobag_pred_mode = function(){

   self$oobag_pred_mode <- self$oobag_pred_type != "none"

   if(self$oobag_pred_mode && self$sample_fraction == 1){
    stop(
     "cannot compute out-of-bag predictions if no samples are out-of-bag.",
     " Try setting sample_fraction < 1 or oobag_pred_type = 'none'.",
     call. = FALSE
    )
   }

  },
  init_mtry = function(){

   n_col_x <- length(private$data_names$x_ref_code)

   if(is.null(self$mtry)){

    self$mtry <-  ceiling(sqrt(n_col_x))

   } else {

    check_arg_lteq(
     arg_value = self$mtry,
     arg_name = 'mtry',
     bound = n_col_x,
     append_to_msg = "(number of columns in the one-hot encoded x-matrix)"
    )

   }

   if(is.null(self$control$lincomb_df_target)){

    self$control$lincomb_df_target <- self$mtry

   } else {

    check_arg_lteq(
     arg_value = self$control$lincomb_df_target,
     arg_name = 'df_target',
     bound = self$mtry,
     append_to_msg = "(number of randomly selected predictors)"
    )

   }

  },
  init_weights = function(){

   # set weights as 1 if user did not supply them.
   # length of weights depends on how missing are handled.
   if(is.null(self$weights)){

    self$weights <- rep(1, self$n_obs)

   } else {

    private$check_weights()

   }

  },
  init_oobag_eval_function = function(){

   if(is.null(self$oobag_eval_function)){

    self$oobag_eval_function <- function(x) x

   } else {

    private$check_oobag_eval_function()

   }

  },

  # checkers

  check_data = function(){

   check_arg_is(arg_value = self$data,
                arg_name = 'data',
                expected_class = 'data.frame')

   # Minimal checks for now, other checks occur later
   if(nrow(self$data) == 0 || ncol(self$data) ==  0){
    stop("training data are empty",
         call. = FALSE)
   }

   # check for blanks first b/c the check for non-standard symbols
   # will detect blanks with >1 empty characters

   blank_names <- grepl(pattern = '^\\s*$', x = names(self$data))

   if(any(blank_names)){

    s_if_plural_blank_otherwise <- ""

    to_list <- which(blank_names)

    if(length(to_list) > 1) s_if_plural_blank_otherwise <- "s"

    last <- ifelse(length(to_list) == 2, ' and ', ', and ')

    stop("Blank or empty names detected in training data: see column",
         s_if_plural_blank_otherwise, " ",
         paste_collapse(x = to_list, last = last),
         call. = FALSE)

   }

   ns_names <- grepl(pattern = '[^a-zA-Z0-9\\.\\_]+', x = names(self$data))

   if(any(ns_names)){

    last <- ifelse(sum(ns_names) == 2, ' and ', ', and ')

    stop("Non-standard names detected in training data: ",
         paste_collapse(x = names(self$data)[ns_names],
                        last = last),
         call. = FALSE)

   }

  },
  check_var_names = function(.names){

   names_not_found <- setdiff(c(.names), names(self$data))

   if(!is_empty(names_not_found)){
    msg <- paste0(
     "variables in formula were not found in data: ",
     paste_collapse(names_not_found, last = ' and ')
    )
    stop(msg, call. = FALSE)
   }

  },
  check_var_types = function(var_types, valid_types){

   good_vars <- var_types %in% valid_types

   if(!all(good_vars)){

    bad_vars <- which(!good_vars)

    vars_to_list <- .names[bad_vars]
    types_to_list <- var_types[bad_vars]

    meat <- paste0(' <', vars_to_list, '> has type <',
                   types_to_list, '>', collapse = '\n')

    msg <- paste0("some variables have unsupported type:\n",
                  meat, '\nsupported types are ',
                  paste_collapse(valid_types, last = ' and '))

    stop(msg, call. = FALSE)

   }

   var_types


  },
  check_var_missing = function(){

   if(self$na_action == 'fail'){

    if(any(is.na(select_cols(self$data, private$data_names$x_original)))){
     stop("Please remove missing values from data, or impute them.",
          call. = FALSE)
    }

   }

   if(any(is.na(select_cols(self$data, private$data_names$y))))
    stop("Please remove missing values from the outcome variable(s).",
         call. = FALSE)

   for(i in private$data_names$x_original){

    if(collapse::allNA(self$data[[i]])){
     stop("column ", i, " has no observed values.",
          call. = FALSE)
    }

    if(any(is.infinite(self$data[[i]]))){
     stop("Please remove infinite values from ", i, ".",
          call. = FALSE)
    }

   }

  },
  check_formula = function(){

   check_arg_is(arg_value = self$formula,
                arg_name = 'formula',
                expected_class = 'formula')

   if(length(self$formula) != 3){
    stop("formula must be two sided, i.e. left side ~ right side",
         call. = FALSE)
   }

   formula_deparsed <- as.character(self$formula)[[3]]

   for( symbol in c("*", "^", ":", "(", ")", "["," ]", "|", "%") ){

    if(grepl(symbol, formula_deparsed, fixed = TRUE)){

     stop("unrecognized symbol in formula: ", symbol,
          "\norsf recognizes '+', '-', and '.' symbols.",
          call. = FALSE)

    }

   }

  },
  check_control = function(){

   check_arg_is(arg_value = self$control,
                arg_name = 'control',
                expected_class = 'orsf_control')

   if(self$control$lincomb_type == 'net'){

    if (!requireNamespace("glmnet", quietly = TRUE)) {
     stop(
      "Package \"glmnet\" must be installed to use",
      " orsf_control_net() with orsf().",
      call. = FALSE
     )
    }

   }

  },
  check_weights = function(){

   check_arg_type(arg_value = self$weights,
                  arg_name = 'weights',
                  expected_type = 'numeric')

   check_arg_gteq(arg_value = self$weights,
                  arg_name = 'weights',
                  bound = 0)

   check_arg_length(
    arg_value = self$weights,
    arg_name  = 'weights',
    expected_length = self$n_obs
   )

  },
  check_n_tree = function(){

   check_arg_type(arg_value = self$n_tree,
                  arg_name = 'n_tree',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$n_tree,
                        arg_name = 'n_tree')

   check_arg_gteq(arg_value = self$n_tree,
                  arg_name = 'n_tree',
                  bound = 1)

   check_arg_length(arg_value = self$n_tree,
                    arg_name = 'n_tree',
                    expected_length = 1)

  },
  check_n_split = function(){

   check_arg_type(arg_value = self$n_split,
                  arg_name = 'n_split',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$n_split,
                        arg_name = 'n_split')

   check_arg_gteq(arg_value = self$n_split,
                  arg_name = 'n_split',
                  bound = 1)

   check_arg_length(arg_value = self$n_split,
                    arg_name = 'n_split',
                    expected_length = 1)

  },
  check_n_retry = function(){

   check_arg_type(arg_value = self$n_retry,
                  arg_name = 'n_retry',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$n_retry,
                        arg_name = 'n_retry')

   check_arg_gteq(arg_value = self$n_retry,
                  arg_name = 'n_retry',
                  bound = 0)

   check_arg_length(arg_value = self$n_retry,
                    arg_name = 'n_retry',
                    expected_length = 1)

  },
  check_n_thread = function(){

   check_arg_type(arg_value = self$n_thread,
                  arg_name = 'n_thread',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$n_thread,
                        arg_name = 'n_thread')

   check_arg_gteq(arg_value = self$n_thread,
                  arg_name = 'n_thread',
                  bound = 0)

   check_arg_length(arg_value = self$n_thread,
                    arg_name = 'n_thread',
                    expected_length = 1)

  },
  check_mtry = function(){

   if(!is.null(self$mtry)){

    check_arg_type(arg_value = self$mtry,
                   arg_name = 'mtry',
                   expected_type = 'numeric')

    check_arg_is_integer(arg_value = self$mtry,
                         arg_name = 'mtry')

    check_arg_gteq(arg_value = self$mtry,
                   arg_name = 'mtry',
                   bound = 1)

    check_arg_length(arg_value = self$mtry,
                     arg_name = 'mtry',
                     expected_length = 1)

   }

  },
  check_sample_with_replacement = function(){

   check_arg_type(arg_value = self$sample_with_replacement,
                  arg_name = 'sample_with_replacement',
                  expected_type = 'logical')

   check_arg_length(arg_value = self$sample_with_replacement,
                    arg_name = 'sample_with_replacement',
                    expected_length = 1)

  },
  check_sample_fraction = function(){

   check_arg_type(arg_value = self$sample_fraction,
                  arg_name = 'sample_fraction',
                  expected_type = 'numeric')

   check_arg_gt(arg_value = self$sample_fraction,
                arg_name = 'sample_fraction',
                bound = 0)

   check_arg_lteq(arg_value = self$sample_fraction,
                  arg_name = 'sample_fraction',
                  bound = 1)

   check_arg_length(arg_value = self$sample_fraction,
                    arg_name = 'sample_fraction',
                    expected_length = 1)

  },
  check_leaf_min_obs = function(){

   check_arg_type(arg_value = self$leaf_min_obs,
                  arg_name = 'leaf_min_obs',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$leaf_min_obs,
                        arg_name = 'leaf_min_obs')

   check_arg_gteq(arg_value = self$leaf_min_obs,
                  arg_name = 'leaf_min_obs',
                  bound = 1)

   check_arg_length(arg_value = self$leaf_min_obs,
                    arg_name = 'leaf_min_obs',
                    expected_length = 1)

   check_arg_lteq(arg_value = self$leaf_min_obs,
                  arg_name = "leaf_min_obs",
                  bound = round(self$n_obs / 2),
                  append_to_msg = "(number of observations divided by 2)")

  },
  check_split_rule = function(){

   check_arg_type(arg_value = self$split_rule,
                  arg_name = 'split_rule',
                  expected_type = 'character')

   check_arg_length(arg_value = self$split_rule,
                    arg_name = 'split_rule',
                    expected_length = 1)

  },
  check_split_rule_internal = function(){

   NULL

  },
  check_split_min_obs = function(){

   check_arg_type(arg_value = self$split_min_obs,
                  arg_name = 'split_min_obs',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$split_min_obs,
                        arg_name = 'split_min_obs')

   check_arg_gteq(arg_value = self$split_min_obs,
                  arg_name = 'split_min_obs',
                  bound = 1)

   check_arg_length(arg_value = self$split_min_obs,
                    arg_name = 'split_min_obs',
                    expected_length = 1)

   check_arg_lt(arg_value = self$split_min_obs,
                arg_name = "split_min_obs",
                bound = self$n_obs,
                append_to_msg = "(number of observations)")

  },
  check_split_min_stat = function(){

   check_arg_type(arg_value = self$split_min_stat,
                  arg_name = 'split_min_stat',
                  expected_type = 'numeric')

   check_arg_gteq(arg_value = self$split_min_stat,
                  arg_name = 'split_min_stat',
                  bound = 0)

   if(self$split_rule %in% c('cstat', 'gini')){

    check_arg_lt(arg_value = self$split_min_stat,
                 arg_name = 'split_min_stat',
                 bound = 1,
                 append_to_msg = paste0("(split stat <",
                                        self$split_rule,
                                        "> is always < 1)"))

   }

   check_arg_length(arg_value = self$split_min_stat,
                    arg_name = 'split_min_stat',
                    expected_length = 1)

  },
  check_oobag_pred_type = function(){

   check_arg_type(arg_value = self$oobag_pred_type,
                  arg_name = 'oobag_pred_type',
                  expected_type = 'character')

   check_arg_length(arg_value = self$oobag_pred_type,
                    arg_name = 'oobag_pred_type',
                    expected_length = 1)

  },
  check_oobag_pred_type_internal = function(){
   NULL
  },
  check_oobag_eval_every = function(){

   check_arg_type(arg_value = self$oobag_eval_every,
                  arg_name = 'oobag_eval_every',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$oobag_eval_every,
                        arg_name = 'oobag_eval_every')

   check_arg_gteq(arg_value = self$oobag_eval_every,
                  arg_name = 'oobag_eval_every',
                  bound = 1)

   check_arg_lteq(arg_value = self$oobag_eval_every,
                  arg_name = 'oobag_eval_every',
                  bound = self$n_tree)

   check_arg_length(arg_value = self$oobag_eval_every,
                    arg_name = 'oobag_eval_every',
                    expected_length = 1)

  },
  check_importance_type = function(){

   check_arg_type(arg_value = self$importance_type,
                  arg_name = 'importance',
                  expected_type = 'character')

   check_arg_length(arg_value = self$importance_type,
                    arg_name = 'importance',
                    expected_length = 1)

   check_arg_is_valid(arg_value = self$importance_type,
                      arg_name = 'importance',
                      valid_options = c("none",
                                        "anova",
                                        "negate",
                                        "permute"))

  },
  check_importance_max_pvalue = function(){

   check_arg_type(arg_value = self$importance_max_pvalue,
                  arg_name = 'importance_max_pvalue',
                  expected_type = 'numeric')

   check_arg_gt(arg_value = self$importance_max_pvalue,
                arg_name = 'importance_max_pvalue',
                bound = 0)

   check_arg_lt(arg_value = self$importance_max_pvalue,
                arg_name = 'importance_max_pvalue',
                bound = 1)

   check_arg_length(arg_value = self$importance_max_pvalue,
                    arg_name = 'importance_max_pvalue',
                    expected_length = 1)


  },
  check_importance_group_factors = function(){

   check_arg_type(arg_value = self$importance_group_factors,
                  arg_name = 'group_factors',
                  expected_type = 'logical')

   check_arg_length(arg_value = self$importance_group_factors,
                    arg_name = 'group_factors',
                    expected_length = 1)

  },
  check_tree_seeds = function(){

   if(!is.null(self$tree_seeds)){
    check_arg_type(arg_value = self$tree_seeds,
                   arg_name = 'tree_seed',
                   expected_type = 'numeric')

    check_arg_is_integer(arg_value = self$tree_seeds,
                         arg_name = 'tree_seeds')

    if(length(self$tree_seeds) > 1 &&
       length(self$tree_seeds) != self$n_tree){

     stop('tree_seeds should have length = 1 or length = ",
          "n_tree <', self$n_tree,
          "> (the number of trees) but instead has length <",
          length(self$tree_seeds), ">", call. = FALSE)

    }
   }

  },
  check_na_action = function(){

   check_arg_type(arg_value = self$na_action,
                  arg_name = 'na_action',
                  expected_type = 'character')

   check_arg_length(arg_value = self$na_action,
                    arg_name = 'na_action',
                    expected_length = 1)

   check_arg_is_valid(arg_value = self$na_action,
                      arg_name = 'na_action',
                      valid_options = c("fail", "omit", "impute_meanmode"))

  },
  check_oobag_eval_function = function(){

   if(!is.null(self$oobag_eval_function)){

    oobag_fun_args <- names(formals(self$oobag_eval_function))

    if(length(oobag_fun_args) != 3) stop(
     "oobag_fun should have 3 input arguments but instead has ",
     length(oobag_fun_args),
     call. = FALSE
    )

    if(oobag_fun_args[1] != 'y_mat') stop(
     "the first input argument of oobag_fun should be named 'y_mat' ",
     "but is instead named '", oobag_fun_args[1], "'",
     call. = FALSE
    )

    if(oobag_fun_args[2] != 'w_vec') stop(
     "the second input argument of oobag_fun should be named 'w_vec' ",
     "but is instead named '", oobag_fun_args[1], "'",
     call. = FALSE
    )

    if(oobag_fun_args[3] != 's_vec') stop(
     "the third input argument of oobag_fun should be named 's_vec' ",
     "but is instead named '", oobag_fun_args[2], "'",
     call. = FALSE
    )

    test_time <- seq(from = 1, to = 5, length.out = 100)
    test_status <- rep(c(0,1), each = 50)

    .y_mat <- cbind(time = test_time, status = test_status)
    .w_vec <- rep(1, times = 100)
    .s_vec <- seq(0.9, 0.1, length.out = 100)

    test_output <- try(self$oobag_eval_function(y_mat = .y_mat,
                                                w_vec = .w_vec,
                                                s_vec = .s_vec),
                       silent = FALSE)

    if(is_error(test_output)){

     stop("oobag_fun encountered an error when it was tested. ",
          "Please make sure your oobag_fun works for this case:\n\n",
          "test_time <- seq(from = 1, to = 5, length.out = 100)\n",
          "test_status <- rep(c(0,1), each = 50)\n\n",
          "y_mat <- cbind(time = test_time, status = test_status)\n",
          "w_vec <- rep(1, times = 100)\n",
          "s_vec <- seq(0.9, 0.1, length.out = 100)\n\n",
          "test_output <- oobag_fun(y_mat = y_mat, w_vec = w_vec, s_vec = s_vec)\n\n",
          "test_output should be a numeric value of length 1",
          call. = FALSE)

    }

    if(!is.numeric(test_output)) stop(
     "oobag_fun should return a numeric output but instead returns ",
     "output of type ", class(test_output)[1],
     call. = FALSE
    )

    if(length(test_output) != 1) stop(
     "oobag_fun should return output of length 1 instead returns ",
     "output of length ", length(test_output),
     call. = FALSE
    )

   }

  },
  check_verbose_progress = function(){

   check_arg_type(arg_value = self$verbose_progress,
                  arg_name = 'verbose_progress',
                  expected_type = 'logical')

   check_arg_length(arg_value = self$verbose_progress,
                    arg_name = 'verbose_progress',
                    expected_length = 1)

  }

 )

)

# ForestSurvival class ----

ForestSurvival <- R6::R6Class(
 "ForestSurvival",
 inherit = Forest,
 cloneable = FALSE,
 public = list(

  leaf_min_events = NULL,
  split_min_events = NULL,
  pred_horizon = NULL

 ),

 private = list(

  row_sort = NULL,
  max_time = NULL,
  n_events = NULL,

  check_split_rule_internal= function(){

   check_arg_is_valid(arg_value = self$split_rule,
                      arg_name = 'split_rule',
                      valid_options = c("logrank", "cstat"))

  },
  check_oobag_pred_type_internal = function(){


   check_arg_is_valid(arg_value = self$oobag_pred_type,
                      arg_name = 'oobag_pred_type',
                      valid_options = c("none", "surv", "risk",
                                        "chf", "mort", "leaf"))

  },
  check_pred_horizon = function(){

   if(self$oobag_pred_mode)
    arg_name <- 'oobag_pred_horizon'
   else
    arg_name <- 'pred_horizon'

   check_arg_type(arg_value = self$pred_horizon,
                  arg_name = 'pred_horizon',
                  expected_type = 'numeric')

   for(i in seq_along(pred_horizon)){

    check_arg_gteq(arg_value = self$pred_horizon[i],
                   arg_name = 'pred_horizon',
                   bound = 0)

   }

  },
  check_leaf_min_events = function(){

    check_arg_type(arg_value = self$leaf_min_events,
                   arg_name = 'leaf_min_events',
                   expected_type = 'numeric')

    check_arg_is_integer(arg_value = self$leaf_min_events,
                         arg_name = 'leaf_min_events')

    check_arg_gteq(arg_value = self$leaf_min_events,
                   arg_name = 'leaf_min_events',
                   bound = 1)

    check_arg_length(arg_value = self$leaf_min_events,
                     arg_name = 'leaf_min_events',
                     expected_length = 1)

    check_arg_lteq(
     arg_value = self$leaf_min_events,
     arg_name = 'leaf_min_events',
     bound = round(private$n_events / 2),
     append_to_msg = "(number of events divided by 2)"
    )

  },
  check_split_min_events = function(){

   check_arg_type(arg_value = self$split_min_events,
                  arg_name = 'split_min_events',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = self$split_min_events,
                        arg_name = 'split_min_events')

   check_arg_gteq(arg_value = self$split_min_events,
                  arg_name = 'split_min_events',
                  bound = 1)

   check_arg_length(arg_value = self$split_min_events,
                    arg_name = 'split_min_events',
                    expected_length = 1)

   check_arg_lt(
    arg_value = self$split_min_events,
    arg_name = "split_min_events",
    bound = private$n_events,
    append_to_msg = "(number of events)"
   )

  },

  init_internal = function(){

   self$tree_type <- "survival"

   y <- self$data[, private$data_names$y]

   # if pred_horizon is unspecified, provide sensible default
   # if it is specified, check for correctness
   if(is.null(self$pred_horizon)){
    self$pred_horizon <- collapse::fmedian(y[, 1])
   } else {
    private$check_pred_horizon()
   }

   # order observations by event time and event status
   private$row_sort <- collapse::radixorder(y[, 1], -y[, 2])
   private$max_time <- y[last_value(private$row_sort), 1]
   private$n_events <- collapse::fsum(y[, 2])

   private$check_split_rule_internal()
   private$check_oobag_pred_type_internal()
   private$check_leaf_min_events()
   private$check_split_min_events()



  }

 )

)

# ForestClassification class ----

ForestClassification <- R6::R6Class(
 "ForestClassification",
 inherit = Forest,
 cloneable = FALSE,
 public = list(),
 private = list(

  check_split_rule_internal = function(){

   check_arg_is_valid(arg_value = self$split_rule,
                      arg_name = 'split_rule',
                      valid_options = c("gini", "cstat"))

  },
  check_oobag_pred_type_internal = function(){


   check_arg_is_valid(arg_value = self$oobag_pred_type,
                      arg_name = 'oobag_pred_type',
                      valid_options = c("none", "prob", "class", "leaf"))

  },
  init_internal = function(){

   self$tree_type <- "Classification"
   check_split_rule_internal()
   check_oobag_pred_type_internal()

  }

 )
)
