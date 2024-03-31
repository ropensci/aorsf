
# TODO:
# - add nocov to cpp
# - tests for survival forest w/no censored


# ObliqueForest class ----

ObliqueForest <- R6::R6Class(
 "ObliqueForest",


 public = list(

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

  # this is used if length(tree_seeds) = 1
  forest_seed = NULL,

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
  pred_type = NULL,
  pred_horizon = NULL,
  oobag_eval_type = NULL,
  oobag_eval_every = NULL,
  oobag_eval_function = NULL,

  # prediction fields
  pred_aggregate = NULL,

  # variable importance fields
  importance_type = NULL,
  importance_max_pvalue = NULL,
  importance_group_factors = NULL,

  # missing data fields
  na_action = NULL,

  # miscellaneous
  verbose_progress = NULL,


  # Assign incoming member variables, then run input checks.
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
                        pred_type,
                        oobag_pred_horizon,
                        oobag_eval_every = NULL,
                        oobag_fun = NULL,
                        importance_type,
                        importance_max_pvalue,
                        importance_group_factors,
                        tree_seeds,
                        na_action,
                        verbose_progress) {

   # if these arguments were specified in the initial
   # call or in updates, then they should be checked.
   # Note the arguments in this list are NULL in orsf()
   # by default, which means their defaults are dynamic.
   # other arguments with non-dynamic defaults should
   # always be checked if a user wants to use update().

   private$user_specified <- list(
    control             = !is.null(control),
    lincomb_df_target   = !is.null(control$lincomb_df_target),
    weights             = !is.null(weights),
    mtry                = !is.null(mtry),
    split_rule          = !is.null(split_rule),
    split_min_stat      = !is.null(split_min_stat),
    pred_type           = !is.null(pred_type),
    pred_horizon        = !is.null(oobag_pred_horizon),
    oobag_eval_function = !is.null(oobag_fun),
    tree_seeds          = !is.null(tree_seeds),
    oobag_eval_every    = !is.null(oobag_eval_every)
   )

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
   self$pred_type                <- pred_type
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

  # Update: re-initialize if dynamic args were unspecified
  update = function(args) {

   # for args with default values of NULL, keep track of whether
   # user specified them or not. If they were un-specified, then
   # the standard init() function is called. If they were specified,
   # the standard check function is called.
   # this allows someone to set control to NULL and revert
   # to using a default control for the given forest.

   data <- args$data %||% self$data

   if("formula" %in% names(args)){
    formula <- args[['formula']]
    terms_old <- terms(self$formula, data = data)
    self$formula <- stats::update(as.formula(terms_old), new = formula)
   }

   null_defaults <- c(
    control = "control",
    weights = "weights",
    mtry = "mtry",
    split_rule = "split_rule",
    split_min_stat = "split_min_stat",
    pred_type = "oobag_pred_type",
    pred_horizon = "oobag_pred_horizon",
    oobag_eval_every = "oobag_eval_every",
    oobag_eval_function = "oobag_fun",
    tree_seeds = "tree_seeds"
   )

   hard_defaults <- c(
    n_tree = "n_tree",
    n_split = "n_split",
    n_retry = "n_retry",
    n_thread = "n_thread",
    sample_with_replacement = "sample_with_replacement",
    sample_fraction = "sample_fraction",
    leaf_min_events = "leaf_min_events",
    leaf_min_obs = "leaf_min_obs",
    split_min_events = "split_min_events",
    split_min_obs = "split_min_obs",
    importance_type = "importance",
    importance_max_pvalue = "importance_max_pvalue",
    importance_group_factors = "group_factors",
    na_action = "na_action",
    verbose_progress = "verbose_progress"
   )

   for( i in seq_along(null_defaults) ){

    input_name <- null_defaults[i]

    if(input_name %in% names(args)){

     input <- args[[input_name]]

     r6_name <- names(null_defaults)[i]

     if(is.null(input)){
      private$user_specified[[r6_name]] <- FALSE
     } else {
      self[[r6_name]] <- input
      private$user_specified[[r6_name]] <- TRUE
     }

    }

   }

   # browser()

   for(i in seq_along(hard_defaults)){

    input_name <- hard_defaults[i]

    if(input_name %in% names(args)){

     input <- args[[input_name]]
     r6_name <- names(hard_defaults)[i]

     self[[r6_name]] <- input

    }

   }

   private$init(data = data)

  },

  # Print object with data dependent on tree_type and lincomb_type
  print = function(){

   # TODO: make tree-specific get_ function for this
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

   oobag_stat <- "none"
   oobag_type <- self$oobag_eval_type

   if(!is_empty(self$eval_oobag$stat_values)){

    oobag_stat <- last_value(self$eval_oobag$stat_values)
    oobag_stat <- round_magnitude(oobag_stat)

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
        paste0('                N events: ', private$n_events             ),

       if(self$tree_type == 'classification')
        paste0('               N classes: ', self$n_class             ),

       paste0('                 N trees: ', self$n_tree                   ),
       paste0('      N predictors total: ', n_predictors                  ),
       paste0('   N predictors per node: ', self$mtry                     ),
       paste0(' Average leaves per tree: ', private$mean_leaves           ),
       paste0('Min observations in leaf: ', self$leaf_min_obs             ),

       if(self$tree_type == 'survival')
        paste0('      Min events in leaf: ', self$leaf_min_events         ),

       paste0('          OOB stat value: ', oobag_stat                    ),
       paste0('           OOB stat type: ', oobag_type                    ),
       paste0('     Variable importance: ', self$importance_type          ),
       '\n-----------------------------------------',
       sep = '\n')

   invisible(self)

  },


  # train an untrained ObliqueForest
  #
  #  ... parameters to override defaults in orsf_cpp. Review
  #   prep_cpp_args() before trying to use this, as some input names
  #   are slightly different in the call to cpp.
  #
  # note: ... functionality is not used in exported
  # orsf_train() currently. I think it's better to
  # keep updating functionality in the orsf_update
  # function.
  train = function(...){

   private$compute_means()
   private$compute_stdev()
   private$compute_modes()
   private$compute_bounds()

   private$prep_x()
   private$prep_y()
   private$prep_w()

   # for survival, inputs should be sorted by time
   private$sort_inputs()

   # allow re-training
   if(self$trained){ self$forest <- list() }

   cpp_args <- private$prep_cpp_args(...)

   cpp_output <- do.call(orsf_cpp, args = cpp_args)

   cpp_output$eval_oobag$stat_type <- self$oobag_eval_type

   .dots <- list(...)
   if(length(.dots) > 0)
    for(i in names(.dots)) self[[i]] <- .dots[[i]]

   self$forest <- cpp_output$forest
   self$importance <- cpp_output$importance
   self$pred_oobag <- cpp_output$pred_oobag
   self$eval_oobag <- cpp_output$eval_oobag

   # don't let rows_oobag contain an empty vector, otherwise it
   # will crash R when cpp tries to load the tree later
   empty_oob_rows <- vapply(self$forest$rows_oobag, is_empty, logical(1))

   if(any(empty_oob_rows)) {
    for(i in which(empty_oob_rows)){
     self$forest$rows_oobag[[i]] <- numeric(1)
    }
   }

   if(self$importance_type != 'none'){
    private$clean_importance()
   }

   if(self$pred_type != 'none'){
    private$clean_pred_oobag()
   }

   self$trained <- TRUE

   private$compute_mean_leaves()

   # free up space
   private$x <- NULL
   private$y <- NULL
   private$w <- NULL

  },

  # this is the inverse of train.
  untrain = function(){
   self$forest <- NULL
   self$importance <- NULL
   self$pred_oobag <- NULL
   self$eval_oobag <- NULL
   self$trained <- FALSE
   private$mean_leaves <- 0
   private$data_means <- NULL
   private$data_modes <- NULL
   private$data_stdev <- NULL
   private$data_bounds <- NULL
  },


  # Prediction always returns matrix or array of predictions
  predict = function(new_data,
                     pred_type,
                     pred_horizon,
                     pred_aggregate,
                     pred_simplify,
                     na_action,
                     boundary_checks,
                     n_thread,
                     verbose_progress){

   public_state <- list(data             = self$data,
                        pred_type        = self$pred_type,
                        pred_horizon     = self$pred_horizon,
                        pred_aggregate   = self$pred_aggregate,
                        na_action        = self$na_action,
                        n_thread         = self$n_thread,
                        verbose_progress = self$verbose_progress)

   private_state <- list(data_rows_complete = private$data_rows_complete)

   # run checks before you assign new values to object.
   # otherwise, if a check throws an error, the object will
   # not be restored to its normal state.

   private$check_data(new = TRUE, data = new_data)
   private$check_na_action(new = TRUE, na_action = na_action)
   private$check_var_missing(new = TRUE, data = new_data, na_action)
   private$check_var_values(new = TRUE, data = new_data)
   private$check_units(data = new_data)
   private$check_boundary_checks(boundary_checks)
   private$check_n_thread(n_thread)
   private$check_verbose_progress(verbose_progress)
   private$check_pred_aggregate(pred_aggregate)
   private$check_pred_simplify(pred_simplify)

   # check and/or set self$pred_horizon and self$pred_type
   # with defaults depending on tree type
   private$init_pred(pred_horizon, pred_type, boundary_checks)
   # set the rest
   self$data             <- new_data
   self$na_action        <- na_action
   self$n_thread         <- n_thread
   self$verbose_progress <- verbose_progress
   self$pred_aggregate   <- pred_aggregate

   out <- try(
    expr = {

     private$init_data_rows_complete()
     private$prep_x()
     # y and w do not need to be prepped for prediction,
     # but they need to match orsf_cpp()'s expectations
     private$y <- matrix(0, nrow = 1, ncol = 1)
     private$w <- rep(1, nrow(private$x))

     private$predict_internal(simplify = pred_simplify)

    },
    silent = TRUE
   )

   # return fields we modified to their original state
   private$restore_state(public_state, private_state)

   # free up space
   private$x <- NULL
   private$y <- NULL
   private$w <- NULL

   # errors within the try statement are not expected,
   # but if they do occur we want to make sure the
   # object is restored to its original state.
   if(is_error(out)) stop(out, call. = FALSE)

   return(out)


  },

  # Variable importance
  # returns a named numeric vector with importance values

  compute_vi = function(type_vi,
                        oobag_fun,
                        n_thread,
                        verbose_progress){

   if(!is.null(oobag_fun)){

    private$check_oobag_eval_function(oobag_fun)
    oobag_eval_function <- oobag_fun
    oobag_eval_type <- "user-specified function"

   } else {

    oobag_eval_function <- self$oobag_eval_function
    oobag_eval_type <- self$oobag_eval_type

   }

   private$prep_x()
   private$prep_y()
   private$prep_w()

   private$sort_inputs()

   # oobag should be FALSE for computing importance
   # even though it is called 'oobag' importance.
   # this is because we subset the x/y/w mats to only
   # contain oobag observations in cpp.

   cpp_args <-
    private$prep_cpp_args(
     oobag_eval_function = oobag_eval_function,
     oobag_eval_type = oobag_eval_type,
     importance_type = type_vi,
     pred_type = self$get_pred_type_vi(),
     pred_mode = FALSE,
     pred_aggregate = TRUE,
     oobag = FALSE,
     write_forest = FALSE,
     run_forest = TRUE,
     n_thread = n_thread %||% self$n_thread,
     verbosity = verbose_progress %||% self$verbose_progress
    )

   out <- do.call(orsf_cpp, args = cpp_args)$importance
   rownames(out) <- colnames(private$x)
   out

  },

  # Partial dependence
  # returns a data.table with dependence values

  compute_dependence = function(pd_data,
                                pred_spec,
                                pred_horizon,
                                pred_type,
                                na_action,
                                expand_grid,
                                prob_values,
                                prob_labels,
                                boundary_checks,
                                n_thread,
                                verbose_progress,
                                oobag,
                                type_output){

   na_action <- na_action %||% self$na_action

   public_state <- list(data             = self$data,
                        na_action        = self$na_action,
                        pred_horizon     = self$pred_horizon,
                        n_thread         = self$n_thread,
                        verbose_progress = self$verbose_progress)

   private_state <- list(data_rows_complete = private$data_rows_complete)

   # run checks before you assign new values to object.
   private$check_boundary_checks(boundary_checks)

   type_input <- 'user'

   if(inherits(pred_spec, 'pspec_auto')){

    private$check_var_names(.names = pred_spec,
                         data = private$data_names$x_original,
                         location = "pred_spec")

    pred_spec <- list_init(pred_spec)

    for(i in names(pred_spec)){

     pred_spec[[i]] <- unique(self$get_var_bounds(i))

    }

   } else if (inherits(pred_spec, 'pspec_intr')){

    type_input <- 'intr'

    pairs <- utils::combn(pred_spec, m = 2)

    pred_spec <- vector(mode = 'list', length = ncol(pairs))


    for(i in seq_along(pred_spec)){

     pred_spec[[i]] <- expand.grid(
      self$get_var_bounds(pairs[1,i]),
      self$get_var_bounds(pairs[2,i])
     )

     colnames(pred_spec[[i]]) <- pairs[, i, drop = TRUE]

    }

   } else {

    private$check_pred_spec(pred_spec, boundary_checks)

   }

   private$check_n_thread(n_thread)
   private$check_expand_grid(expand_grid)
   private$check_verbose_progress(verbose_progress)
   private$check_oobag_pred_mode(oobag, label = 'oobag')

   prob_values <- prob_values %||% c(0.025, 0.50, 0.975)
   prob_labels <- prob_labels %||% c('lwr', 'medn', 'upr')

   private$check_prob_values(prob_values)
   private$check_prob_labels(prob_labels)

   if(length(prob_values) != length(prob_labels)){
    stop("prob_values and prob_labels must have the same length.",
         call. = FALSE)
   }

   # oobag=FALSE to match the format of arg in orsf_pd().
   private$check_pred_type(pred_type, oobag = FALSE,
                           context = 'partial dependence')

   pred_type <- pred_type %||% self$pred_type

   private$check_pred_horizon(pred_horizon, boundary_checks, pred_type)

   if(!oobag){
    private$check_data(new = TRUE, data = pd_data)
    # say new = FALSE to prevent na_action = 'pass'
    private$check_na_action(new = FALSE, na_action = na_action)
    private$check_var_missing(new = TRUE, data = pd_data, na_action)
    private$check_units(data = pd_data)
    self$data <- pd_data
   }

   self$na_action <- na_action
   self$n_thread <- n_thread
   self$verbose_progress <- verbose_progress

   out <- try(
    private$compute_dependence_internal(pred_spec = pred_spec,
                                        pred_type = pred_type,
                                        pred_horizon = pred_horizon,
                                        type_input  = type_input,
                                        type_output = type_output,
                                        expand_grid = expand_grid,
                                        prob_labels = prob_labels,
                                        prob_values = prob_values,
                                        oobag = oobag),
    silent = FALSE
   )

   private$restore_state(public_state, private_state)

   # free up space
   private$x <- NULL
   private$y <- NULL
   private$w <- NULL

   # errors within the try statement are not expected,
   # but if they do occur we want to make sure the
   # object is restored to its original state.
   if(is_error(out)) stop(out, call. = FALSE)

   out

  },

  # Variable selection
  # returns a data.table with variable selection info
  select_variables = function(n_predictor_min, verbose_progress){

   public_state <- list(verbose_progress = self$verbose_progress,
                        forest           = self$forest,
                        control          = self$control)

   object_trained <- self$trained

   out <- try(
    private$select_variables_internal(n_predictor_min, verbose_progress)
   )

   private$restore_state(public_state, private_state = NULL)

   # you can pass a specification to orsf_vs, but it will be trained
   # during the variable selection process. Use untrain() to prevent
   # a specification from being modified when it is passed to orsf_vs
   if(!object_trained) self$untrain()

   if(is_error(out)) stop(out, call. = FALSE)

   out


  },

  # Univariate summary of predictors
  # calls both orsf_vi and orsf_pd_oob to create a summary object
  summarize_uni = function(n_variables = NULL,
                           pred_horizon = NULL,
                           pred_type = NULL,
                           importance_type = NULL,
                           verbose_progress = FALSE){

   # check incoming values if they were specified.
   private$check_n_variables(n_variables)
   private$check_verbose_progress(verbose_progress)

   if(!is.null(pred_horizon)){
    private$check_pred_horizon(pred_horizon, boundary_checks = TRUE)
   }

   private$check_pred_type(pred_type, oobag = FALSE,
                        context = 'partial dependence')
   private$check_importance_type(importance_type)

   names_x <- private$data_names$x_original

   # use existing values if incoming ones were not specified
   n_variables     <- n_variables     %||% length(names_x)
   pred_horizon    <- pred_horizon    %||% self$pred_horizon
   pred_type       <- pred_type       %||% self$pred_type
   importance_type <- importance_type %||% self$importance_type

   # bindings for CRAN check
   value <- NULL
   level <- NULL

   # TODO: make this go away. Just sort alphabetically if no importance
   if(importance_type == 'none' && is_empty(self$importance_type))
    stop("importance cannot be 'none' if object does not have variable",
         " importance values.", call. = FALSE)

   vi <- switch(
    importance_type,
    'anova' = orsf_vi_anova(self, group_factors = TRUE,
                            verbose_progress = verbose_progress),
    'negate' = orsf_vi_negate(self, group_factors = TRUE,
                              verbose_progress = verbose_progress),
    'permute' = orsf_vi_permute(self, group_factors = TRUE,
                                verbose_progress = verbose_progress),
    'none' = NULL
   )

   bounds <- private$data_bounds
   fctrs  <- private$data_fctrs
   n_obs  <- self$n_obs

   names_vi <- names(vi) %||% names_x

   pred_spec <- list_init(names_vi)[seq(n_variables)]

   for(i in names(pred_spec)){

    if(i %in% colnames(bounds)){

     pred_spec[[i]] <- self$get_var_bounds(i)
      # unique(
      #  as.numeric(bounds[c('25%','50%','75%'), i])
      # )

    } else if (i %in% fctrs$cols) {

     pred_spec[[i]] <- fctrs$lvls[[i]]

    }

   }

   pd_output <- orsf_pd_oob(object = self,
                            pred_spec = pred_spec,
                            expand_grid = FALSE,
                            pred_type = pred_type,
                            prob_values = c(0.25, 0.50, 0.75),
                            pred_horizon = pred_horizon,
                            boundary_checks = FALSE,
                            verbose_progress = verbose_progress)

   fctrs_unordered <- c()

   # did the orsf have factor variables?
   if(!is_empty(fctrs$cols)){
    fctrs_unordered <- fctrs$cols[!fctrs$ordr]
   }

   # some cart-wheels here for backward compatibility.
   f <- as.factor(pd_output$variable)

   name_rep <- rle(as.integer(f))

   pd_output$importance <- rep(vi[levels(f)[name_rep$values]],
                               times = name_rep$lengths)

   pd_output[, value := fifelse(test = is.na(value),
                                yes = as.character(level),
                                no = round_magnitude(value))]

   # if a := is used inside a function with no DT[] before the end of the
   # function, then the next time DT or print(DT) is typed at the prompt,
   # nothing will be printed. A repeated DT or print(DT) will print.
   # To avoid this: include a DT[] after the last := in your function.
   pd_output[]

   new_order <- c('variable', 'importance', 'value',
                  'mean', 'medn', 'lwr', 'upr')

   if(self$tree_type == 'classification'){
    new_order <- insert_vals(new_order, 2, 'class')
   }

   setcolorder(pd_output, new_order)

   structure(
    .Data = list(dt = pd_output,
                 pred_type = pred_type,
                 pred_horizon = pred_horizon),
    class = 'orsf_summary_uni'
   )


  },

  # setter

  set_field = function(...){

   .dots <- list(...)

   valid_fields <- c("data",
                     "formula",
                     "control",
                     "weights",
                     "n_tree",
                     "n_split",
                     "n_retry",
                     "n_thread",
                     "mtry",
                     "sample_with_replacement",
                     "sample_fraction",
                     "leaf_min_events",
                     "leaf_min_obs",
                     "split_rule",
                     "split_min_events",
                     "split_min_obs",
                     "split_min_stat",
                     "pred_type",
                     "oobag_pred_horizon",
                     "oobag_eval_every",
                     "oobag_fun",
                     "importance_type",
                     "importance_max_pvalue",
                     "importance_group_factors",
                     "tree_seeds",
                     "na_action",
                     "verbose_progress")

   if(!is_empty(.dots)){

    for(i in names(.dots)){

     if(!(i %in% valid_fields))
      stop("field ", i, " is unrecognized")

     self[[i]] <- .dots[[i]]

    }

   }

  },

  # getters
  # the get_ functions return a private member variable in
  # standard or requested format, allowing for info about
  # specific variables to be returned.
  get_names_x = function(ref_coded = FALSE){

   if(ref_coded) return(private$data_names$x_ref_code)

   return(private$data_names$x_original)

  },

  get_names_y = function(){
   return(private$data_names$y)
  },

  get_var_bounds = function(.name){

   if(.name %in% private$data_names$x_numeric){

    out <- unique(as.numeric(private$data_bounds[, .name]))

    if(length(out) < 5){
     # too few unique values to use quantiles,
     # so use the most common unique values instead.
     unis <- sort(table(self$data[[.name]]), decreasing = TRUE)
     n_items <- min(5, length(unis))
     out <- sort(as.numeric(names(unis)[seq(n_items)]))
    }

    return(out)

   } else {

    return(private$data_fctrs$lvls[[.name]])

   }

  },

  get_var_type = function(.name){
   return(class(self$data[[.name]])[1])
  },

  get_var_unit = function(.name){
   return(private$data_units[[.name]])
  },

  get_fctr_info = function(){
   return(private$data_fctrs)
  },

  get_importance_raw = function(){
   return(private$importance_raw)
  },

  get_importance_clean = function(importance_raw = NULL,
                                  group_factors = NULL){


   input <- importance_raw %||% private$importance_raw
   group_factors <- group_factors %||% self$group_factors
   private$clean_importance(input, group_factors)

   return(self$importance)

  },

  get_mean_leaves_per_tree = function(){
   return(private$mean_leaves)
  },

  get_means = function(){
   return(private$data_means)
  },

  get_modes = function(){
   return(private$data_modes)
  },

  get_stdev = function(){
   return(private$data_stdev)
  },

  get_bounds = function(){
   return(private$data_bounds)
  }

 ),

 # private ----
 private = list(

  user_specified = NULL,
  data_rows_complete = NULL,
  data_types = NULL,
  data_names = NULL,
  data_fctrs = NULL,
  data_units = NULL,
  data_means = NULL,
  data_stdev = NULL,
  data_modes = NULL,
  data_bounds = NULL,

  x = NULL,
  y = NULL,
  w = NULL,

  importance_raw = NULL,

  mean_leaves = 0,

  # checkers
  check_data = function(data = NULL, new = FALSE){

   # additional data checks are run during initialization.
   input <- data %||% self$data

   check_arg_is(arg_value = input,
                arg_name = 'data',
                expected_class = 'data.frame')

   # Minimal checks for now, other checks occur later
   if(nrow(input) == 0 || ncol(input) ==  0){
    data_label <- if(new) "new" else "training"
    stop(data_label, " data are empty",
         call. = FALSE)

   }

   if(!new){

    # check for blanks first b/c the check for non-standard symbols
    # will detect blanks with >1 empty characters
    blank_names <- grepl(pattern = '^\\s*$', x = names(input))

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

    ns_names <- grepl(pattern = '[^a-zA-Z0-9\\.\\_]+', x = names(input))

    if(any(ns_names)){

     last <- ifelse(sum(ns_names) == 2, ' and ', ', and ')

     stop("Non-standard names detected in training data: ",
          paste_collapse(x = names(input)[ns_names],
                         last = last),
          call. = FALSE)

    }

   } else {

    private$check_var_names_new(new_data = input)
    private$check_var_types_new(new_data = input)
    private$check_fctrs_new(new_data = input)

   }


  },


  check_var_names = function(.names,
                             data = NULL,
                             location = "formula"){

   data <- data %||% self$data

   if(is.character(data)){
    data_names <- data
   } else {
    data_names <- names(data)
   }

   names_not_found <- setdiff(c(.names), data_names)

   if(!is_empty(names_not_found)){
    msg <- paste0(
     "variables in ", location, " were not found in data: ",
     paste_collapse(names_not_found, last = ' and ')
    )
    stop(msg, call. = FALSE)
   }

  },

  check_var_names_new = function(new_data,
                                 check_new_in_ref = FALSE,
                                 check_ref_in_new = TRUE){

   check_new_data_names(new_data,
                        ref_names = private$data_names$x_original,
                        label_new = "new_data",
                        label_ref = 'training data',
                        check_new_in_ref = check_new_in_ref,
                        check_ref_in_new = check_ref_in_new)


  },

  # Check variable types
  #
  # orsf() should only be run with certain types of variables. This function
  #   checks input data to make sure all variables have a primary (i.e., first)
  #   class that is within the list of valid options.
  #
  # return an error if something is wrong.
  check_var_types = function(var_types, var_names, valid_types){

   good_vars <- var_types %in% valid_types

   if(!all(good_vars)){

    bad_vars <- which(!good_vars)

    vars_to_list <- var_names[bad_vars]
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

  check_var_types_new = function(new_data){

   check_new_data_types(new_data,
                        ref_names = private$data_names$x_original,
                        ref_types = private$data_types$x,
                        label_new = "new_data",
                        label_ref = "training data")

  },

  check_var_missing = function(data = NULL,
                               new = FALSE,
                               na_action = NULL){

   input <- data %||% self$data
   na_action <- na_action %||% self$na_action

   if(na_action == 'fail'){

    if(any(is.na(select_cols(input, private$data_names$x_original)))){
     data_label <- if(new) "new data" else "training data"
     stop("Please remove missing values from ", data_label, ", or impute them.",
          call. = FALSE)
    }

   }

   if(na_action == 'omit'){

    if(length(private$data_rows_complete) == 0){

     data_label <- if(new) "new data" else "training data"
     stop("There are no observations in ",
          data_label, " with complete data ",
          "for the predictors used by this orsf object.",
          call. = FALSE)

    }

   }

   if(!new){

    if(any(is.na(select_cols(input, private$data_names$y))))
     stop("Please remove missing values from the outcome variable(s).",
          call. = FALSE)

   }

  },

  check_var_values = function(data = NULL, new = FALSE){

   input <- data %||% self$data

   for(i in private$data_names$x_original){

    if(collapse::allNA(input[[i]])){
     stop("column ", i, " has no observed values.",
          call. = FALSE)
    }

    if(any(is.infinite(input[[i]]))){
     stop("Please remove infinite values from ", i, ".",
          call. = FALSE)
    }

    if(!new){
     if(collapse::fnunique(collapse::na_omit(input[[i]])) == 1L){
      stop("column ", i, " is constant.",
           call. = FALSE)
     }
    }

   }

  },

  check_fctrs_new = function(new_data){

   check_new_data_fctrs(new_data,
                        names_x = private$data_names$x_original,
                        fi_ref = private$data_fctrs,
                        label_new = "new_data")

  },

  check_formula = function(formula = NULL){

   input <- formula %||% self$formula

   check_arg_is(arg_value = input,
                arg_name = 'formula',
                expected_class = 'formula')

   if(length(input) != 3){
    stop("formula must be two sided, i.e. left side ~ right side",
         call. = FALSE)
   }

   formula_deparsed <- as.character(input)[[3]]

   for( symbol in c("*", "^", ":", "(", ")", "["," ]", "|", "%") ){

    if(grepl(symbol, formula_deparsed, fixed = TRUE)){

     stop("unrecognized symbol in formula: ", symbol,
          "\norsf recognizes '+', '-', and '.' symbols.",
          call. = FALSE)

    }

   }

  },
  check_control = function(control = NULL){

   input <- control %||% self$control

   check_arg_is(arg_value = input,
                arg_name = 'control',
                expected_class = 'orsf_control')

   if(input$lincomb_type == 'net'){

    if (!requireNamespace("glmnet", quietly = TRUE)) {
     stop(
      "Package \"glmnet\" must be installed to use",
      " orsf_control_net() with orsf().",
      call. = FALSE
     )
    }

   }

  },
  check_weights = function(weights = NULL){

   input <- weights %||% self$weights

   check_arg_type(arg_value = input,
                  arg_name = 'weights',
                  expected_type = 'numeric')

   check_arg_gteq(arg_value = input,
                  arg_name = 'weights',
                  bound = 0)

   check_arg_length(arg_value = input,
                    arg_name  = 'weights',
                    expected_length = self$n_obs)

  },
  check_n_tree = function(n_tree = NULL){

   input <- n_tree %||% self$n_tree

   check_arg_type(arg_value = input,
                  arg_name = 'n_tree',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'n_tree')

   check_arg_gteq(arg_value = input,
                  arg_name = 'n_tree',
                  bound = 1)

   check_arg_lteq(arg_value = input,
                  arg_name = 'n_tree',
                  bound = 10000)

   check_arg_length(arg_value = input,
                    arg_name = 'n_tree',
                    expected_length = 1)

  },
  check_n_split = function(n_split = NULL){

   input <- n_split %||% self$n_split

   check_arg_type(arg_value = input,
                  arg_name = 'n_split',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'n_split')

   check_arg_gteq(arg_value = input,
                  arg_name = 'n_split',
                  bound = 1)

   check_arg_length(arg_value = input,
                    arg_name = 'n_split',
                    expected_length = 1)

  },
  check_n_retry = function(n_retry = NULL){

   input <- n_retry %||% self$n_retry

   check_arg_type(arg_value = input,
                  arg_name = 'n_retry',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'n_retry')

   check_arg_gteq(arg_value = input,
                  arg_name = 'n_retry',
                  bound = 0)

   check_arg_length(arg_value = input,
                    arg_name = 'n_retry',
                    expected_length = 1)

  },
  check_n_thread = function(n_thread = NULL){

   input <- n_thread %||% self$n_thread

   check_arg_type(arg_value = input,
                  arg_name = 'n_thread',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'n_thread')

   check_arg_gteq(arg_value = input,
                  arg_name = 'n_thread',
                  bound = 0)

   check_arg_length(arg_value = input,
                    arg_name = 'n_thread',
                    expected_length = 1)

  },

  check_n_variables = function(n_variables = NULL){

   # n_variables is not a field of ObliqueForest,
   # so it is only checked as an incoming input.

   if(!is.null(n_variables)){

    check_arg_type(arg_value = n_variables,
                   arg_name = 'n_variables',
                   expected_type = 'numeric')

    check_arg_is_integer(arg_value = n_variables,
                         arg_name = 'n_variables')

    check_arg_gteq(arg_value = n_variables,
                   arg_name = 'n_variables',
                   bound = 1)

    check_arg_lteq(arg_value = n_variables,
                   arg_name = 'n_variables',
                   bound = length(private$data_names$x_original),
                   append_to_msg = "(total number of predictors)")


    check_arg_length(arg_value = n_variables,
                     arg_name = 'n_variables',
                     expected_length = 1)

   }

  },

  check_mtry = function(mtry = NULL){

   input <- mtry %||% self$mtry

   n_predictors <- length(private$data_names$x_ref_code)

   # okay for this to be unspecified at startup
   if(!is.null(input)){

    check_arg_type(arg_value = input,
                   arg_name = 'mtry',
                   expected_type = 'numeric')

    check_arg_is_integer(arg_value = input,
                         arg_name = 'mtry')

    check_arg_gteq(arg_value = input,
                   arg_name = 'mtry',
                   bound = 1)

    check_arg_length(arg_value = input,
                     arg_name = 'mtry',
                     expected_length = 1)

    if(!is.null(n_predictors)){
     check_arg_lteq(
      arg_value = input,
      arg_name = 'mtry',
      bound = n_predictors,
      append_to_msg = "(number of columns in the reference coded x-matrix)"
     )
    }

   }

  },
  check_sample_with_replacement = function(sample_with_replacement = NULL){

   input <- sample_with_replacement %||% self$sample_with_replacement

   check_arg_type(arg_value = input,
                  arg_name = 'sample_with_replacement',
                  expected_type = 'logical')

   check_arg_length(arg_value = input,
                    arg_name = 'sample_with_replacement',
                    expected_length = 1)

  },
  check_sample_fraction = function(sample_fraction = NULL){

   input <- sample_fraction %||% self$sample_fraction

   check_arg_type(arg_value = input,
                  arg_name = 'sample_fraction',
                  expected_type = 'numeric')

   check_arg_gt(arg_value = input,
                arg_name = 'sample_fraction',
                bound = 0)

   check_arg_lteq(arg_value = input,
                  arg_name = 'sample_fraction',
                  bound = 1)

   check_arg_length(arg_value = input,
                    arg_name = 'sample_fraction',
                    expected_length = 1)

  },
  check_leaf_min_obs = function(leaf_min_obs = NULL,
                                n_obs = NULL){

   input <- leaf_min_obs %||% self$leaf_min_obs
   n_obs <- n_obs %||% self$n_obs

   # users should never get this error but it may help me debug
   if(is.null(n_obs))
    stop("cannot check leaf_min_obs when n_obs is unspecified")

   check_arg_type(arg_value = input,
                  arg_name = 'leaf_min_obs',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'leaf_min_obs')

   check_arg_gteq(arg_value = input,
                  arg_name = 'leaf_min_obs',
                  bound = 1)

   check_arg_length(arg_value = input,
                    arg_name = 'leaf_min_obs',
                    expected_length = 1)

   check_arg_lteq(arg_value = input,
                  arg_name = "leaf_min_obs",
                  bound = round(n_obs / 2),
                  append_to_msg = "(number of observations divided by 2)")

  },
  check_split_rule = function(split_rule = NULL){

   input <- split_rule %||% self$split_rule

   # okay to pass a NULL value on startup
   if(!is.null(input)){

    check_arg_type(arg_value = input,
                   arg_name = 'split_rule',
                   expected_type = 'character')

    check_arg_length(arg_value = input,
                     arg_name = 'split_rule',
                     expected_length = 1)

    private$check_split_rule_internal()

   }


  },
  check_split_rule_internal = function(){

   stop("this method should be defined in a derived class.")

  },
  check_split_min_obs = function(split_min_obs = NULL,
                                 n_obs = NULL){

   input <- split_min_obs %||% self$split_min_obs
   n_obs <- n_obs %||% self$n_obs

   # users should never get this error but it may help me debug
   if(is.null(n_obs))
    stop("cannot check split_min_obs when n_obs is unspecified")

   check_arg_type(arg_value = input,
                  arg_name = 'split_min_obs',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'split_min_obs')

   check_arg_gteq(arg_value = input,
                  arg_name = 'split_min_obs',
                  bound = 1)

   check_arg_length(arg_value = input,
                    arg_name = 'split_min_obs',
                    expected_length = 1)

   check_arg_lt(arg_value = input,
                arg_name = "split_min_obs",
                bound = n_obs,
                append_to_msg = "(number of observations)")

  },
  check_split_min_stat = function(split_min_stat = NULL,
                                  split_rule = NULL){

   input <- split_min_stat %||% self$split_min_stat

   split_rule <- split_rule %||% self$split_rule

   # okay to pass a NULL value on startup
   if(!is.null(input)){

    check_arg_type(arg_value = input,
                   arg_name = 'split_min_stat',
                   expected_type = 'numeric')

    check_arg_gteq(arg_value = input,
                   arg_name = 'split_min_stat',
                   bound = 0)

    check_arg_length(arg_value = input,
                     arg_name = 'split_min_stat',
                     expected_length = 1)

    # also okay for split_rule to be NULL on startup
    if(!is.null(split_rule)){

     if(split_rule %in% c('cstat', 'gini')){

      append <- paste0("(split stat <", self$split_rule, "> is always < 1)")

      check_arg_lt(arg_value = input,
                   arg_name = 'split_min_stat',
                   bound = 1,
                   append_to_msg = append)

     }

    }

   }

  },

  # must specify oobag when you call this to make sure it isn't forgotten
  check_pred_type = function(pred_type = NULL, oobag, context = NULL){

   input <- pred_type %||% self$pred_type

   # okay to pass a NULL value on startup
   if(!is.null(input)){

    arg_name <- if(oobag) 'oobag_pred_type' else "pred_type"

    check_arg_type(arg_value = input,
                   arg_name = arg_name,
                   expected_type = 'character')

    check_arg_length(arg_value = input,
                     arg_name = arg_name,
                     expected_length = 1)

    private$check_pred_type_internal(oobag = oobag,
                                  pred_type = pred_type,
                                  context = context)

   }



  },

  check_pred_aggregate = function(pred_aggregate = NULL){

   input <- pred_aggregate %||% self$pred_aggregate

   check_arg_type(arg_value = input,
                  arg_name = 'pred_aggregate',
                  expected_type = 'logical')

   check_arg_length(arg_value = input,
                    arg_name = 'pred_aggregate',
                    expected_length = 1)

  },

  check_pred_simplify = function(pred_simplify){

   check_arg_type(arg_value = pred_simplify,
                  arg_name = 'pred_simplify',
                  expected_type = 'logical')

   check_arg_length(arg_value = pred_simplify,
                    arg_name = 'pred_simplify',
                    expected_length = 1)

  },

  check_oobag_eval_every = function(oobag_eval_every = NULL,
                                    n_tree = NULL){

   input <- oobag_eval_every %||% self$oobag_eval_every

   n_tree <- n_tree %||% self$n_tree

   check_arg_type(arg_value = input,
                  arg_name = 'oobag_eval_every',
                  expected_type = 'numeric')

   check_arg_is_integer(arg_value = input,
                        arg_name = 'oobag_eval_every')

   check_arg_gteq(arg_value = input,
                  arg_name = 'oobag_eval_every',
                  bound = 1)

   check_arg_length(arg_value = input,
                    arg_name = 'oobag_eval_every',
                    expected_length = 1)

   check_arg_lteq(arg_value = self$oobag_eval_every,
                  arg_name = 'oobag_eval_every',
                  bound = n_tree)


  },

  check_importance_type = function(importance_type = NULL){

   input <- importance_type <- self$importance_type

   check_arg_type(arg_value = input,
                  arg_name = 'importance',
                  expected_type = 'character')

   check_arg_length(arg_value = input,
                    arg_name = 'importance',
                    expected_length = 1)

   check_arg_is_valid(arg_value = input,
                      arg_name = 'importance',
                      valid_options = c("none",
                                        "anova",
                                        "negate",
                                        "permute"))

  },
  check_importance_max_pvalue = function(importance_max_pvalue = NULL){

   input <- importance_max_pvalue %||% self$importance_max_pvalue

   check_arg_type(arg_value = input,
                  arg_name = 'importance_max_pvalue',
                  expected_type = 'numeric')

   check_arg_gt(arg_value = input,
                arg_name = 'importance_max_pvalue',
                bound = 0)

   check_arg_lt(arg_value = input,
                arg_name = 'importance_max_pvalue',
                bound = 1)

   check_arg_length(arg_value = input,
                    arg_name = 'importance_max_pvalue',
                    expected_length = 1)


  },

  check_importance_group_factors = function(importance_group_factors = NULL){

   input <- importance_group_factors %||% self$importance_group_factors

   check_arg_type(arg_value = input,
                  arg_name = 'group_factors',
                  expected_type = 'logical')

   check_arg_length(arg_value = input,
                    arg_name = 'group_factors',
                    expected_length = 1)

  },
  check_tree_seeds = function(tree_seeds = NULL,
                              n_tree = NULL){

   input <- tree_seeds %||% self$tree_seeds
   n_tree <- n_tree %||% self$n_tree

   # okay for this to be unspecified at start-up
   if(!is.null(input)){

    check_arg_type(arg_value = input,
                   arg_name = 'tree_seed',
                   expected_type = 'numeric')

    check_arg_is_integer(arg_value = input,
                         arg_name = 'tree_seeds')

    if(!is.null(n_tree)){

     if(length(input) > 1 && length(input) != n_tree){

      stop('tree_seeds should have length = 1 or length = ",
          "n_tree <', self$n_tree,
           "> (the number of trees) but instead has length <",
           length(input), ">", call. = FALSE)

     }

    }

   }

  },
  check_na_action = function(na_action = NULL, new = FALSE){

   input <- na_action %||% self$na_action

   check_arg_type(arg_value = input,
                  arg_name = 'na_action',
                  expected_type = 'character')

   check_arg_length(arg_value = input,
                    arg_name = 'na_action',
                    expected_length = 1)

   valid_options <- c("fail", "omit", "impute_meanmode")

   if(new) valid_options <- c(valid_options, 'pass')

   check_arg_is_valid(arg_value = input,
                      arg_name = 'na_action',
                      valid_options = valid_options)

  },
  check_oobag_eval_function = function(oobag_eval_function = NULL){

   input <- oobag_eval_function %||% self$oobag_eval_function

   if(!is.null(input)){

    oobag_fun_args <- names(formals(input))

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

    test_output <- private$check_oobag_eval_function_internal(input)

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


  check_lincomb_R_function = function(lincomb_R_function = NULL){

   input <- lincomb_R_function %||% self$control$lincomb_R_function

   args <- names(formals(input))

   if(length(args) != 3) stop(
    "input should have 3 input arguments but instead has ",
    length(args),
    call. = FALSE
   )

   arg_names_expected <- c("x_node",
                           "y_node",
                           "w_node")

   arg_names_refer <- c('first', 'second', 'third')

   for(i in seq_along(arg_names_expected)){
    if(args[i] != arg_names_expected[i])
     stop(
      "the ", arg_names_refer[i], " input argument of input ",
      "should be named '", arg_names_expected[i],"' ",
      "but is instead named '", args[i], "'",
      call. = FALSE
     )
   }

   test_output <- private$check_lincomb_R_function_internal(input)

   if(!is.matrix(test_output)) stop(
    "user-supplied function should return a matrix output ",
    "but instead returns output of type ", class(test_output)[1],
    call. = FALSE
   )

   if(ncol(test_output) != 1) stop(
    "user-supplied function should return a matrix with 1 column ",
    "but instead returns a matrix with ", ncol(test_output), " columns.",
    call. = FALSE
   )

   if(nrow(test_output) != 3L) stop(
    "user-supplied function should return a matrix with 1 row for each ",
    " column in x_node but instead returns a matrix with ",
    nrow(test_output), " rows ", "in a testing case where x_node has ",
    3L, " columns",
    call. = FALSE
   )

  },

  check_verbose_progress = function(verbose_progress = NULL){

   input <- verbose_progress %||% self$verbose_progress

   # check_arg_type(arg_value = input,
   #                arg_name = 'verbose_progress',
   #                expected_type = 'logical')
   #
   # check_arg_length(arg_value = input,
   #                  arg_name = 'verbose_progress',
   #                  expected_length = 1)

  },
  check_units = function(data = NULL){

   input <- data %||% self$data

   ui_new <- unit_info(data = input, .names = names(private$data_units))

   ui_missing <- setdiff(names(private$data_units), names(ui_new))

   if(!is_empty(ui_missing)){

    if(length(ui_missing) == 1){
     stop(ui_missing, " had unit attributes in training data but",
          " did not have unit attributes in testing data.",
          " Please ensure that variables in new data have the same",
          " units as their counterparts in the training data.",
          call. = FALSE)
    }

    stop(length(ui_missing),
         " variables (",
         paste_collapse(ui_missing, last = ' and '),
         ") had unit attributes in training",
         " data but did not have unit attributes in new data.",
         " Please ensure that variables in new data have the same",
         " units as their counterparts in the training data.",
         call. = FALSE)

   }

   for(i in names(private$data_units)){

    if(private$data_units[[i]]$label != ui_new[[i]]$label){

     msg <- paste0("variable <", i, "> has unit '",
                   private$data_units[[i]]$label,
                   "' in training data but has unit '",
                   ui_new[[i]]$label, "' in new data")

     stop(msg, call. = FALSE)

    }

   }

  },

  check_pred_spec = function(pred_spec, boundary_checks){

   if(is_empty(pred_spec)){
    stop("pred_spec is empty", call. = FALSE)
   }

   if(is_empty(names(pred_spec))){
    stop("pred_spec is unnamed", call. = FALSE)
   }

   bad_name_index <- which(
    is.na(match(names(pred_spec), private$data_names$x_original))
   )

   if(!is_empty(bad_name_index)){

    bad_names <- names(pred_spec)[bad_name_index]

    stop("variables in pred_spec are not recognized as predictors in object: ",
         paste_collapse(bad_names, last = ' and '),
         call. = FALSE)

   }

   numeric_bounds <- private$data_bounds
   numeric_names <- intersect(colnames(numeric_bounds), names(pred_spec))

   if(!is_empty(numeric_names) && boundary_checks){

    for(.name in numeric_names){

     vals_above_stop <- which(pred_spec[[.name]] > numeric_bounds['90%', .name])
     vals_below_stop <- which(pred_spec[[.name]] < numeric_bounds['10%', .name])

     boundary_error <- FALSE
     vals_above_list <- vals_below_list <- " "

     if(!is_empty(vals_above_stop)){
      vals_above_list <- paste_collapse(
       round_magnitude(pred_spec[[.name]][vals_above_stop]),
       last = ' and '
      )

      boundary_error <- TRUE
      vals_above_list <-
       paste0(" (",vals_above_list," > ", numeric_bounds['90%', .name],") ")

     }

     if(!is_empty(vals_below_stop)){

      vals_below_list <- paste_collapse(
       round_magnitude(pred_spec[[.name]][vals_below_stop]),
       last = ' and '
      )

      boundary_error <- TRUE

      vals_below_list <-
       paste0(" (",vals_below_list," < ", numeric_bounds['10%', .name],") ")

     }

     if(boundary_error)
      stop("Some values for ",
           .name,
           " in pred_spec are above",
           vals_above_list,
           "or below",
           vals_below_list,
           "90th or 10th percentiles in training data.",
           " Change pred_spec or set boundary_checks = FALSE",
           " to prevent this error",
           call. = FALSE)

    }

   }



  },

  check_boundary_checks = function(boundary_checks){

   # not a field so boundary_checks should never be null

   check_arg_type(arg_value = boundary_checks,
                  arg_name = 'boundary_checks',
                  expected_type = 'logical')

   check_arg_length(arg_value = boundary_checks,
                    arg_name = 'boundary_checks',
                    expected_length = 1)


  },

  check_expand_grid = function(expand_grid){

   # not a field so boundary_checks should never be null
   check_arg_type(arg_value = expand_grid,
                  arg_name = 'expand_grid',
                  expected_type = 'logical')

   check_arg_length(arg_value = expand_grid,
                    arg_name = 'expand_grid',
                    expected_length = 1)

  },

  check_prob_values = function(prob_values){

   check_arg_type(arg_value = prob_values,
                  arg_name = 'prob_values',
                  expected_type = 'numeric')

   check_arg_gteq(arg_value = prob_values,
                  arg_name = 'prob_values',
                  bound = 0)

   check_arg_lteq(arg_value = prob_values,
                  arg_name = 'prob_values',
                  bound = 1)

  },

  check_prob_labels = function(prob_labels){


   check_arg_type(arg_value = prob_labels,
                  arg_name = 'prob_labels',
                  expected_type = 'character')

  },

  check_oobag_pred_mode = function(oobag_pred_mode, label,
                                   sample_fraction = NULL){

   sample_fraction <- sample_fraction %||% self$sample_fraction

   check_arg_type(arg_value = oobag_pred_mode,
                  arg_name = label,
                  expected_type = 'logical')

   check_arg_length(arg_value = oobag_pred_mode,
                    arg_name = label,
                    expected_length = 1)

   if(!is.null(sample_fraction)){

    if(oobag_pred_mode && sample_fraction == 1){
     stop(
      "cannot compute out-of-bag predictions if no samples are out-of-bag.",
      " Try setting sample_fraction < 1 or oobag_pred_type = 'none'.",
      call. = FALSE
     )
    }

   }

  },

  check_lincomb_df_target = function(lincomb_df_target = NULL,
                                     mtry = NULL){

   input <- lincomb_df_target %||% self$control$lincomb_df_target
   mtry <- mtry %||% self$mtry

   check_arg_lteq(
    arg_value = input,
    arg_name = 'df_target',
    bound = mtry,
    append_to_msg = "(number of randomly selected predictors)"
   )

  },

  # runs checks and sets defaults where needed.
  # data is NULL when we are creating a new forest,
  # but may be non-NULL if we update an existing one
  init = function(data = NULL) {

   # look for odd symbols in formula before you check variables in data
   private$check_formula()
   # check & init data should be near first bc they set up other checks
   private$check_data(data)
   private$init_data(data)

   # if data is not null, it means we are updating an orsf spec
   # and in that process applying it to a new dataset, so:
   if(!is.null(data)) self$data <- data

   if(private$user_specified$control){
    private$check_control()
   } else {
    private$init_control()
   }


   if(private$user_specified$mtry){
    private$check_mtry()
   } else {
    private$init_mtry()
   }

   if(private$user_specified$lincomb_df_target){
    private$check_lincomb_df_target()
   } else {
    private$init_lincomb_df_target()
   }

   if(private$user_specified$weights){
    private$check_weights()
   } else {
    private$init_weights()
   }

   if(private$user_specified$pred_type){
    private$check_pred_type(oobag = TRUE)
   } else {
    private$init_pred_type()
   }

   if(private$user_specified$split_rule){
    private$check_split_rule()
   } else {
    private$init_split_rule()
   }

   if(private$user_specified$split_min_stat){
    private$check_split_min_stat()
   } else {
    private$init_split_min_stat()
   }

   if(private$user_specified$oobag_eval_function){
    private$check_oobag_eval_function()
    self$oobag_eval_type <- "User-specified function"
   } else {
    private$init_oobag_eval_function()
   }

   if(private$user_specified$oobag_eval_every){
    private$check_oobag_eval_every()
   } else {
    private$init_oobag_eval_every()
   }

   if(self$control$lincomb_type == 'custom'){
    private$check_lincomb_R_function()
   } else if (is.null(self$control$lincomb_R_function)){
    private$init_lincomb_R_function()
   }

   # arguments with hard defaults do not need an init option
   private$check_n_tree()
   private$check_n_split()
   private$check_n_retry()
   private$check_n_thread()
   private$check_sample_with_replacement()
   private$check_sample_fraction()
   private$check_leaf_min_obs()
   private$check_split_min_obs()
   private$check_importance_type()
   private$check_importance_max_pvalue()
   private$check_importance_group_factors()
   private$check_na_action()
   private$check_verbose_progress()

   # args below depend on at least one upstream arg

   if(private$user_specified$tree_seeds && is.null(self$forest_seed)){
    private$check_tree_seeds()
   } else if (!is.null(self$forest_seed)){
    # this only happens when the forest is updated
    private$plant_tree_seeds(self$forest_seed)
   } else {
    private$init_tree_seeds()
   }

   if(length(self$tree_seeds) == 1 && self$n_tree > 1){
    private$plant_tree_seeds(self$tree_seeds)
   }

   # oobag_pred_mode depends on pred_type, which is checked above,
   # so there is no reason to check it here.
   private$init_oobag_pred_mode()
   # check if sample_fraction conflicts with oobag_pred_mode
   private$check_oobag_pred_mode(self$oobag_pred_mode,
                              label = 'oobag_pred_mode',
                              sample_fraction = self$sample_fraction)

   private$init_internal()


  },
  init_data = function(data = NULL){

   if(!is.null(data)) self$data <- data

   formula_terms <- suppressWarnings(
    stats::terms(x = self$formula, data = self$data)
   )

   if(attr(formula_terms, 'response') == 0)
    stop("formula should have a response", call. = FALSE)

   names_x_data <- attr(formula_terms, 'term.labels')
   names_y_data <- all.vars(self$formula[[2]])

   # check factors in x data
   fctr_check(self$data, names_x_data)
   fctr_id_check(self$data, names_x_data)

   private$check_var_names(c(names_x_data, names_y_data), data = self$data)

   private$data_names <- list(y = names_y_data,
                              x_original = names_x_data)


   types_y_data <- vapply(names_y_data, self$get_var_type, character(1))
   types_x_data <- vapply(names_x_data, self$get_var_type, character(1))

   valid_y_types <- c('numeric', 'integer', 'units', 'factor', "Surv")
   valid_x_types <- c(valid_y_types, 'ordered')

   private$check_var_types(types_y_data, names_y_data, valid_y_types)
   private$check_var_types(types_x_data, names_x_data, valid_x_types)

   private$data_types <- list(y = types_y_data, x = types_x_data)
   private$data_fctrs <- fctr_info(self$data, private$data_names$x_original)

   private$init_numeric_names()
   private$init_ref_code_names()
   private$init_data_rows_complete()

   self$n_obs <- ifelse(test = self$na_action == 'omit',
                        yes  = length(private$data_rows_complete),
                        no   = nrow(self$data))

   private$check_var_missing()
   private$check_var_values()

   unit_names <- c(names_y_data[types_y_data == 'units'],
                   names_x_data[types_x_data == 'units'])

   private$data_units <- unit_info(data = self$data, .names = unit_names)

  },
  init_data_rows_complete = function(){

   private$data_rows_complete <- collapse::whichv(
    collapse::missing_cases(self$data, private$data_names$x_original),
    value = FALSE
   )

  },
  init_tree_seeds = function(){

   if(is.null(self$tree_seeds)){
    self$tree_seeds <- sample(1e6, size = 1)
   }

  },
  init_numeric_names = function(){

   pattern <- "^integer$|^numeric$|^units$"

   index <- grep(pattern = pattern, x = private$data_types$x)

   private$data_names[["x_numeric"]] = private$data_names$x_original[index]

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

   # if pred_type is null when this is run, it means
   # the user did not specify pred_type, which means the
   # family-specific default will be used, which means
   # pred_type will not be 'none', so it is safe to assume
   # oobag_pred_mode is TRUE if pred_type is currently null

   if(is.null(self$pred_type)){
    self$oobag_pred_mode <- TRUE
   } else {
    self$oobag_pred_mode <- self$pred_type != "none"
   }

  },
  init_mtry = function(){

   n_col_x <- length(private$data_names$x_ref_code)

   self$mtry <- ceiling(sqrt(n_col_x))

  },




  init_lincomb_df_target = function(mtry = NULL){

   mtry <- mtry %||% self$mtry

   self$control$lincomb_df_target <- mtry

  },

  init_weights = function(){

   # set weights as 1 if user did not supply them.
   # length of weights depends on how missing are handled.
   self$weights <- rep(1, self$n_obs)

  },


  init_oobag_eval_function = function(){

   self$oobag_eval_function <- function(y_mat, w_vec, s_vec){
    return(1)
   }

  },

  init_oobag_eval_every = function(n_tree = NULL){
   self$oobag_eval_every = n_tree %||% self$n_tree
  },

  init_lincomb_R_function = function(){

   self$control$lincomb_R_function <- function(x) x

  },

  # use a starter seed to create n_tree seeds
  plant_tree_seeds = function(start_seed){

   self$forest_seed <- start_seed
   set.seed(start_seed)
   self$tree_seeds <- sample(1e5, size = self$n_tree)

  },

  # computers

  compute_means = function(){

   numeric_data <- select_cols(self$data, private$data_names$x_numeric)

   if(self$na_action == 'omit'){
    numeric_data <- collapse::fsubset(numeric_data, private$data_rows_complete)
   }

   private$data_means <- collapse::fmean(numeric_data, w = self$weights)

  },


  compute_modes = function(){

   nominal_data <- select_cols(self$data, private$data_fctrs$cols)

   if(self$na_action == 'omit'){
    nominal_data <- collapse::fsubset(nominal_data, private$data_rows_complete)
   }

   private$data_modes <- vapply(nominal_data,
                                collapse::fmode,
                                FUN.VALUE = integer(1),
                                w = self$weights)

  },

  compute_stdev = function(){

   numeric_data <- select_cols(self$data, private$data_names$x_numeric)

   if(self$na_action == 'omit'){
    numeric_data <- collapse::fsubset(numeric_data, private$data_rows_complete)
   }

   private$data_stdev <- collapse::fsd(numeric_data, w = self$weights)

  },

  compute_bounds = function(){

   numeric_data <- select_cols(self$data, private$data_names$x_numeric)

   if(self$na_action == 'omit'){
    numeric_data <- collapse::fsubset(numeric_data, private$data_rows_complete)
   }

   private$data_bounds <- matrix(
    data = c(
     collapse::fnth(numeric_data, 0.10, w = self$weights),
     collapse::fnth(numeric_data, 0.25, w = self$weights),
     collapse::fnth(numeric_data, 0.50, w = self$weights),
     collapse::fnth(numeric_data, 0.75, w = self$weights),
     collapse::fnth(numeric_data, 0.90, w = self$weights)
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(c('10%', '25%', '50%', '75%', '90%'),
                    names(numeric_data))
   )

  },

  compute_mean_leaves = function(){

   leaf_counts <- vapply(X = self$forest$leaf_summary,
                         FUN = function(x) sum(x != 0),
                         FUN.VALUE = integer(1))

   private$mean_leaves <- collapse::fmean(leaf_counts)

  },

  compute_dependence_internal = function(pred_spec,
                                         pred_type,
                                         pred_horizon = NULL,
                                         type_input,
                                         type_output,
                                         prob_labels,
                                         prob_values,
                                         expand_grid,
                                         oobag){

   # make a visible binding for CRAN
   id_variable = NULL

   pred_horizon <- pred_horizon %||% self$pred_horizon %||% 1

   # just use first value for these pred_types
   if(pred_type %in% c("mort", "time")) pred_horizon <- pred_horizon[1]

   pred_horizon_order <- order(pred_horizon)
   pred_horizon_ordered <- pred_horizon[pred_horizon_order]


   private$init_data_rows_complete()
   private$prep_x()
   # y and w do not need to be prepped for prediction,
   # but they need to match orsf_cpp()'s expectations
   private$prep_y(placeholder = TRUE)
   private$w <- rep(1, nrow(private$x))

   if(oobag){ private$sort_inputs(sort_y = FALSE) }

   # the values in pred_spec need to be centered & scaled to match x,
   # which is also centered and scaled
   means <- private$data_means
   stdev <- private$data_stdev

   if(type_input == 'user'){

    for(i in intersect(names(means), names(pred_spec))){
     pred_spec[[i]] <- (pred_spec[[i]] - means[i]) / stdev[i]
    }

   } else {

    for(j in seq_along(pred_spec)){

     for(i in intersect(names(means), names(pred_spec[[j]]))){
      pred_spec[[j]][[i]] <- (pred_spec[[j]][[i]] - means[i]) / stdev[i]
     }

    }

    expand_grid <- FALSE

   }

   fi <- private$data_fctrs

   if(expand_grid){

    if(!is.data.frame(pred_spec))
     pred_spec <- expand.grid(pred_spec, stringsAsFactors = TRUE)

    for(i in seq_along(fi$cols)){

     ii <- fi$cols[i]

     if(is.character(pred_spec[[ii]]) && !fi$ordr[i]){

      pred_spec[[ii]] <- factor(pred_spec[[ii]], levels = fi$lvls[[ii]])

     }

    }

    check_new_data_fctrs(new_data  = pred_spec,
                         names_x   = private$data_names$x_original,
                         fi_ref    = fi,
                         label_new = "pred_spec")

    pred_spec_new <- ref_code(x_data = pred_spec, fi = fi,
                              names_x_data = names(pred_spec))

    x_cols <- list(match(names(pred_spec_new), colnames(private$x))-1)

    pred_spec_new <- list(as.matrix(pred_spec_new))

    pd_bind <- list(pred_spec)

   } else if(type_input == 'user'){

    pred_spec_new <- pd_bind <- x_cols <- list()

    for(i in seq_along(pred_spec)){

     pred_spec_new[[i]]  <- as.data.frame(pred_spec[i])
     pd_name <- names(pred_spec)[i]

     pd_bind[[i]] <- data.frame(
      variable = pd_name,
      value = rep(NA_real_, length(pred_spec[[i]])),
      level = rep(NA_character_, length(pred_spec[[i]]))
     )

     if(pd_name %in% fi$cols) {

      pd_bind[[i]]$level <- as.character(pred_spec[[i]])

      pred_spec_new[[i]] <- ref_code(pred_spec_new[[i]],
                                     fi = fi,
                                     names_x_data = pd_name)

     } else {

      pd_bind[[i]]$value <- pred_spec[[i]]

     }

     x_cols[[i]] <- match(names(pred_spec_new[[i]]), colnames(private$x)) - 1
     pred_spec_new[[i]] <- as.matrix(pred_spec_new[[i]])

    }

   } else {

    pd_bind <- pred_spec_new <- pred_spec

    for(i in seq_along(pred_spec_new)){
     pred_spec_new[[i]] <- ref_code(pred_spec_new[[i]], fi = fi,
                                    names_x_data = names(pred_spec_new[[i]]))
    }

    pred_spec_new <- lapply(pred_spec_new, collapse::qM)

    x_cols <- lapply(
     pred_spec_new,
     function(x){
      match(colnames(x), colnames(private$x)) - 1
     }
    )

   }

   cpp_args <- private$prep_cpp_args(x = private$x,
                                     y = private$y,
                                     w = private$w,
                                     importance_type = 'none',
                                     pred_type = pred_type,
                                     pred_mode = FALSE,
                                     pred_aggregate = TRUE,
                                     pred_horizon = pred_horizon_ordered,
                                     oobag = oobag,
                                     oobag_eval_type = 'none',
                                     pd_type_R = switch(type_output,
                                                        "smry" = 1L,
                                                        "ice" = 2L),
                                     pd_x_vals = pred_spec_new,
                                     pd_x_cols = x_cols,
                                     pd_probs = prob_values,
                                     write_forest = FALSE,
                                     run_forest = TRUE)

   pd_vals <- do.call(orsf_cpp, cpp_args)$pd_values

   row_delim <- switch(self$tree_type,
                       "survival" = pred_horizon_ordered,
                       "regression" = 1,
                       "classification" = self$class_levels)

   row_delim_label <- switch(self$tree_type,
                             "survival" = "pred_horizon",
                             "regression" = "pred_row",
                             "classification" = "class")

   for(i in seq_along(pd_vals)){

    pd_bind[[i]]$id_variable <- seq(nrow(pd_bind[[i]]))

    for(j in seq_along(pd_vals[[i]])){

     # nans <- which(is.nan(pd_vals[[i]][[j]]))
     nans <- which(pd_vals[[i]][[j]]==0)

     if(!is_empty(nans)){
      pd_vals[[i]][[j]][nans] <- NA_real_
     }

     pd_vals[[i]][[j]] <- matrix(pd_vals[[i]][[j]],
                                 nrow=length(row_delim),
                                 byrow = T)

     rownames(pd_vals[[i]][[j]]) <- row_delim


     if(type_output=='smry'){

      pd_vals[[i]][[j]] <- t(
       apply(
        X = pd_vals[[i]][[j]],
        MARGIN = 1,
        function(x, p){
         c(collapse::fmean(x), stats::quantile(x, p, na.rm=TRUE))
        },
        prob_values
       )
      )

      colnames(pd_vals[[i]][[j]]) <- c('mean', prob_labels)

     } else {

      colnames(pd_vals[[i]][[j]]) <- c(paste(1:nrow(private$x)))

     }

     pd_vals[[i]][[j]] <- as.data.table(pd_vals[[i]][[j]],
                                        keep.rownames = row_delim_label)

     if(type_output == 'ice'){

      measure.vars <- setdiff(names(pd_vals[[i]][[j]]), row_delim_label)

      pd_vals[[i]][[j]] <- melt_aorsf(data = pd_vals[[i]][[j]],
                                      id.vars = row_delim_label,
                                      variable.name = 'id_row',
                                      value.name = 'pred',
                                      measure.vars = measure.vars)

     }

    }

    pd_vals[[i]] <- rbindlist(pd_vals[[i]], idcol = 'id_variable')

    # this seems awkward but the reason I convert back to data.frame
    # here is to avoid a potential memory leak from forder & bmerge.
    # I have no idea why this memory leak may be occurring but it does
    # not if I apply merge.data.frame instead of merge.data.table
    pd_vals[[i]] <- merge(as.data.frame(pd_vals[[i]]),
                          as.data.frame(pd_bind[[i]]),
                          by = 'id_variable')

    if(type_input == 'intr'){

     v1 <- colnames(pred_spec[[i]])[1]
     v2 <- colnames(pred_spec[[i]])[2]

     pd_vals[[i]][['var_1_name']] <- v1
     pd_vals[[i]][['var_2_name']] <- v2

     names(pd_vals[[i]])[names(pd_vals[[i]]) == v1] <- "var_1_value"
     names(pd_vals[[i]])[names(pd_vals[[i]]) == v2] <- "var_2_value"

     pd_vals[[i]]$var_1_value <- as.numeric(pd_vals[[i]]$var_1_value)
     pd_vals[[i]]$var_2_value <- as.numeric(pd_vals[[i]]$var_2_value)

    }

   }

   out <- rbindlist(pd_vals)

   # missings may occur when oobag=TRUE and n_tree is small
   if(type_output == 'ice') {
    out <- collapse::na_omit(out, cols = 'pred')
   }

   ids <- c('id_variable')

   if(type_output == 'ice') ids <- c(ids, 'id_row')

   mid <- setdiff(names(out), c(ids, 'mean', prob_labels, 'pred'))

   end <- setdiff(names(out), c(ids, mid))

   setcolorder(out, neworder = c(ids, mid, end))

   if(self$tree_type == 'classification'){
    out[, class := factor(class, levels = self$class_levels)]
    setkey(out, class)
   }


   if(self$tree_type == 'survival' && !(pred_type %in% c('mort', 'time')))
    out[, pred_horizon := as.numeric(pred_horizon)]

   if(self$tree_type == 'regression'){
    out[, pred_row := NULL]
   }

   if(pred_type %in% c('mort', 'time'))
    out[, pred_horizon := NULL]

   # not needed for summary
   if(type_output == 'smry')
    out[, id_variable := NULL]

   # put data back into original scale
   if(type_input == 'intr'){

    for(j in collapse::funique(out$var_1_name)){

     if(j %in% names(means)){
      var_index <- out$var_1_name %==% j
      var_value <- (out$var_1_value[var_index] * stdev[j]) + means[j]
      set(out, i = var_index, j = 'var_1_value', value = var_value)
     }


    }

    for(j in collapse::funique(out$var_2_name)){

     if(j %in% names(means)){
      var_index <- out$var_2_name %==% j
      var_value <- (out$var_2_value[var_index] * stdev[j]) + means[j]
      set(out, i = var_index, j = 'var_2_value', value = var_value)
     }

    }

   } else {

    for(j in intersect(names(means), names(pred_spec))){

     if(j %in% names(out)){

      var_index <- collapse::seq_row(out)
      var_value <- (out[[j]] * stdev[j]) + means[j]
      var_name  <- j

     } else {

      var_index <- out$variable %==% j
      var_value <- (out$value[var_index] * stdev[j]) + means[j]
      var_name  <- 'value'

     }

     set(out, i = var_index, j = var_name, value = var_value)

    }

   }


   # silent print after modify in place
   out[]

   out

  },

  select_variables_internal = function(n_predictor_min, verbose_progress){

   n_predictors <- length(private$data_names$x_original)

   # verbose progress on the forest should always be FALSE
   # because for orsf_vs, verbosity is coordinated in R
   self$verbose_progress <- FALSE

   oob_data <- data.table(
    n_predictors = seq(n_predictors),
    stat_value = rep(NA_real_, n_predictors),
    predictors_included = vector(mode = 'list', length = n_predictors),
    predictor_dropped = rep(NA_character_, n_predictors)
   )

   # if the forest was not trained prior to variable selection
   if(!self$trained){
    private$compute_means()
    private$compute_stdev()
    private$compute_modes()
    private$compute_bounds()
   }

   private$prep_x()
   private$prep_y()
   private$prep_w()

   # for survival, inputs should be sorted by time
   private$sort_inputs()

   # allow re-training.
   self$forest <- list()

   pred_type <- switch(self$tree_type,
                       'survival' = 'mort',
                       'classification' = 'prob',
                       'regression' = 'mean')

   cpp_args <- private$prep_cpp_args(pred_type = pred_type,
                                     oobag_pred = TRUE,
                                     importance_group_factors = TRUE,
                                     write_forest = FALSE)

   max_progress <- n_predictors - n_predictor_min
   current_progress <- 0
   start_time <- last_time <- Sys.time()

   while(n_predictors >= n_predictor_min){

    if(verbose_progress){

     time_passed <- as.numeric(
      as.difftime(Sys.time() - start_time),
      units = "secs"
     )

     time_lapsed <- as.numeric(
      as.difftime(Sys.time() - last_time),
      units = "secs"
     )

     if(current_progress > 0 && time_lapsed > 3){

      relative_progress <- current_progress / max_progress

      remaining_time = round((1 / relative_progress - 1) * time_passed)

      cat(
       "Selecting variables:",
       paste0( round_magnitude(relative_progress*100), "%" ),
       "~ time remaining:", beautifyTime(remaining_time),
       "\n")

      last_time <- Sys.time()

     }

    }

    mtry_safe <- ceiling(sqrt(n_predictors))

    if(self$control$lincomb_df_target > mtry_safe){
     self$control$lincomb_df_target <- mtry_safe
    }

    cpp_args$mtry <- mtry_safe
    cpp_output <- do.call(orsf_cpp, args = cpp_args)

    worst_index <- which.min(cpp_output$importance)
    worst_predictor <- colnames(cpp_args$x)[worst_index]

    oob_data[n_predictors,
             `:=`(n_predictors = n_predictors,
                  stat_value = cpp_output$eval_oobag$stat_values[1,1],
                  predictors_included = colnames(cpp_args$x),
                  predictor_dropped = worst_predictor)]

    cpp_args$x <- cpp_args$x[, -worst_index]
    n_predictors <- n_predictors - 1
    current_progress <- current_progress + 1

   }

   if(verbose_progress){
    cat("Selecting variables: 100%\n")
   }

   collapse::na_omit(oob_data)

  },

  # preppers

  prep_x = function(){

   data_x <- select_cols(self$data, private$data_names$x_original)

   if(any(is.na(data_x))){

    switch(
     self$na_action,

     'fail' = {
      # this should already have been checked but just in case
      stop("Please remove missing values from data, or impute them.",
           call. = FALSE)
     },

     'omit' = {
      data_x <- data_x[private$data_rows_complete, ]
     },

     'pass' = {
      data_x <- data_x[private$data_rows_complete, ]
     },

     'impute_meanmode' = {

      data_x <- data_impute(data_x,
                            cols = names(data_x),
                            values = c(as.list(private$data_means),
                                       as.list(private$data_modes)))
     }
    )

   }

   x <- ref_code(data_x, private$data_fctrs, names(data_x))

   for(i in private$data_names$x_numeric){

    if(has_units(x[[i]])) x[[i]] <- as.numeric(x[[i]])

    # can't modify by reference here, it would modify the user's data
    x[[i]] <- (x[[i]] - private$data_means[i]) / private$data_stdev[i]

   }

   private$x = as_matrix(x)

  },

  # Prep the outcome variable
  #
  # Coerce the outcome to be compatible with C++ routines.
  # If there is no feasible way to make it work, throw an error.
  # If there aren't any problems, return the outcome.
  #
  # placeholder (*logical*) whether to engage with the
  #   outcome member variable or just make a version that is safe
  #   to pass to orsf_cpp().
  #
  # For survival, y is a matrix with two columns:
  #   - first column: time values
  #   - second column: status values
  #
  # For classification, y is a factor or character vector
  #
  # For regression, y is a numeric vector
  #
  # @return stored the outcome in a matrix format in private$y
  #
  # - Survival outcomes are transformed to have status values of 0 and 1
  # - Classification outcomes are transformed to a one-hot matrix
  # - Regression outcomes are converted to a matrix
  #
  prep_y = function(placeholder = FALSE){

   private$y <- select_cols(self$data, private$data_names$y)

   if(self$na_action == 'omit' && !placeholder)
    private$y <- private$y[private$data_rows_complete, ]

   private$prep_y_internal(placeholder)

  },

  prep_w = function(){

   # re-scale so that sum(w) == nrow(data)
   private$w <- self$weights * length(self$weights) / sum(self$weights)

  },

  prep_cpp_args = function(...){

   .dots <- list(...)

   args <- list(
    x = private$x,
    y = private$y,
    w = private$w,
    tree_type_R = switch(self$tree_type,
                         'classification' = 1,
                         'regression'= 2,
                         'survival' = 3),
    tree_seeds = self$tree_seeds,
    loaded_forest = self$forest,
    n_tree = .dots$n_tree %||% self$n_tree,
    mtry = .dots$mtry %||% self$mtry,
    sample_with_replacement = .dots$sample_with_replacement %||% self$sample_with_replacement,
    sample_fraction = .dots$sample_fraction %||% self$sample_fraction,
    vi_type_R = switch(.dots$importance_type %||% self$importance_type,
                       "none"    = 0,
                       "negate"  = 1,
                       "permute" = 2,
                       "anova"   = 3),
    vi_max_pvalue = .dots$importance_max_pvalue %||% self$importance_max_pvalue,
    leaf_min_events = .dots$leaf_min_events %||% self$leaf_min_events %||% 1,
    leaf_min_obs = .dots$leaf_min_obs %||% self$leaf_min_obs,
    split_rule_R = switch(self$split_rule,
                          "logrank"  = 1,
                          "cstat"    = 2,
                          "gini"     = 3,
                          "variance" = 4),
    split_min_events = .dots$split_min_events %||% self$split_min_events %||% 1,
    split_min_obs = .dots$split_min_obs %||% self$split_min_obs,
    split_min_stat = .dots$split_min_stat %||% self$split_min_stat,
    split_max_cuts = .dots$split_max_cuts %||% self$n_split,
    split_max_retry = .dots$split_max_retry %||% self$n_retry,
    lincomb_R_function = self$control$lincomb_R_function,
    lincomb_type_R = switch(self$control$lincomb_type,
                            'glm'    = 1,
                            'random' = 2,
                            'net'    = 3,
                            'custom' = 4),
    lincomb_eps = self$control$lincomb_eps,
    lincomb_iter_max = self$control$lincomb_iter_max,
    lincomb_scale = self$control$lincomb_scale,
    lincomb_alpha = self$control$lincomb_alpha,
    lincomb_df_target = self$control$lincomb_df_target,
    lincomb_ties_method = switch(tolower(self$control$lincomb_ties_method),
                                 'breslow' = 0,
                                 'efron'   = 1),
    pred_type_R = switch(.dots$pred_type %||% self$pred_type,
                         "none"  = 0,
                         "risk"  = 1,
                         "surv"  = 2,
                         "chf"   = 3,
                         "mort"  = 4,
                         "mean"  = 5,
                         "prob"  = 6,
                         "class" = 7,
                         "leaf"  = 8,
                         "time"  = 2), # time=2 is not a typo
    pred_mode = .dots$pred_mode %||% FALSE,
    pred_aggregate = .dots$pred_aggregate %||% (self$pred_type != 'leaf'),
    pred_horizon = if(self$pred_type == 'time'){
     self$event_times
    } else {
     .dots$pred_horizon %||% self$pred_horizon %||% 1
    },
    oobag = .dots$oobag %||% self$oobag_pred_mode,
    oobag_R_function = .dots$oobag_eval_function %||% self$oobag_eval_function,
    oobag_eval_type_R = switch(
     tolower(.dots$oobag_eval_type %||% self$oobag_eval_type),
     "none" = 0,
     "harrell's c-index" = 1,
     "auc-roc" = 1,
     "user-specified function" = 2,
     "mse" = 3,
     "rsq" = 4
    ),
    oobag_eval_every = .dots$oobag_eval_every %||% self$oobag_eval_every,
    # switch(pd_type, "none" = 0L, "smry" = 1L, "ice" = 2L)
    pd_type_R = .dots$pd_type %||% 0L,
    pd_x_vals = .dots$pd_x_vals %||% list(matrix(0, ncol=0, nrow=0)),
    pd_x_cols = .dots$pd_x_cols %||% list(matrix(0, ncol=0, nrow=0)),
    pd_probs  = .dots$pd_probs %||% 0,
    n_thread  = .dots$n_thread %||% self$n_thread,
    write_forest = .dots$write_forest %||% TRUE,
    run_forest   = .dots$run_forest %||% TRUE,
    verbosity    = as.integer(.dots$verbosity %||% self$verbose_progress)
   )

  },

  sort_inputs = function(sort_y = NULL,
                         sort_x = NULL,
                         sort_w = NULL){
   NULL
  },

  # cleaners

  clean_importance = function(importance = NULL, group_factors = NULL){

   out <- importance %||% self$importance
   group_factors <- group_factors %||% self$importance_group_factors

   # nan indicates a variable was never used
   out[is.nan(out)] <- 0

   rownames(out) <- private$data_names$x_ref_code

   private$importance_raw <- out

   if(group_factors){

    fi <- private$data_fctrs

    if(!is_empty(fi$cols)){

     for(f in fi$cols[!fi$ordr]){

      f_lvls <- fi$lvls[[f]]
      f_rows <- match(paste(f, f_lvls[-1], sep = '_'), rownames(out))
      f_wts <- 1

      if(length(f_lvls) > 2){
       f_wts <- prop.table(x = table(self$data[[f]], useNA = 'no')[-1])
      }

      f_vi <- sum(out[f_rows] * f_wts, na.rm = TRUE)

      out[f_rows] <- f_vi
      rownames(out)[f_rows] <- f

     }

     if(!is_empty(fi$cols[!fi$ordr])) {
      # take extra care in case there are duplicate vi values.
      out <- data.frame(variable = rownames(out), value = out)
      out <- unique(out)
      out$variable <- NULL
      out <- as.matrix(out)
      colnames(out) <- NULL
     }

    }

   }

   self$importance <- rev(out[order(out), , drop=TRUE])

  },

  clean_pred_oobag = function(){

   if(self$pred_type == 'leaf'){

    all_rows <- seq(self$n_obs)

    for(i in seq(self$n_tree)){

     rows_inbag <- setdiff(all_rows, self$forest$rows_oobag[[i]]+1)
     self$pred_oobag[rows_inbag, i] <- NA_real_

    }

   }

   self$pred_oobag[is.nan(self$pred_oobag)] <- NA_real_

   private$clean_pred_oobag_internal()

  },

  clean_pred_oobag_internal = function(){
   NULL
  },

  clean_pred_new = function(preds){

   if(self$na_action == "pass"){

    out <- matrix(nrow = nrow(self$data),
                  ncol = ncol(preds))

    out[private$data_rows_complete, ] <- preds

   } else {

    out <- preds

   }

   if(self$pred_type == "leaf" || !self$pred_aggregate) return(out)

   private$clean_pred_new_internal(out)

  },

  restore_state = function(public_state = list(),
                           private_state = list()){

   if(!is_empty(public_state)){
    for(i in names(public_state)){
     self[[i]] <- public_state[[i]]
    }
   }

   if(!is_empty(private_state)){
    for(i in names(private_state)){
     private[[i]] <- private_state[[i]]
    }
   }

  }


 )

)

# ObliqueForestSurvival class ----

ObliqueForestSurvival <- R6::R6Class(
 "ObliqueForestSurvival",
 inherit = ObliqueForest,
 public = list(

  get_max_time = function(){
   return(private$max_time)
  },

  get_pred_type_vi = function(){
   return("mort")
  },

  leaf_min_events = NULL,
  split_min_events = NULL,
  pred_horizon = NULL,
  event_times = NULL

 ),

 # private ----
 private = list(

  pred_horizon_order = NULL,
  data_row_sort = NULL,
  max_time = NULL,
  n_events = NULL,

  check_split_rule_internal= function(){

   check_arg_is_valid(arg_value = self$split_rule,
                      arg_name = 'split_rule',
                      valid_options = c("logrank", "cstat"))

  },
  check_pred_type_internal = function(oobag,
                                      pred_type = NULL,
                                      context = NULL){

   input <- pred_type %||% self$pred_type

   arg_name <- if(oobag) 'oobag_pred_type' else 'pred_type'

   if(is.null(context)){
    valid_options <- c("none", "surv", "risk", "chf", "mort", "leaf", "time")
   } else {
    valid_options <- switch(
     context,
     'partial dependence' = c("surv", "risk", "chf", "mort", "time"),
     'prediction' = c("surv", "risk", "chf", "mort", "leaf", "time")
    )
    context <- paste(context, 'with survival forests')
   }

   check_arg_is_valid(arg_value = input,
                      arg_name = arg_name,
                      valid_options = valid_options,
                      context = context)

  },

  check_pred_horizon = function(pred_horizon = NULL,
                                boundary_checks = TRUE,
                                pred_type = NULL){

   pred_type <- pred_type %||% self$pred_type
   input <- pred_horizon %||% self$pred_horizon

   if(is.null(input) && pred_type %in% c('risk', 'surv', 'chf')){

    stop("pred_horizon must be specified for ",
         pred_type, " predictions.", call. = FALSE)

   }


   if(self$oobag_pred_mode)
    arg_name <- 'oobag_pred_horizon'
   else
    arg_name <- 'pred_horizon'

   check_arg_type(arg_value = input,
                  arg_name = 'pred_horizon',
                  expected_type = 'numeric')

   for(i in seq_along(input)){

    check_arg_gteq(arg_value = input[i],
                   arg_name = 'pred_horizon',
                   bound = 0)

   }

   if(boundary_checks){

    if(any(input > private$max_time)){

     stop("prediction horizon should ",
          "be <= max follow-up time ",
          "observed in training data: ",
          private$max_time,
          call. = FALSE)

    }

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

  check_oobag_eval_function_internal = function(oobag_fun){

   test_time <- seq(from = 1, to = 5, length.out = 100)
   test_status <- rep(c(0,1), each = 50)

   .y_mat <- cbind(time = test_time, status = test_status)
   .w_vec <- rep(1, times = 100)
   .s_vec <- seq(0.9, 0.1, length.out = 100)

   test_output <- try(oobag_fun(y_mat = .y_mat,
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

   test_output

  },

  check_lincomb_R_function_internal = function(lincomb_R_function = NULL){

   input <- lincomb_R_function %||% self$lincomb_R_function

   test_time <- seq(from = 1, to = 5, length.out = 100)
   test_status <- rep(c(0,1), each = 50)

   .x_node <- matrix(rnorm(300), ncol = 3)
   .y_node <- cbind(time = test_time, status = test_status)
   .w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)


   out <- try(input(.x_node, .y_node, .w_node), silent = FALSE)

   if(is_error(out)){

    stop("user-supplied function encountered an error when it was tested. ",
         "Please make sure the function works for this case:\n\n",
         "test_time <- seq(from = 1, to = 5, length.out = 100)\n",
         "test_status <- rep(c(0,1), each = 50)\n\n",
         ".x_node <- matrix(seq(-1, 1, length.out = 300), ncol = 3)\n",
         ".y_node <- cbind(time = test_time, status = test_status)\n",
         ".w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)\n\n",
         "test_output <- user_function(.x_node, .y_node, .w_node)\n\n",
         "test_output should be a numeric matrix with 1 column and",
         " with nrow(test_output) = ncol(.x_node)",
         call. = FALSE)

   }

   out

  },

  sort_inputs = function(sort_x = TRUE,
                         sort_y = TRUE,
                         sort_w = TRUE){

   if(sort_x)
    private$x <- private$x[private$data_row_sort, , drop = FALSE]
   if(sort_y)
    private$y <- private$y[private$data_row_sort, , drop = FALSE]
   if(sort_w)
    private$w <- private$w[private$data_row_sort]

  },

  init_control = function(){

   self$control <- orsf_control_survival(method = 'glm',
                                         scale_x = FALSE,
                                         max_iter = 1)

  },

  init_pred_type = function(){
   self$pred_type <- 'risk'
  },

  init_split_rule = function(){
   self$split_rule <- 'logrank'
  },

  init_split_min_stat = function(){

   if(is.null(self$split_rule))
    stop("cannot init split_min_stat without split_rule", call. = FALSE)

   self$split_min_stat <- switch(self$split_rule,
                                 'logrank' = 3.841459,
                                 'cstat' = 0.55)
  },

  init_internal = function(){

   self$tree_type <- "survival"

   if(!is.function(self$control$lincomb_R_function) &&
      self$control$lincomb_type == 'net'){
    self$control$lincomb_R_function <- penalized_cph
   }

   y <- select_cols(self$data, private$data_names$y)

   if(inherits(y[[1]], 'Surv')){
    y <- as.matrix(y[[1]])
    y <- as.data.frame(y)
   }

   # strip units for these tests
   if(has_units(y[[1]])) y[[1]] <- as.numeric(y[[1]])

   check_arg_type(arg_value = y[[1]],
                  arg_name = "time to event",
                  expected_type = 'numeric')

   check_arg_gt(arg_value = y[[1]],
                arg_name = "time to event",
                bound = 0)

   # strip extra attributes for these tests
   if(is.factor(y[[2]])  ||
      is.logical(y[[2]]) ||
      has_units(y[[2]])  ){
    y[[2]] <- as.numeric(y[[2]])
   }

   check_arg_type(arg_value = y[[2]],
                  arg_name = "status indicator",
                  expected_type = 'numeric')

   if(self$na_action == 'omit')
    y <- y[private$data_rows_complete, ]


   # order observations by event time and event status
   private$data_row_sort <- collapse::radixorder(y[[1]], -y[[2]])
   # boundary check for pred horizon
   private$max_time <- y[last_value(private$data_row_sort), 1]
   # boundary check for event-based tree parameters
   private$n_events <- collapse::fsum(y[, 2])
   # unique event times sorted in ascending order
   self$event_times <- sort(collapse::funique(
    y[[1]][collapse::whichv(x = y[[2]], value = 1)]
   ), decreasing = FALSE)

   # if pred_horizon is unspecified, provide sensible default
   # if it is specified, check for correctness
   if(private$user_specified$pred_horizon){
    private$check_pred_horizon(self$pred_horizon, boundary_checks = TRUE)
   } else {
    self$pred_horizon <- collapse::fmedian(y[, 1])
   }

   private$check_leaf_min_events()
   private$check_split_min_events()

   if(!self$oobag_pred_mode) self$oobag_eval_type <- "none"

   # use default if eval type was not specified by user
   if(self$oobag_pred_mode && is.null(self$oobag_eval_type)){
    self$oobag_eval_type <- "Harrell's C-index"
   }

  },

  init_pred = function(pred_horizon = NULL, pred_type = NULL,
                       boundary_checks = TRUE){

   pred_type_supplied <- !is.null(pred_type)
   pred_horizon_supplied <- !is.null(pred_horizon)

   if(pred_type_supplied){
    private$check_pred_type(oobag = FALSE, pred_type = pred_type)
   } else {
    pred_type <- self$pred_type %||% "risk"
   }

   if(pred_horizon_supplied){
    private$check_pred_horizon(pred_horizon, boundary_checks, pred_type)
   } else {
    pred_horizon <- self$pred_horizon
    if(is.null(pred_horizon)){
     stop("pred_horizon was not specified and could not be found in object.",
          call. = FALSE)
    }
   }

   if(pred_type_supplied &&
      pred_horizon_supplied &&
      pred_type %in% c('leaf', 'mort', 'time')){

    extra_text <- if(length(pred_horizon)>1){
     " Predictions at each value of pred_horizon will be identical."
    } else {
     ""
    }

    warning("pred_horizon does not impact predictions",
            " when pred_type is '", pred_type, "'.",
            extra_text, call. = FALSE)

   }

   self$pred_horizon <- pred_horizon
   self$pred_type <- pred_type

  },

  prep_y_internal = function(placeholder = FALSE){


   if(placeholder){
    private$y <- matrix(0, ncol = 2, nrow = 1)
    return()
   }

   y <- private$y
   cols <- names(y)

   # if a surv object was included in the formula, it probably
   # doesn't need to be checked here, but it's still checked
   # just to be safe b/c if there is a problem with y it is
   # likely that it will cause R to crash in orsf_cpp().
   if(length(cols) == 1 && inherits(y[[1]], 'Surv')){
    y <- as.data.frame(as.matrix(y))
    cols <- names(y)
   }

   for(i in cols){
    if(has_units(y[[i]])) y[[i]] <- as.numeric(y[[i]])
   }

   # make sure y gets passed to CPP with status values of 0 and 1
   if(is.factor(y[[2]]) || is.logical(y[[2]])){
    y[[2]] <- as.numeric(y[[2]])
   }

   status_uni <- collapse::funique(y[[2]])

   # if nothing is censored
   if(all(status_uni == 1)) {

    private$y <- as_matrix(y)
    return()

   }

   # status values are modified if they are not all 0 and 1
   if(!is_equivalent(c(0, 1), status_uni)){

    # assume the lowest value of status indicates censoring
    censor_indicator <- collapse::fmin(status_uni)

    if(censor_indicator != as.integer(censor_indicator)){
     stop("the censoring indicator is not integer valued. ",
          "This can occur when the time column is swapped ",
          "with the status column by mistake. ",
          "Did you enter `status + time` when you meant ",
          "to put `time + status`?", call. = FALSE)
    }

    # assume the integer value 1 above censor is an event
    event_indicator <- censor_indicator + 1

    if( !(event_indicator %in% status_uni) ){
     stop("there does not appear to be an event indicator ",
          "in the status column. The censoring indicator appears to ",
          "be ", censor_indicator, " but there are no values of ",
          1 + censor_indicator, ", which we would assume to be an ",
          " event indicator. This can occur when the time column is ",
          "swapped with the status column by mistake. Did you enter ",
          "`status + time` when you meant  to enter `time + status`? ",
          "If this problem persists, try setting status to 0 for ",
          "censored observations and 1 for events.", call. = FALSE)
    }

    other_events <- which(y[[2]] > event_indicator)

    # we don't handle competing risks yet
    if(!is_empty(other_events)){

     stop("detected >1 event type in status variable. ",
          "Currently only 1 event type is supported. ",
          call. = FALSE)

     # y[other_events, 2] <- censor_indicator

    }

    # The status column needs to contain only 0s and 1s when it
    # gets passed into the C++ routines.
    y[[2]] <- y[[2]] - censor_indicator

    # check to make sure we don't send something into C++
    # that is going to make the R session crash.
    status_uni <- collapse::funique(y[[2]])

    if(!is_equivalent(c(0, 1), status_uni)){
     stop("could not coerce the status column to values of 0 and 1. ",
          "This can occur when the time column is ",
          "swapped with the status column by mistake. Did you enter ",
          "`status + time` when you meant  to enter `time + status`? ",
          "Please modify your data so that the status values are ",
          "0 and 1, with 0 indicating censored observations and 1 ",
          "indicating events.", call. = FALSE)
    }

   }

   private$y <- as_matrix(y)

  },

  clean_pred_oobag_internal = function(){


   # put the oob predictions into the same order as the training data.
   unsorted <- collapse::radixorder(private$data_row_sort)
   self$pred_oobag <- self$pred_oobag[unsorted, , drop = FALSE]

   if(self$pred_type == 'time'){

    self$eval_oobag$stat_values <-
     self$eval_oobag$stat_values[, ncol(self$pred_oobag), drop = FALSE]

    self$pred_oobag <- apply(self$pred_oobag,
                             MARGIN = 1,
                             FUN = rmst,
                             times = self$event_times)

    self$pred_oobag <- collapse::qM(self$pred_oobag)

   }

   # these predictions should always be 1 column
   # b/c they do not depend on the prediction horizon
   if(self$pred_type == 'mort'){

    self$eval_oobag$stat_values <-
     self$eval_oobag$stat_values[, 1L, drop = FALSE]

    self$pred_oobag <- self$pred_oobag[, 1L, drop = FALSE]

   }

  },
  clean_pred_new_internal = function(preds){

   if(self$pred_type == 'time'){

    preds <- apply(preds,
                   MARGIN = 1,
                   FUN = rmst,
                   times = self$event_times)
    preds <- collapse::qM(preds)

   }

   # don't let multiple pred horizon values through for mort
   if(self$pred_type %in% c('mort', 'time')){
    return(preds[, 1, drop = FALSE])
   }

   # output in the same order as user's pred_horizon vector
   preds <- preds[, order(private$pred_horizon_order), drop = FALSE]

   preds

  },

  predict_internal = function(simplify){

   private$pred_horizon_order <- order(self$pred_horizon)
   pred_horizon_ordered <- self$pred_horizon[private$pred_horizon_order]

   cpp_args = private$prep_cpp_args(x = private$x,
                                    y = private$y,
                                    w = private$w,
                                    importance_type = 'none',
                                    pred_type = self$pred_type,
                                    pred_aggregate = self$pred_aggregate,
                                    pred_horizon = pred_horizon_ordered,
                                    oobag_pred = FALSE,
                                    pred_mode = TRUE,
                                    write_forest = FALSE,
                                    run_forest = TRUE)

   if(length(self$pred_horizon) > 1 && !self$pred_aggregate){

    results <- vector(mode = 'list', length = length(self$pred_horizon))

    for(i in seq_along(results)){

     cpp_args$pred_horizon <- self$pred_horizon[i]

     results[[i]] <- do.call(orsf_cpp, args = cpp_args)$pred_new

     results[[i]] <- private$clean_pred_new(results[[i]])

    }

    # all components are the same if pred type is mort
    # (user also gets a warning if they ask for this)
    if(self$pred_type %in% c('mort', 'leaf', 'time')) return(results[[1]])

    if(simplify){
     results <- simplify2array(results)
    }

    return(results)

   }

   out_values <- do.call(orsf_cpp, args = cpp_args)$pred_new

   private$clean_pred_new(out_values)

  }

 )

)

# ObliqueForestClassification class ----

ObliqueForestClassification <- R6::R6Class(
 "ObliqueForestClassification",
 inherit = ObliqueForest,
 public = list(

  n_class = NULL,

  class_levels = NULL,

  get_pred_type_vi = function(){
   return("prob")
  }

 ),

 # private ----
 private = list(

  check_split_rule_internal = function(){

   check_arg_is_valid(arg_value = self$split_rule,
                      arg_name = 'split_rule',
                      valid_options = c("gini", "cstat"))

  },
  check_pred_type_internal = function(oobag,
                                      pred_type = NULL,
                                      context = NULL){

   input <- pred_type %||% self$pred_type

   arg_name <- if(oobag) 'oobag_pred_type' else 'pred_type'

   if(is.null(context)){
    valid_options <- c("none", "prob", "class", "leaf")
   } else {
    valid_options <- switch(
     context,
     'partial dependence' = c("prob"),
     'prediction' = c("prob", "class", "leaf")
    )
    context <- paste(context, 'with classification forests')
   }

   check_arg_is_valid(arg_value = input,
                      arg_name = arg_name,
                      valid_options = valid_options,
                      context = context)

  },

  check_pred_horizon = function(pred_horizon = NULL,
                                boundary_checks = TRUE,
                                pred_type = NULL){

   # nothing to check
   NULL

  },

  check_oobag_eval_function_internal = function(oobag_fun){

   test_y <- rep(c(0,1), each = 50)

   .y_mat <- matrix(test_y, ncol = 1)
   .w_vec <- rep(1, times = 100)
   .s_vec <- seq(0.9, 0.1, length.out = 100)

   test_output <- try(oobag_fun(y_mat = .y_mat,
                                w_vec = .w_vec,
                                s_vec = .s_vec),
                      silent = FALSE)

   if(is_error(test_output)){

    stop("oobag_fun encountered an error when it was tested. ",
         "Please make sure your oobag_fun works for this case:\n\n",
         "test_y <- rep(c(0,1), each = 50)\n",
         "y_mat <- matrix(test_y, ncol = 1)\n",
         "w_vec <- rep(1, times = 100)\n",
         "s_vec <- seq(0.9, 0.1, length.out = 100)\n\n",
         "test_output <- oobag_fun(y_mat = y_mat, w_vec = w_vec, s_vec = s_vec)\n\n",
         "test_output should be a numeric value of length 1",
         call. = FALSE)

   }

   test_output

  },

  check_lincomb_R_function_internal = function(lincomb_R_function = NULL){

   input <- lincomb_R_function %||% self$lincomb_R_function

   .x_node <- matrix(rnorm(300), ncol = 3)
   .y_node <- matrix(rbinom(100, size = 1, prob = 1/2), ncol = 1)
   .w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)

   out <- try(input(.x_node, .y_node, .w_node), silent = FALSE)

   if(is_error(out)){

    stop("user-supplied function encountered an error when it was tested. ",
         "Please make sure the function works for this case:\n\n",
         ".x_node <- matrix(rnorm(300), ncol = 3)\n",
         ".y_node <- matrix(rbinom(100, size = 1, prob = 1/2), ncol = 1)\n",
         ".w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)\n",
         "test_output <- your_function(.x_node, .y_node, .w_node)\n\n",
         "test_output should be a numeric matrix with 1 column and",
         " with nrow(test_output) = ncol(.x_node)",
         call. = FALSE)

   }

   out

  },

  init_control = function(){

   self$control <- orsf_control_classification(method = 'glm',
                                               scale_x = FALSE,
                                               max_iter = 1)

  },

  init_pred_type = function(){
   self$pred_type <- 'prob'
  },

  init_split_rule = function(){
   self$split_rule <- 'gini'
  },

  init_split_min_stat = function(){

   if(is.null(self$split_rule))
    stop("cannot init split_min_stat without split_rule", call. = FALSE)

   self$split_min_stat <- switch(self$split_rule,
                                 'gini' = 0,
                                 'cstat' = 0.55)

  },

  init_internal = function(){

   self$tree_type <- "classification"

   if(!is.function(self$control$lincomb_R_function) &&
      self$control$lincomb_type == 'net'){
    self$control$lincomb_R_function <- penalized_logreg
   }

   if(!self$oobag_pred_mode) self$oobag_eval_type <- "none"

   # use default if eval type was not specified by user
   if(self$oobag_pred_mode && is.null(self$oobag_eval_type)){
    self$oobag_eval_type <- "AUC-ROC"
   }

   y <- self$data[[private$data_names$y]]

   if(is.factor(y)){
    self$class_levels <- levels(y)
    self$n_class <- length(self$class_levels)
   } else {
    self$class_levels <- unique(y)
    self$n_class <- length(self$class_levels)
   }



  },

  init_pred = function(pred_horizon = NULL, pred_type = NULL,
                       boundary_checks = TRUE){

   if(!is.null(pred_horizon)){
    warning("pred_horizon does not impact predictions",
            " for classification forests", call. = FALSE)
   }

   if(!is.null(pred_type)){
    private$check_pred_type(oobag = FALSE, pred_type = pred_type)
   } else {
    pred_type <- self$pred_type %||% "prob"
   }

   self$pred_type <- pred_type

  },

  prep_y_internal = function(placeholder = FALSE){

   if(placeholder){
    private$y <- matrix(0, ncol = self$n_class, nrow = 1)
    return()
   }

   # y is always 1 column for classification (right?)
   y <- private$y[[1]]

   input_was_numeric <- !is.factor(y)

   if(input_was_numeric) y <- as.factor(y)

   n_class <- length(levels(y))

   if(n_class > 5 && input_was_numeric){
    stop("The outcome is numeric and has > 5 unique values.",
         " Did you mean to use orsf_control_regression()? If not,",
         " please convert ", private$data_names$y, " to a factor and re-run",
         call. = FALSE)
   }

   y <- as.numeric(y) - 1

   if(min(y) > 0) stop("y is less than 0")

   private$y <- expand_y_clsf(as_matrix(y), n_class)

  },

  predict_internal = function(simplify){

   # resize y to have the right number of columns
   private$y <- matrix(0, ncol = self$n_class)

   cpp_args = private$prep_cpp_args(x = private$x,
                                    y = private$y,
                                    w = private$w,
                                    importance_type = 'none',
                                    pred_type = self$pred_type,
                                    pred_aggregate = self$pred_aggregate,
                                    oobag_pred = FALSE,
                                    pred_mode = TRUE,
                                    write_forest = FALSE,
                                    run_forest = TRUE)

   if(!self$pred_aggregate && self$n_class > 2 && self$pred_type == 'prob'){
    stop("unaggregated probability predictions for outcomes with >2 classes",
         " is not currently supported.", call. = FALSE)
   }

   out <- do.call(orsf_cpp, args = cpp_args)$pred_new

   if(self$pred_type == 'prob' && self$pred_aggregate){

    colnames(out) <- self$class_levels

    if(simplify && self$n_class == 2)
     out <- out[, 2L, drop = TRUE]

   }

   if(self$pred_type == 'class'){

    # cpp class levels start at 0, R levels start at 1
    out <- out + 1

    # convert to factor (and coerce to vector) if asked
    if(simplify)
     out <- factor(out,
                   levels = seq(self$n_class),
                   labels = self$class_levels)
   }

   out

  }

 )
)


# ObliqueForestRegression class ----

ObliqueForestRegression <- R6::R6Class(
 "ObliqueForestRegression",
 inherit = ObliqueForest,
 public = list(

  get_pred_type_vi = function(){
   return("mean")
  }

 ),

 # private ----
 private = list(

  check_split_rule_internal = function(){

   check_arg_is_valid(arg_value = self$split_rule,
                      arg_name = 'split_rule',
                      valid_options = c("variance"))

  },

  check_pred_type_internal = function(oobag,
                                      pred_type = NULL,
                                      context = NULL){

   input <- pred_type %||% self$pred_type

   arg_name <- if(oobag) 'oobag_pred_type' else 'pred_type'

   if(is.null(context)){
    valid_options <- c("none", "mean", "leaf")
   } else {
    valid_options <- switch(
     context,
     'partial dependence' = c("mean"),
     'prediction' = c("mean", "leaf")
    )
    context <- paste(context, 'with regression forests')
   }

   check_arg_is_valid(arg_value = input,
                      arg_name = arg_name,
                      valid_options = valid_options,
                      context = context)

  },

  check_pred_horizon = function(pred_horizon = NULL,
                                boundary_checks = TRUE,
                                pred_type = NULL){

   # nothing to check
   NULL

  },

  check_oobag_eval_function_internal = function(oobag_fun){


   test_y <- seq(0, 1, length.out = 100)

   .y_mat <- matrix(test_y, ncol = 1)
   .w_vec <- rep(1, times = 100)
   .s_vec <- seq(0.9, 0.1, length.out = 100)

   test_output <- try(oobag_fun(y_mat = .y_mat,
                                w_vec = .w_vec,
                                s_vec = .s_vec),
                      silent = FALSE)

   if(is_error(test_output)){

    stop("oobag_fun encountered an error when it was tested. ",
         "Please make sure your oobag_fun works for this case:\n\n",
         "test_y <- seq(0, 1, length.out = 100)\n",
         "y_mat <- matrix(test_y, ncol = 1)\n",
         "w_vec <- rep(1, times = 100)\n",
         "s_vec <- seq(0.9, 0.1, length.out = 100)\n\n",
         "test_output <- oobag_fun(y_mat = y_mat, w_vec = w_vec, s_vec = s_vec)\n\n",
         "test_output should be a numeric value of length 1",
         call. = FALSE)

   }

   test_output

  },

  check_lincomb_R_function_internal = function(lincomb_R_function = NULL){

   input <- lincomb_R_function %||% self$lincomb_R_function

   .x_node <- matrix(rnorm(300), ncol = 3)
   .y_node <- matrix(rnorm(100), ncol = 1)
   .w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)

   out <- try(input(.x_node, .y_node, .w_node), silent = FALSE)

   if(is_error(out)){

    stop("user-supplied function encountered an error when it was tested. ",
         "Please make sure the function works for this case:\n\n",
         ".x_node <- matrix(rnorm(300), ncol = 3)\n",
         ".y_node <- matrix(rnorm(100), ncol = 1)\n",
         ".w_node <- matrix(rep(c(1,2,3,4), each = 25), ncol = 1)\n",
         "test_output <- your_function(.x_node, .y_node, .w_node)\n\n",
         "test_output should be a numeric matrix with 1 column and",
         " with nrow(test_output) = ncol(.x_node)",
         call. = FALSE)

   }

   out

  },

  init_control = function(){

   self$control <- orsf_control_regression(method = 'glm',
                                           scale_x = FALSE,
                                           max_iter = 1)

  },

  init_pred_type = function(){
   self$pred_type <- 'mean'
  },

  init_split_rule = function(){
   self$split_rule <- 'variance'
  },

  init_split_min_stat = function(){

   if(is.null(self$split_rule))
    stop("cannot init split_min_stat without split_rule", call. = FALSE)

   self$split_min_stat <- switch(self$split_rule, 'variance' = 0)

  },

  init_internal = function(){

   self$tree_type <- "regression"

   if(is.factor(self$data[[private$data_names$y]])){
    stop("Cannot fit regression trees to outcome ",
         private$data_names$y, " because it is a factor.",
         " Did you mean to use orsf_control_classification()?",
         call. = FALSE)
   }

   if(!is.function(self$control$lincomb_R_function) &&
      self$control$lincomb_type == 'net'){
    self$control$lincomb_R_function <- penalized_linreg
   }

   if(!self$oobag_pred_mode) self$oobag_eval_type <- "none"

   # use default if eval type was not specified by user
   if(self$oobag_pred_mode && is.null(self$oobag_eval_type)){
    self$oobag_eval_type <- "RSQ"
   }

  },

  init_pred = function(pred_horizon = NULL, pred_type = NULL,
                       boundary_checks = TRUE){

   if(!is.null(pred_horizon)){
    warning("pred_horizon does not impact predictions",
            " for regression forests", call. = FALSE)
   }

   if(!is.null(pred_type)){
    private$check_pred_type(oobag = FALSE, pred_type = pred_type)
   } else {
    pred_type <- self$pred_type %||% "mean"
   }

   self$pred_type <- pred_type

  },

  prep_y_internal = function(placeholder = FALSE){

   if(placeholder){
    private$y <- matrix(0, ncol = 1, nrow = 1)
    return()
   }

   # y is always 1 column for regression (for now)
   private$y <- as_matrix(private$y[[1]])
   colnames(private$y) <- private$data_names$y

  },

  predict_internal = function(simplify){

   # resize y to have the right number of columns
   private$y <- matrix(0, ncol = 1)

   cpp_args = private$prep_cpp_args(x = private$x,
                                    y = private$y,
                                    w = private$w,
                                    importance_type = 'none',
                                    pred_type = self$pred_type,
                                    pred_aggregate = self$pred_aggregate,
                                    oobag_pred = FALSE,
                                    pred_mode = TRUE,
                                    write_forest = FALSE,
                                    run_forest = TRUE)


   out <- do.call(orsf_cpp, args = cpp_args)$pred_new

   if(simplify) dim(out) <- NULL

   out

  }

 )
)
