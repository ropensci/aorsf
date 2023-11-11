
# TODO:
# - write virtual function to get prediction type when calculating importance
# - defaults for split rule split min stat and oobag pred type in init_internal
# ObliqueForest class ----

ObliqueForest <- R6::R6Class(
 "ObliqueForest",
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
                        pred_type,
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

   # allow re-training (b/c why not).
   if(self$trained){ self$forest <- list() }

   cpp_args <- private$prep_cpp_args(...)

   cpp_output <- do.call(orsf_cpp, args = cpp_args)

   self$forest <- cpp_output$forest
   self$importance <- cpp_output$importance
   self$pred_oobag <- cpp_output$pred_oobag
   self$eval_oobag <- cpp_output$eval_oobag

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

  predict = function(new_data,
                     pred_horizon,
                     pred_type,
                     na_action,
                     boundary_checks,
                     n_thread,
                     verbose_progress,
                     pred_aggregate){

   public_state <- list(data             = self$data,
                        pred_horizon     = self$pred_horizon,
                        pred_type        = self$pred_type,
                        na_action        = self$na_action,
                        n_thread         = self$n_thread,
                        verbose_progress = self$verbose_progress,
                        pred_aggregate   = self$pred_aggregate)

   private_state <- list(data_rows_complete = private$data_rows_complete)

   # run checks before you assign new values to object.
   # otherwise, if a check throws an error, the object will
   # not be restored to its normal state.

   private$check_data(new = TRUE, data = new_data)
   private$check_na_action(new = TRUE, na_action = na_action)
   private$check_var_missing(new = TRUE, data = new_data, na_action)
   private$check_units(data = new_data)
   private$check_pred_type(oobag = FALSE, pred_type = pred_type)
   private$check_n_thread(n_thread)
   private$check_verbose_progress(verbose_progress)
   private$check_pred_aggregate(pred_aggregate)

   if(self$tree_type == 'survival')
    private$check_pred_horizon(pred_horizon, boundary_checks)

   self$data             <- new_data
   self$pred_horizon     <- pred_horizon
   self$pred_type        <- pred_type
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

     private$predict_internal()

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

  compute_vi = function(type_vi,
                        oobag_fun,
                        n_thread,
                        verbose_progress){

   # TODO: add checks here


   if(!is.null(oobag_fun)){
    check_oobag_fun(oobag_fun)
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

   cpp_args <-
    private$prep_cpp_args(oobag_eval_function = oobag_eval_function,
                          oobag_eval_type = oobag_eval_type,
                          importance_type = type_vi,
                          pred_type = switch(self$tree_type,
                                             'survival' = 'mort',
                                             'classification' = 'prob'),
                          pred_mode = FALSE,
                          pred_aggregate = TRUE,
                          # oobag should be FALSE for computing importance
                          # even though it is called 'oobag' importance.
                          oobag = FALSE,
                          write_forest = FALSE,
                          run_forest = TRUE,
                          n_thread = n_thread)

   out <- do.call(orsf_cpp, args = cpp_args)$importance
   rownames(out) <- colnames(private$x)
   out

  },

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
                                oobag,
                                type_output){

   public_state <- list(data         = self$data,
                        na_action    = self$na_action,
                        pred_horizon = self$pred_horizon)

   private_state <- list(data_rows_complete = private$data_rows_complete)

   private$check_boundary_checks(boundary_checks)
   private$check_pred_spec(pred_spec, boundary_checks)
   private$check_n_thread(n_thread)
   private$check_expand_grid(expand_grid)
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
   private$check_pred_type(pred_type, oobag = FALSE)

   pred_type <- pred_type %||% self$pred_type

   private$check_pred_horizon(pred_horizon, boundary_checks, pred_type)

   pred_horizon <- pred_horizon %||% self$pred_horizon %||% 1

   pred_horizon_order <- order(pred_horizon)
   pred_horizon_ordered <- pred_horizon[pred_horizon_order]

   # run checks before you assign new values to object.
   # otherwise, if a check throws an error, the object will
   # not be restored to its normal state.


   if(!oobag){
    private$check_data(new = TRUE, data = pd_data)
    # say new = FALSE to prevent na_action = 'pass'
    private$check_na_action(new = FALSE, na_action = na_action)
    private$check_var_missing(new = TRUE, data = pd_data, na_action)
    private$check_units(data = pd_data)
    self$data <- pd_data
   }

   self$pred_horizon <- pred_horizon
   self$na_action <- na_action

   # make a visible binding for CRAN
   id_variable = NULL

   private$init_data_rows_complete()
   private$prep_x()
   # y and w do not need to be prepped for prediction,
   # but they need to match orsf_cpp()'s expectations
   private$prep_y(placeholder = TRUE)
   private$w <- rep(1, nrow(private$x))


   if(oobag){ private$sort_inputs(sort_y = FALSE) }

   # the values in pred_spec need to be centered & scaled to match x_new,
   # which is also centered and scaled
   means <- private$data_means
   stdev <- private$data_stdev

   for(i in intersect(names(means), names(pred_spec))){
    pred_spec[[i]] <- (pred_spec[[i]] - means[i]) / stdev[i]
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

    x_cols <- list(match(names(pred_spec_new), colnames(private$x)))

    pred_spec_new <- list(as.matrix(pred_spec_new))

    pd_bind <- list(pred_spec)

   } else {

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

     x_cols[[i]] <- match(names(pred_spec_new[[i]]), colnames(private$x))
     pred_spec_new[[i]] <- as.matrix(pred_spec_new[[i]])

    }

   }



   cpp_args <- private$prep_cpp_args(x = private$x,
                                     y = private$y,
                                     w = private$w,
                                     importance_type = 'none',
                                     pred_type = pred_type,
                                     pred_mode = TRUE,
                                     pred_aggregate = TRUE,
                                     pred_horizon = pred_horizon_ordered,
                                     oobag = oobag,
                                     oobag_eval_type = 'none',
                                     n_thread = n_thread,
                                     write_forest = FALSE,
                                     run_forest = TRUE,
                                     verbosity = 0)


   pd_vals <- list()

   for(i in seq_along(pred_spec_new)){

    pd_vals_i <- list()

    x_pd <- private$x

    for(j in seq(nrow(pred_spec_new[[i]]))){

     x_pd[, x_cols[[i]]] <- pred_spec_new[[i]][j, ]

     cpp_args$x <- x_pd

     pd_vals_i[[j]] <- do.call(orsf_cpp, cpp_args)$pred_new

    }

    if(type_output == 'smry'){
     pd_vals_i <- lapply(
      pd_vals_i,
      function(x) {
       apply(x, 2, function(x_col){
        as.numeric(
         c(mean(x_col, na.rm = TRUE),
           quantile(x_col, probs = prob_values, na.rm = TRUE))
        )
       })
      }
     )
    }


    pd_vals[[i]] <- pd_vals_i

   }

   for(i in seq_along(pd_vals)){

    pd_bind[[i]]$id_variable <- seq(nrow(pd_bind[[i]]))

    for(j in seq_along(pd_vals[[i]])){


     pd_vals[[i]][[j]]

     if(self$tree_type == 'survival'){

      pd_vals[[i]][[j]] <- matrix(pd_vals[[i]][[j]],
                                  nrow=length(pred_horizon),
                                  byrow = T)

      rownames(pd_vals[[i]][[j]]) <- pred_horizon

     } else {

      pd_vals[[i]][[j]] <- t(pd_vals[[i]][[j]])

      if(self$tree_type == 'classification'){
       rownames(pd_vals[[i]][[j]]) <- self$class_levels
      }

     }

     if(type_output=='smry')
      colnames(pd_vals[[i]][[j]]) <- c('mean', prob_labels)
     else
      colnames(pd_vals[[i]][[j]]) <- c(paste(1:nrow(private$x)))

     # this will be null for non-survival objects
     ph <- rownames(pd_vals[[i]][[j]])

     pd_vals[[i]][[j]] <- as.data.frame(pd_vals[[i]][[j]])

     rownames(pd_vals[[i]][[j]]) <- NULL

     pd_vals[[i]][[j]][['pred_horizon']] <- ph

     if(type_output == 'ice'){

      pd_vals[[i]][[j]] <- melt_aorsf(
       data = pd_vals[[i]][[j]],
       id.vars = 'pred_horizon',
       variable.name = 'id_row',
       value.name = 'pred',
       measure.vars = setdiff(names(pd_vals[[i]][[j]]), 'pred_horizon'))

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

   }

   out <- rbindlist(pd_vals)

   # # missings may occur when oobag=TRUE and n_tree is small
   # if(type_output == 'ice') {
   #  out <- collapse::na_omit(out, cols = 'pred')
   # }

   ids <- c('id_variable')

   if(type_output == 'ice') ids <- c(ids, 'id_row')

   mid <- setdiff(names(out), c(ids, 'mean', prob_labels, 'pred'))

   end <- setdiff(names(out), c(ids, mid))

   setcolorder(out, neworder = c(ids, mid, end))

   if(self$tree_type == 'classification'){
    setnames(out, old = 'pred_horizon', new = 'class')
    out[, class := factor(class, levels = self$class_levels)]
    setkey(out, class)
   }

   if(self$tree_type == 'survival' && pred_type != 'mort')
    out[, pred_horizon := as.numeric(pred_horizon)]

   if(pred_type == 'mort'){
    out[, pred_horizon := NULL]
   }

   # not needed for summary
   if(type_output == 'smry')
    out[, id_variable := NULL]

   # put data back into original scale
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

   # silent print after modify in place
   out[]

   private$restore_state(public_state, private_state)

   # free up space
   private$x <- NULL
   private$y <- NULL
   private$w <- NULL

   out


  },

  select_variables = function(n_predictor_min, verbose_progress){

   public_state <- list(verbose_progress = self$verbose_progress,
                        forest           = self$forest,
                        control          = self$control)

   n_predictors <- length(private$data_names$x_original)

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

   cpp_args <- private$prep_cpp_args(pred_type = 'mort',
                                     oobag_pred = TRUE,
                                     importance_group_factors = TRUE,
                                     write_forest = FALSE)

   mtry_safe <- self$mtry


   while(n_predictors >= n_predictor_min){

    if(mtry_safe >= n_predictors){
     mtry_safe <- max(mtry_safe - 1, 1)
    }

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

   }

   collapse::na_omit(oob_data)

  },

  summarize_uni = function(n_variables = NULL,
                           pred_horizon = NULL,
                           pred_type = NULL,
                           importance_type = NULL){

   # check incoming values if they were specified.
   private$check_n_variables(n_variables)

   if(!is.null(pred_horizon)){
    private$check_pred_horizon(pred_horizon, boundary_checks = TRUE)
   }

   private$check_pred_type(pred_type, oobag = FALSE)
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
    'anova' = orsf_vi_anova(self, group_factors = TRUE),
    'negate' = orsf_vi_negate(self, group_factors = TRUE),
    'permute' = orsf_vi_permute(self, group_factors = TRUE),
    'none' = NULL
   )

   bounds <- private$data_bounds
   fctrs  <- private$data_fctrs
   n_obs  <- self$n_obs

   names_vi <- names(vi) %||% names_x

   pred_spec <- list_init(names_vi)[seq(n_variables)]

   for(i in names(pred_spec)){

    if(i %in% colnames(bounds)){

     pred_spec[[i]] <- unique(
      as.numeric(bounds[c('25%','50%','75%'), i])
     )

    } else if (i %in% fctrs$cols) {

     pred_spec[[i]] <- fctrs$lvls[[i]]

    }

   }

   pd_output <- orsf_pd_oob(object = self,
                            pred_spec = pred_spec,
                            expand_grid = FALSE,
                            pred_type = pred_type,
                            prob_values = c(0.25, 0.50, 0.75),
                            pred_horizon = pred_horizon)

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

  # getters

  get_names_x = function(ref_coded = FALSE){

   if(ref_coded) return(private$data_names$x_ref_code)

   return(private$data_names$x_original)

  },

  get_names_y = function(){
   return(private$data_names$y)
  },

  get_var_bounds = function(.name){

   if(.name %in% private$data_names$x_numeric)
    return(private$data_bounds[, .name])
   else
    return(private$data_fctrs$lvls[[.name]])

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

  x = NULL,
  y = NULL,
  w = NULL,

  importance_raw = NULL,

  mean_leaves = 0,

  # runs checks and sets defaults where needed
  init = function() {

   private$check_data()
   private$check_formula()

   if(is.null(self$control)){
    private$init_control()
   } else {
    private$check_control()
   }

   private$init_data()
   private$init_mtry()
   private$init_weights()

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
   private$check_pred_type(oobag = TRUE)
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
   stop("this method should only be called from derived classes")
  },
  init_data = function(){

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


   private$check_var_names(c(names_x_data, names_y_data))

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

   if(length(self$tree_seeds) == 1 && self$n_tree > 1){
     set.seed(self$tree_seeds)
     self$tree_seeds <- sample(self$n_tree*10, size = self$n_tree)
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

   if(is.null(self$pred_type)){
    self$oobag_pred_mode <- TRUE
   } else {
    self$oobag_pred_mode <- self$pred_type != "none"
   }

   if(!self$oobag_pred_mode) self$oobag_eval_type <- "none"

   if(self$oobag_pred_mode && self$sample_fraction == 1){
    stop(
     "cannot compute out-of-bag predictions if no samples are out-of-bag.",
     " Try setting sample_fraction < 1 or pred_type = 'none'.",
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

    private$w <- rep(1, self$n_obs)

   } else {

    private$check_weights()
    private$w <- self$weights

   }

  },
  init_oobag_eval_function = function(){

   if(is.null(self$oobag_eval_function)){

    self$oobag_eval_function <- function(x) x

   } else {

    private$check_oobag_eval_function()
    self$oobag_eval_type <- "User-specified function"

   }

  },

  # checkers
  check_data = function(data = NULL, new = FALSE){

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

    check_new_data_names(new_data  = input,
                         ref_names = private$data_names$x_original,
                         label_new = "new_data",
                         label_ref = 'training data')

    check_new_data_types(new_data  = input,
                         ref_names = private$data_names$x_original,
                         ref_types = private$data_types$x,
                         label_new = "new_data",
                         label_ref = 'training data')

    check_new_data_fctrs(new_data  = input,
                         names_x   = private$data_names$x_original,
                         fi_ref    = private$data_fctrs,
                         label_new = "new_data")

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

   for(i in private$data_names$x_original){

    if(collapse::allNA(input[[i]])){
     stop("column ", i, " has no observed values.",
          call. = FALSE)
    }

    if(any(is.infinite(input[[i]]))){
     stop("Please remove infinite values from ", i, ".",
          call. = FALSE)
    }

   }

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

   check_arg_length(
    arg_value = input,
    arg_name  = 'weights',
    expected_length = self$n_obs
   )

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
  check_pred_type = function(pred_type = NULL, oobag){

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

    private$check_pred_type_internal(oobag, pred_type)

   }



  },

  check_pred_type_internal = function(oobag, pred_type = NULL){

   stop("this method should be defined in a derived class.")

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

    test_time <- seq(from = 1, to = 5, length.out = 100)
    test_status <- rep(c(0,1), each = 50)

    .y_mat <- cbind(time = test_time, status = test_status)
    .w_vec <- rep(1, times = 100)
    .s_vec <- seq(0.9, 0.1, length.out = 100)

    test_output <- try(input(y_mat = .y_mat, w_vec = .w_vec, s_vec = .s_vec),
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
  check_verbose_progress = function(verbose_progress = NULL){

   input <- verbose_progress %||% self$verbose_progress

   check_arg_type(arg_value = input,
                  arg_name = 'verbose_progress',
                  expected_type = 'logical')

   check_arg_length(arg_value = input,
                    arg_name = 'verbose_progress',
                    expected_length = 1)

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

  check_oobag_pred_mode = function(oobag_pred_mode, label){

   check_arg_type(arg_value = oobag_pred_mode,
                  arg_name = label,
                  expected_type = 'logical')

   check_arg_length(arg_value = oobag_pred_mode,
                    arg_name = label,
                    expected_length = 1)


  },

  # computers

  compute_means = function(){

   numeric_data <- select_cols(self$data, private$data_names$x_numeric)
   private$data_means <- collapse::fmean(numeric_data, w = self$weights)

  },


  compute_modes = function(){

   private$data_modes <- vapply(
    select_cols(self$data, private$data_fctrs$cols),
    collapse::fmode,
    FUN.VALUE = integer(1),
    w = self$weights
   )

  },

  compute_stdev = function(){

   numeric_data <- select_cols(self$data, private$data_names$x_numeric)
   private$data_stdev <- collapse::fsd(numeric_data, w = self$weights)

  },

  compute_bounds = function(){

   numeric_data <- select_cols(self$data, private$data_names$x_numeric)

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

  prep_y = function(placeholder = FALSE){

   private$y <- select_cols(self$data, private$data_names$y)

   if(self$na_action == 'omit' && !placeholder)
    private$y <- private$y[private$data_rows_complete, ]

   private$prep_y_internal(placeholder)

  },

  prep_w = function(){

   # re-initialize
   private$init_weights()

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
                       "none" = 0,
                       "negate" = 1,
                       "permute" = 2,
                       "anova" = 3),
    vi_max_pvalue = .dots$importance_max_pvalue %||% self$importance_max_pvalue,
    leaf_min_events = .dots$leaf_min_events %||% self$leaf_min_events %||% 1,
    leaf_min_obs = .dots$leaf_min_obs %||% self$leaf_min_obs,
    split_rule_R = switch(self$split_rule,
                          "logrank" = 1,
                          "cstat" = 2,
                          "gini" = 3),
    split_min_events = .dots$split_min_events %||% self$split_min_events %||% 1,
    split_min_obs = .dots$split_min_obs %||% self$split_min_obs,
    split_min_stat = .dots$split_min_stat %||% self$split_min_stat,
    split_max_cuts = .dots$split_max_cuts %||% self$n_split,
    split_max_retry = .dots$split_max_retry %||% self$n_retry,
    lincomb_R_function = self$control$lincomb_R_function,
    lincomb_type_R = switch(self$control$lincomb_type,
                            'glm' = 1,
                            'random' = 2,
                            'net' = 3,
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
                         "none" = 0,
                         "risk" = 1,
                         "surv" = 2,
                         "chf"  = 3,
                         "mort" = 4,
                         "prob" = 6,
                         "class" = 7,
                         "leaf" = 8),
    pred_mode = .dots$pred_mode %||% FALSE,
    pred_aggregate = .dots$pred_aggregate %||% (self$pred_type != 'leaf'),
    pred_horizon = .dots$pred_horizon %||% self$pred_horizon %||% 1,
    oobag = .dots$oobag %||% self$oobag_pred_mode,
    oobag_R_function = .dots$oobag_eval_function %||% self$oobag_eval_function,
    oobag_eval_type_R = switch(
     tolower(.dots$oobag_eval_type %||% self$oobag_eval_type),
     "none" = 0,
     "harrell's c-index" = 1,
     "auc-roc" = 1,
     "user-specified function" = 2
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

  sort_inputs = function(){
   NULL
  },

  # cleaners

  clean_importance = function(){

   out <- self$importance

   # nan indicates a variable was never used
   out[is.nan(out)] <- 0

   rownames(out) <- private$data_names$x_ref_code

   private$importance_raw <- out

   if(self$importance_group_factors){

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
 cloneable = FALSE,
 public = list(

  get_max_time = function(){
   return(private$max_time)
  },

  leaf_min_events = NULL,
  split_min_events = NULL,
  pred_horizon = NULL

 ),

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
  check_pred_type_internal = function(oobag, pred_type = NULL){

   input <- pred_type %||% self$pred_type

   arg_name <- if(oobag) 'oobag_pred_type' else 'pred_type'

   check_arg_is_valid(arg_value = input,
                      arg_name = arg_name,
                      valid_options = c("none", "surv", "risk",
                                        "chf", "mort", "leaf"))

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

  init_internal = function(){

   self$tree_type <- "survival"

   if(!is.function(self$control$lincomb_R_function) &&
      self$control$lincomb_type == 'net'){
    self$control$lincomb_R_function <- penalized_cph
   }

   self$split_rule <- self$split_rule %||% 'logrank'
   self$pred_type <- self$pred_type %||% 'surv'
   self$split_min_stat <- self$split_min_stat %||%
    switch(self$split_rule, 'logrank' = 3.841459, 'cstat' = 0.50)

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

   # if pred_horizon is unspecified, provide sensible default
   # if it is specified, check for correctness
   if(is.null(self$pred_horizon)){
    self$pred_horizon <- collapse::fmedian(y[, 1])
   } else {
    private$check_pred_horizon(self$pred_horizon, boundary_checks = TRUE)
   }

   private$check_leaf_min_events()
   private$check_split_min_events()

   # use default if eval type was not specified by user
   if(self$oobag_pred_mode && is.null(self$oobag_eval_type)){
    self$oobag_eval_type <- "Harrell's C-index"
   }

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

   # mortality predictions should always be 1 column
   # b/c they do not depend on the prediction horizon
   if(self$pred_type == 'mort'){

    self$eval_oobag$stat_values <-
     self$eval_oobag$stat_values[, 1L, drop = FALSE]

    self$pred_oobag <- self$pred_oobag[, 1L, drop = FALSE]

   }

  },
  clean_pred_new_internal = function(preds){

    # output in the same order as user's pred_horizon vector
   preds <- preds[, order(private$pred_horizon_order), drop = FALSE]

   preds

  },

  predict_internal = function(){

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

    return(simplify2array(results))

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
 cloneable = FALSE,
 public = list(

  n_class = NULL,

  class_levels = NULL

 ),
 private = list(

  check_split_rule_internal = function(){

   check_arg_is_valid(arg_value = self$split_rule,
                      arg_name = 'split_rule',
                      valid_options = c("gini", "cstat"))

  },
  check_pred_type_internal = function(oobag, pred_type = NULL){

   input <- pred_type %||% self$pred_type

   arg_name <- if(oobag) 'oobag_pred_type' else 'pred_type'

   check_arg_is_valid(arg_value = input,
                      arg_name = arg_name,
                      valid_options = c("none", "prob", "class", "leaf"))

  },

  check_pred_horizon = function(pred_horizon = NULL,
                                boundary_checks = TRUE,
                                pred_type = NULL){

   # nothing to check
   NULL

  },

  init_control = function(){

   self$control <- orsf_control_classification(method = 'glm',
                                               scale_x = FALSE,
                                               max_iter = 1)

  },

  init_internal = function(){

   self$tree_type <- "classification"

   if(!is.function(self$control$lincomb_R_function) &&
      self$control$lincomb_type == 'net'){
    self$control$lincomb_R_function <- penalized_logreg
   }

   self$split_rule <- self$split_rule %||% 'gini'
   self$pred_type <- self$pred_type %||% 'prob'
   self$split_min_stat <- self$split_min_stat %||%
    switch(self$split_rule, 'gini' = 0, 'cstat' = 0.50)

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

  prep_y_internal = function(placeholder = FALSE){

   if(placeholder){
    private$y <- matrix(0, ncol = self$n_class-1, nrow = 1)
    return()
   }

   # y is always 1 column for classification (right?)
   y <- private$y[[1]]

   if(!is.factor(y)) y <- as.factor(y)

   n_class <- length(levels(y))

   y <- as.numeric(y) - 1

   private$y <- expand_y_clsf(as_matrix(y), n_class)

  },

  predict_internal = function(){

   # resize y to have the right number of columns
   private$y <- matrix(0, ncol = self$n_class-1)

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

   # no further cleaning needed
   do.call(orsf_cpp, args = cpp_args)$pred_new

  }

 )
)
