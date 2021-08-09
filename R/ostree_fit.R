



ostree_fit <- function(x,
                       y,
                       mtry = 4,
                       n_vars_lc = 2,
                       n_cps = 5,
                       leaf_min_event = 10,
                       leaf_min_obs = 20){

 n_event <- sum(y[,'status'])
 n_obs <- length(y)

 nodes <- tree_leaves <- list()

 parent <- NA_integer_

 weights <- sample(x = 0:10,
                   size = nrow(x),
                   replace = TRUE,
                   prob = dbinom(0:10, size = nrow(x), prob = 1/nrow(x)))

 keep_rows <- weights > 0

 parts <- rep(0, nrow(x))

 to_grow <- n_event > leaf_min_event & n_obs > leaf_min_obs

 if(to_grow) parts_to_grow <- 0 else return(NULL)

 repeat{

  for(p in parts_to_grow){

   grow_rows <- keep_rows & parts == p

   if(p > 0){
    left <- unlist(lapply(nodes, function(x) x$child_left))
    right <- unlist(lapply(nodes, function(x) x$child_right))

    parent <- if(p %in% left)
     as.integer(names(left)[left == p])
    else
     as.integer(names(right)[right == p])
   }

   cnst_cols <- apply(
    x[grow_rows, , drop = FALSE],
    MARGIN = 2,
    FUN = is_constant
   )

   if(any(cnst_cols)) {
    cnames <- colnames(x)[-which(cnst_cols)]
   } else
    cnames <- colnames(x)

   grow_cols <- sample(cnames,
                       size = min(mtry, length(cnames)),
                       replace = FALSE)


   fit <- try(coxph.fit(x = x[grow_rows, grow_cols, drop = FALSE],
                        y = y[grow_rows],
                        strata = NULL,
                        weights = weights[grow_rows],
                        control = coxph.control(toler.inf = Inf),
                        method = 'breslow',
                        rownames = NULL,
                        resid = FALSE,
                        nocenter = c(-1,0,1)))

   if(is_error(fit)) browser()

   coef_na <- is.na(fit$coefficients) | fit$coefficients == 0

   beta <- fit$coefficients[!coef_na]
   se <- sqrt(diag(fit$var)[!coef_na])
   pv <- pchisq((beta/se)^2, df = 1, lower.tail = FALSE)

   keep_cols <- names(sort(pv))[seq(min(length(pv),n_vars_lc))]

   fit <- try(coxph.fit(x = x[grow_rows, keep_cols, drop = FALSE],
                        y = y[grow_rows],
                        strata = NULL,
                        weights = weights[grow_rows],
                        control = coxph.control(toler.inf = Inf),
                        init = coef(fit)[keep_cols],
                        method = 'breslow',
                        rownames = NULL,
                        resid = FALSE,
                        nocenter = c(-1,0,1)))

   if(is_error(fit)) browser()

   lc <- fit$linear.predictors
   lc_uni <- unique(lc)
   lc_range <- range(lc)
   lc_uni_inner <- setdiff(lc_uni, lc_range)
   lc_uni_inner_length <- length(lc_uni_inner)


   if(lc_uni_inner_length == 0){
    cps <- mean(lc_range)
   } else if(lc_uni_inner_length == 1){
    cps <- lc_uni_inner
   } else if (lc_uni_inner_length < n_cps){
    cps <- lc_uni_inner
   } else {
    cps <- setdiff(
     unique(quantile(lc, probs = seq(0.1, 0.9, length.out = n_cps))),
     lc_range
    )
   }

   if(length(cps) > 1){
    lrt <- vapply(
     X = cps,
     FUN = function(cp)
      lrtestC(y[grow_rows, 'time'],
              y[grow_rows, 'status'],
              grp = as.numeric(lc < cp)),
     FUN.VALUE = double(length = 1)
    )
    cp <- cps[which.max(lrt)]
   } else {
    cp <- cps
   }

   # nn is short for new node
   nn_last <- max(parts)
   nn_left <- nn_last + 1L
   nn_right <- nn_last + 2L

   # this sum is subtracted from XB in coxfit's linear.predictors
   # -> add to cut-point that will be used on future linear predictors
   cp <- cp + sum(coef(fit) * fit$means)

   parts[grow_rows] <- fifelse(
    test = x[grow_rows, keep_cols, drop=FALSE] %*% coef(fit) < cp,
    yes = nn_left, no = nn_right
   )


   nodes[[as.character(p)]] <-
    list(beta = coef(fit)[keep_cols],
         parent = parent,
         child_left = nn_left,
         child_right = nn_right,
         cutpoint = cp)

   n_leaf_event <- tapply(y[grow_rows,'status'], parts[grow_rows], sum)

   n_leaf_obs <- table(parts[grow_rows])

   new_part_is_leaf <-
    n_leaf_event < leaf_min_event | n_leaf_obs < leaf_min_obs

   if(any(new_part_is_leaf)){
    for(l in names(new_part_is_leaf)[new_part_is_leaf]){
     nodes[[as.character(l)]] <- list(parent = p)
    }
   }

  }

  n_leaf_event <- tapply(y[keep_rows,'status'], parts[keep_rows], sum)

  n_leaf_obs <- table(parts[keep_rows])

  to_grow <- n_leaf_event>=leaf_min_event & n_leaf_obs>=leaf_min_obs

  parts_to_grow <- as.integer(names(n_leaf_event)[to_grow])

  leaves <- setdiff(as.integer(names(n_leaf_event)), parts_to_grow)

  if(!is_empty(leaves))
   for(l in leaves){
    indx <- parts == l
    tree_leaves[[as.character(l)]] <-
     unclass(survival::survfit(
      Surv(y[indx, 'time'], y[indx, 'status'])~1,
      type = "kaplan-meier",
      se.fit=FALSE
     ))[c('time', 'surv')]
   }

  if(is_empty(parts_to_grow)) break

 }

 for(i in seq_along(tree_leaves)){

  tree_leaves[[i]]$path <- try(
   rev(get_parents(nodes, names(tree_leaves)[i]))
  )

  if(is_error(tree_leaves[[i]]$path)) browser()

 }

 tree_leaves

}
