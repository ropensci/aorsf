

set.seed(1)

formula <- Surv(time, status) ~ . - id

aorsf <- orsf(
 formula = formula,
 data = pbc_orsf,
 n_tree = 10,
 leaf_min_obs = 20
)

formula_terms <- stats::terms(formula, data=pbc_orsf)
names_x_data <- attr(formula_terms, 'term.labels')

fi <- fctr_info(pbc_orsf, names_x_data)
x  <- as.matrix(ref_code(pbc_orsf, fi, names_x_data))

for( tr in seq(10) ){

 tree <- aorsf$forest[[tr]]
 leaf_assigned <- rep(0, nrow(x))

 for (j in seq(0, ncol(tree$betas)-1 ) ) {

  jj <- j+1

  obs_in_node <- which(leaf_assigned==j)

  if(tree$children_left[jj] != 0){

   lc <- x[obs_in_node, (tree$col_indices[, jj] + 1)] %*% tree$betas[, jj, drop = F]

   going_left <- lc <= tree$cut_points[jj]
   going_right <- !going_left

   leaf_assigned[obs_in_node[going_left]] <- tree$children_left[jj]
   leaf_assigned[obs_in_node[going_right]] <- tree$children_left[jj] + 1

  }

 }

 test_that(
  desc = 'check pred_leaf with R script',
  code = {
   expect_equal(as.numeric(ostree_pred_leaf_testthat(tree = tree, x_pred_ = x)),
                as.numeric(leaf_assigned))
  }
 )

}
