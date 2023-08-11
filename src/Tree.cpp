/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include "Tree.h"

 using namespace arma;
 using namespace Rcpp;

 namespace aorsf {

 Tree::Tree(){ }

 void Tree::init(Data* data,
                 int seed,
                 arma::uword mtry,
                 double leaf_min_events,
                 double leaf_min_obs,
                 SplitRule split_rule,
                 double split_min_events,
                 double split_min_obs,
                 double split_min_stat,
                 arma::uword split_max_retry,
                 LinearCombo lincomb_type,
                 double lincomb_eps,
                 arma::uword lincomb_iter_max,
                 bool lincomb_scale,
                 double lincomb_alpha,
                 arma::uword lincomb_df_target){


  this->data = data;
  this->n_cols = data->n_cols;

  this->seed = seed;
  this->mtry = mtry;
  this->leaf_min_events = leaf_min_events;
  this->leaf_min_obs = leaf_min_obs;
  this->split_rule = split_rule;
  this->split_min_events = split_min_events;
  this->split_min_obs = split_min_obs;
  this->split_min_stat = split_min_stat;
  this->split_max_retry = split_max_retry;
  this->lincomb_type = lincomb_type;
  this->lincomb_eps = lincomb_eps;
  this->lincomb_iter_max = lincomb_iter_max;
  this->lincomb_scale = lincomb_scale;
  this->lincomb_alpha = lincomb_alpha;
  this->lincomb_df_target = lincomb_df_target;

 }

 void Tree::bootstrap(){

  // Initialize random number generator and set seed
  random_number_generator.seed(seed);

  uword i, draw, n = data->n_rows;

  // Start with all samples OOB
  uvec boot_wts(n, fill::zeros);

  std::uniform_int_distribution<uword> unif_dist(0, n - 1);

  // sample with replacement
  for (i = 0; i < n; ++i) {
   draw = unif_dist(random_number_generator);
   ++boot_wts[draw];
  }

  // multiply boot_wts by user specified weights.
  if(data->has_weights){
   boot_wts = boot_wts = boot_wts % data->w;
  }

  uvec rows_inbag = find(boot_wts > 0);
  this->x_inbag = data->x_rows(rows_inbag);
  this->y_inbag = data->y_rows(rows_inbag);
  this->w_inbag = data->w_subvec(rows_inbag);
  this->rows_oobag = find(boot_wts == 0);

  // all observations start in node 0
  this->rows_node = linspace<uvec>(0, x_inbag.n_rows-1, x_inbag.n_rows);

  rows_inbag.clear();
  boot_wts.clear();

 }

 bool Tree::is_col_valid(uword j){

  vec status = y_inbag.unsafe_col(1);

  uvec::iterator i;

  // initialize as 0 but do not make comparisons until x_first_value
  // is formally defined at the first instance of status == 1
  double x_first_value=0;

  bool x_first_undef = true;

  for (i = rows_node.begin(); i != rows_node.end(); ++i) {

   if(status[*i] == 1){

    if(x_first_undef){

     x_first_value = x_inbag.at(*i, j);
     x_first_undef = false;

    } else {

     if(x_inbag.at(*i, j) != x_first_value){
      return(true);
     }

    }

   }

  }

  return(false);

 }

 void Tree::sample_cols(){

  // Start empty
  std::vector<uword> cols_assessed, cols_accepted;
  cols_assessed.reserve(n_cols);
  cols_accepted.reserve(mtry);

  std::uniform_int_distribution<uword> unif_dist(0, n_cols - 1);

  uword i, draw;

  for (i = 0; i < n_cols; ++i) {

   draw = unif_dist(random_number_generator);

   // Ensure the drawn number is not already in the sample
   while (std::find(cols_assessed.begin(),
                    cols_assessed.end(),
                    draw) != cols_assessed.end()) {

    draw = unif_dist(random_number_generator);


   }

   cols_assessed.push_back(draw);

   if(is_col_valid(draw)){
    cols_accepted.push_back(draw);
   }

   if(cols_accepted.size() == mtry) break;

  }

  this->cols_node = uvec(cols_accepted.data(),
                         cols_accepted.size(),
                         false,
                         true);

 }

 void Tree::grow(){

  bootstrap();

  sample_cols();

  coef_indices.push_back(cols_node);

  node_assignments.zeros(x_inbag.n_rows);

  uvec nodes_to_grow(1, fill::zeros);
  uvec rows_node;

  if(VERBOSITY > 0){

   uword temp_uword_1, temp_uword_2;

   if(x_inbag.n_rows < 5)
    temp_uword_1 = x_inbag.n_rows-1;
   else
    temp_uword_1 = 5;

   if(x_inbag.n_cols < 5)
    temp_uword_2 = x_inbag.n_cols-1;
   else
    temp_uword_2 = 4;

   Rcout << "   ---- view of x_inbag ---- " << std::endl << std::endl;
   Rcout << round(x_inbag.submat(0, 0, temp_uword_1, temp_uword_2));
   Rcout << std::endl << std::endl;

  }



 } // Tree::grow

 } // namespace aorsf

