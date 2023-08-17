/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#ifndef TREE_H_
#define TREE_H_


#include "Data.h"
#include "globals.h"
#include "utility.h"

 namespace aorsf {

 class Tree {

 public:

  Tree();

  // deleting the copy constructor
  Tree(const Tree&) = delete;
  // deleting the copy assignment operator
  Tree& operator=(const Tree&) = delete;

  void init(Data* data,
            int seed,
            arma::uword mtry,
            double leaf_min_events,
            double leaf_min_obs,
            VariableImportance vi_type,
            SplitRule split_rule,
            double split_min_events,
            double split_min_obs,
            double split_min_stat,
            arma::uword split_max_cuts,
            arma::uword split_max_retry,
            LinearCombo lincomb_type,
            double lincomb_eps,
            arma::uword lincomb_iter_max,
            bool   lincomb_scale,
            double lincomb_alpha,
            arma::uword lincomb_df_target,
            arma::uword lincomb_ties_method);



  void sample_rows();

  void sample_cols();

  bool is_col_splittable(arma::uword j);

  bool is_node_splittable(arma::uword node_id);

  arma::uvec find_cutpoints();

  arma::uword split_node(arma::uvec& cuts_all);

  double score_logrank();

  void grow(arma::vec& vi_numer, arma::uvec& vi_denom);

  std::vector<arma::uvec>& get_coef_indices() {
   return(coef_indices);
  }

  // Pointer to original data
  Data* data;

  arma::uword n_cols_total;
  arma::uword n_rows_total;

  arma::uword n_rows_inbag;

  int seed;

  // views of data
  arma::mat x_inbag;
  arma::mat x_oobag;
  arma::mat x_node;

  arma::mat y_inbag;
  arma::mat y_oobag;
  arma::mat y_node;

  // the 'w' is short for 'weights'
  arma::vec w_inbag;
  arma::vec w_oobag;
  arma::vec w_node;

  // g_node indicates where observations will go when this node splits
  // 0 means go down to left node, 1 means go down to right node
  // the 'g' is short for 'groups'
  arma::uvec g_node;

  // which rows of data are held out while growing the tree
  arma::uvec rows_inbag;
  arma::uvec rows_oobag;
  arma::uvec rows_node;
  arma::uvec cols_node;

  // variable importance
  VariableImportance vi_type;

  // Random number generator
  std::mt19937_64 random_number_generator;

  // tree growing members
  arma::uword mtry;
  double leaf_min_events;
  double leaf_min_obs;

  // node split members
  SplitRule split_rule;
  double split_min_events;
  double split_min_obs;
  double split_min_stat;
  arma::uword split_max_cuts;
  arma::uword split_max_retry;

  // linear combination members
  LinearCombo lincomb_type;
  arma::vec   lincomb;
  arma::uvec  lincomb_sort;
  double      lincomb_eps;
  arma::uword lincomb_iter_max;
  bool        lincomb_scale;
  double      lincomb_alpha;
  arma::uword lincomb_df_target;
  arma::uword lincomb_ties_method;

  // prediction members
  double pred_horizon;


  // which node each inbag observation is currently in.
  arma::uvec node_assignments;

  // coefficients for linear combinations;
  // one row per variable (mtry rows), one column per node
  // leaf nodes have all coefficients=0
  std::vector<arma::vec> coef_values;

  // indices of the predictors used by a node
  std::vector<arma::uvec> coef_indices;

  // cutpoints used to split the nodes
  std::vector<double> cutpoint;

  // directions to the next node (right node = left node + 1)
  std::vector<arma::uword> child_left;

  // leaf values (only in leaf nodes)
  std::vector<arma::mat> leaf_values;

  // predicted values for out-of-bag rows
  arma::vec pred_oobag;

  // contains the node ID of each leaf node
  // e.g., if the first leaf is node 5, then the first value
  // of leaf_index is 5. When new data end in node 5, find which
  // value of leaf_index is 5, and go to that leaf to find the
  // predictions.
  std::vector<size_t> leaf_index;

 protected:



 };

 } // namespace aorsf

#endif /* TREE_H_ */
