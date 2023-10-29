//  Forest.cpp

#include <RcppArmadillo.h>
#include "ForestClassification.h"
#include "TreeClassification.h"

#include <memory>

using namespace arma;
using namespace Rcpp;

namespace aorsf {

ForestClassification::ForestClassification() { }

ForestClassification::ForestClassification(arma::uword n_class){

 this->n_class = n_class;

}

void ForestClassification::resize_pred_mat_internal(arma::mat& p){

 p.zeros(data->n_rows, data->n_cols_y);

}

void ForestClassification::compute_prediction_accuracy_internal(
  arma::mat& y,
  arma::vec& w,
  arma::mat& predictions,
  arma::uword row_fill
) {

 double cstat_sum = 0;

 for(uword i = 0; i < predictions.n_cols; i++){
  vec y_i = y.unsafe_col(i);
  vec p_i = predictions.unsafe_col(i);
  cstat_sum += compute_cstat_clsf(y_i, w, p_i);
 }

 oobag_eval(row_fill, 0) = cstat_sum / predictions.n_cols;

}

// growInternal() in ranger
void ForestClassification::plant() {

 trees.reserve(n_tree);

 // for (arma::uword i = 0; i < n_tree; ++i) {
 //  trees.push_back(std::make_unique<TreeClassification>());
 // }

}

}


