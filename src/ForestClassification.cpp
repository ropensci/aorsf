//  Forest.cpp

#include <RcppArmadillo.h>
#include "ForestClassification.h"
#include "TreeClassification.h"

#include <memory>

using namespace arma;
using namespace Rcpp;

namespace aorsf {

ForestClassification::ForestClassification() { }


void ForestClassification::resize_pred_mat_internal(arma::mat& p){

 p.zeros(data->n_rows, 1);

}

void ForestClassification::compute_prediction_accuracy_internal(
  arma::mat& y,
  arma::vec& w,
  arma::mat& predictions,
  arma::uword row_fill
) {

 // oobag_eval(row_fill, 0) = compute_cstat_clsf(y, w, predictions);

}

}


