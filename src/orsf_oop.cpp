/*-----------------------------------------------------------------------------
 This file is part of aorsf, which is distributed under the MIT license

 You should have received a copy of the MIT License along with aorsf.
 If not, see <http://www.gnu.org/licenses/>.

 Authors:
 - Byron C. Jaeger (http://byronjaeger.com)
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#include "Data.h"
#include "globals.h"
#include "Tree.h"
#include "Forest.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace aorsf;

// @description sample weights to mimic a bootstrap sample
arma::vec bootstrap_sample(const Data* data) {

 // s is the number of times you might get selected into
 // a bootstrap sample. Realistically this won't be >10,
 Rcpp::IntegerVector s = Rcpp::seq(0.0, 10);

 // compute probability of being selected into the bootstrap
 // 0 times, 1, times, ..., 9 times, or 10 times.

 arma::uword n_rows = data->get_n_rows();

 Rcpp::NumericVector probs = Rcpp::dbinom(s, n_rows, 1.0/n_rows, false);

 arma::vec boot_wts = Rcpp::as<arma::vec>(
  Rcpp::RcppArmadillo::sample(s, n_rows, true, probs)
 );

 if(data->has_wts()){

  boot_wts = boot_wts % data->get_weights();

 }

 return(boot_wts);

}


// @description sample weights to mimic a bootstrap sample
// [[Rcpp::export]]
arma::vec bootstrap_sample_testthat(arma::mat& x,
                                    arma::vec& y_ctns,
                                    arma::ivec& y_intg,
                                    arma::vec& weights) {

 Data data = Data(x, y_ctns, y_intg, weights);
 Data* data_ptr = &data;
 return bootstrap_sample(data_ptr);

}


// [[Rcpp::export]]
void orsf_cpp(arma::mat& x,
              arma::vec& y_ctns,
              arma::ivec& y_intg,
              arma::vec& weights,
              const int vi = 0,
              const int sr = 1,
              const int pt = 1){


 const int mtry = 2;
 const int max_retry = DEFAULT_MAX_RETRY;
 const int n_split = DEFAULT_N_SPLIT;
 const int leaf_min_obs = DEFAULT_LEAF_MIN_OBS_SURVIVAL;
 const int split_min_obs = DEFAULT_SPLIT_MIN_OBS;
 const int split_min_stat = DEFAULT_SPLIT_MIN_STAT;
 const int oobag_eval_every = 0;
 const int seed = 0;

 VariableImportance variable_importance = static_cast<VariableImportance>(vi);
 SplitRule split_rule = static_cast<SplitRule>(sr);
 PredType pred_type = static_cast<PredType>(pt);

 // if( variable_importance == VI_NONE ) Rcout << variable_importance << std::endl;

 Data data = Data(x, y_ctns, y_intg, weights);

 Data* data_ptr = &data;

<<<<<<< HEAD
 Rcout << "------------ dimensions ------------" << std::endl;
 Rcout << "N obs total: "     << data.n_rows     << std::endl;
 Rcout << "N columns total: " << data.n_cols     << std::endl;
 Rcout << "mtry: "            << mtry            << std::endl;
 Rcout << "------------------------------------";
 Rcout << std::endl << std::endl << std::endl;

 Tree tree(data_ptr,
           mtry,
           max_retry,
           split_rule,
           n_split,
           leaf_min_obs,
           split_min_obs,
           split_min_stat,
           pred_type,
           oobag_eval_every,
           variable_importance,
           seed);

 Rcout << tree.coef << std::endl;

 tree.guess_max_nodes();

 Rcout << tree.coef << std::endl;


}
