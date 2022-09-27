/*-----------------------------------------------------------------------------
 This file is part of aorsf, which is distributed under the MIT license

 You should have received a copy of the MIT License along with aorsf.
 If not, see <http://www.gnu.org/licenses/>.

 Authors:
 - Byron C. Jaeger (http://byronjaeger.com)
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#include "globals.h"
#include "Data.h"
#include "Tree.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace aorsf;

// @description sample weights to mimic a bootstrap sample
// arma::vec bootstrap_sample_testthat(arma::mat& x,
//                                     arma::mat& y,
//                                     arma::vec& w) {
//
//  Data data = Data(x, y, w);
//  Data* data_ptr = &data;
//  return bootstrap_sample(data_ptr);
//
// }


// [[Rcpp::export]]
void orsf_cpp(arma::mat& x,
              arma::mat& y,
              arma::vec& w,
              int vi = 0,
              int sr = 1,
              int pt = 1){


 int mtry = 2;
 int leaf_min_obs = DEFAULT_LEAF_MIN_OBS_SURVIVAL;
 // int split_min_obs = DEFAULT_SPLIT_MIN_OBS;
 // int split_min_stat = DEFAULT_SPLIT_MIN_STAT;
 // int max_retry = DEFAULT_MAX_RETRY;
 // int n_split = DEFAULT_N_SPLIT;
 // int oobag_eval_every = 0;
 // int seed = 0;

 // VariableImportance variable_importance = static_cast<VariableImportance>(vi);
 // SplitRule split_rule = static_cast<SplitRule>(sr);
 // PredType pred_type = static_cast<PredType>(pt);

 // if( variable_importance == VI_NONE )
 //  Rcout << variable_importance << std::endl;

 Data data = Data(x, y, w);

 if(VERBOSITY > 0){
  Rcout << "------------ dimensions ------------"   << std::endl;
  Rcout << "N obs total: "     << data.get_n_rows() << std::endl;
  Rcout << "N columns total: " << data.get_n_cols() << std::endl;
  Rcout << "------------------------------------";
  Rcout << std::endl << std::endl;
 }

 Data* data_ptr = &data;

 Tree tree(data_ptr, leaf_min_obs, mtry);

 tree.draw_bootstrap_sample();

 tree.grow();



}
