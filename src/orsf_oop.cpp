/*-----------------------------------------------------------------------------
 This file is part of aorsf, which is distributed under the MIT license

 You should have received a copy of the MIT License along with aorsf.
 If not, see <http://www.gnu.org/licenses/>.

 Authors:
 - Byron C. Jaeger (http://byronjaeger.com)
#----------------------------------------------------------------------------*/

#include <RcppArmadillo.h>

#include "Data.h"
#include "globals.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace aorsf;



// [[Rcpp::export]]
void orsf_cpp(arma::mat& x,
              arma::vec& y_ctns,
              arma::ivec& y_intg,
              arma::vec& weights){

 Data data = Data(x, y_ctns, y_intg, weights);

 Rcout << "------------ dimensions ------------"   << std::endl;
 Rcout << "N obs total: "     << data.get_n_rows() << std::endl;
 Rcout << "N columns total: " << data.get_n_cols() << std::endl;
 Rcout << "------------------------------------";
 Rcout << std::endl << std::endl << std::endl;




}
