#include<RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector which(LogicalVector x) {
  IntegerVector v = seq(0, x.length()-1);
  return v[x];
}

// [[Rcpp::export]]
List ostree_fit_cpp(NumericMatrix x_,
                    NumericMatrix y_,
                    Function cox_fit,
                    int mtry = 4,
                    int n_vars_lc = 2,
                    int n_cps = 5,
                    int grow_min_event = 10,
                    int grow_min_obs = 20){

  int n_obs = x_.nrow();
  int n_col = x_.ncol();
  double prob_sampled = 1.0 / n_obs;

  List tree_nodes;
  List tree_leaves;

  arma::mat x(x_.begin(), n_obs, n_col, false);
  arma::mat y(y_.begin(), n_obs, n_col, false);

  // s is the number of times you might get selected into
  // a bootstrap sample. Realistically this won't be >10.
  IntegerVector s = seq(0,10);

  // compute probability of being selected into the bootstrap
  // 0 times, 1, times, ..., 9 times, or 10 times.
  NumericVector probs = dbinom(s, n_obs, prob_sampled, false);

  // do the actual bootstrap sampling (with replacement)
  IntegerVector weights = sample(s, n_obs, true, probs);


  IntegerVector parts(n_obs);
  IntegerVector to_grow {0};

  // keep_rows: a logical vector that says which rows
  // got selected into the bootstrap sample.
  LogicalVector keep_rows = weights > 0;
  // grow_rows: an integer vector that is equivalent to
  // keep_rows, unless we are growing a node beyond the root.
  IntegerVector grow_rows = which(keep_rows);

  List rfit;

  for(IntegerVector::iterator i = to_grow.begin(); i != to_grow.end(); ++i){

    if(to_grow.length() > 1){
      IntegerVector grow_rows = which( (keep_rows) & (parts == *i) );
    }

    // check for constant columns using just the grow rows

    LogicalVector cols_lgl(n_col);

    for(int j = 0; j < n_col; j++){

      // reference to the jth column of x
      NumericMatrix::Column x_j = x_(_, j);

      for(int k=1; k < grow_rows.length(); k++){

        if(x_j[grow_rows[k]] != x_j[grow_rows[0]]){

          cols_lgl[j] = true;
          break;

        }

      }

    }

    IntegerVector cols_to_sample = which(cols_lgl);

    int mtry_temp = min(
      IntegerVector::create(mtry,
                            grow_rows.length(),
                            cols_to_sample.length())
    );

    IntegerVector cols_sampled = sample(cols_to_sample, mtry_temp);

    arma::mat x_sub = x.submat(as<arma::uvec>(grow_rows),
                               as<arma::uvec>(cols_sampled));

    arma::mat y_sub = y.rows(as<arma::uvec>(grow_rows));

    arma::vec wt_sub = as<arma::vec>(weights[grow_rows]);

    List rfit = cox_fit(x_sub,
                        y_sub,
                        wt_sub,
                        mtry_temp,
                        n_vars_lc);

    return(rfit);


    //Rcout << y_sub << std::endl;




  }

  return(
    List::create(
      Named("tree_nodes") = tree_nodes,
      _["tree_leaves"] = tree_leaves,
      _["rfit"] = rfit
    )
  );

}
