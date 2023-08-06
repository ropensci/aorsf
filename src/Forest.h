
//  Forest.h

#ifndef Forest_H
#define Forest_H

#include "Data.h"
#include "globals.h"

namespace aorsf {

class Forest {

public:

 // Constructor

 Forest();

 // deleting the copy constructor
 Forest(const Forest&) = delete;
 // deleting the copy assignment operator
 Forest& operator=(const Forest&) = delete;

 // Methods

 void init(std::unique_ptr<Data> input_data,
           int n_tree,
           Rcpp::IntegerVector& tree_seeds,
           Rcpp::List& tree_params);

 // virtual void initInternal() = 0;

 // Grow or predict
 void run(bool verbose, bool compute_oob_error);

 Rcpp::IntegerVector get_bootstrap_select_times(){
  return bootstrap_select_times;
 }

 Rcpp::NumericVector get_bootstrap_select_probs(){
  return bootstrap_select_probs;
 }


 // Member variables

 Rcpp::IntegerVector bootstrap_select_times;
 Rcpp::NumericVector bootstrap_select_probs;

 int n_tree;

 Rcpp::IntegerVector tree_seeds;
 Rcpp::List tree_objects;

 arma::uword n_split;
 arma::uword mtry;

 double leaf_min_events;
 double leaf_min_obs;

 std::unique_ptr<Data> data;

};

}



#endif /* Forest_H */
