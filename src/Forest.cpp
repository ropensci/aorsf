/*-----------------------------------------------------------------------------
 This file is part of aorsf.
 Author: Byron C Jaeger
 aorsf may be modified and distributed under the terms of the MIT license.
#----------------------------------------------------------------------------*/

#include <armadillo>

#include "globals.h"
#include "Forest.h"

 namespace aorsf {

 Forest::Forest() :
   data(0),
   mtry(0),
   max_retry(DEFAULT_MAX_RETRY),
   split_rule(DEFAULT_SPLITRULE_SURVIVAL),
   n_split(DEFAULT_N_SPLIT),
   leaf_min_obs(0),
   split_min_obs(0),
   split_min_stat(DEFAULT_SPLIT_MIN_STAT),
   pred_type(DEFAULT_PRED_TYPE),
   oobag_eval_every(0),
   variable_importance(DEFAULT_IMPORTANCE),
   seed(0) {

   }



 } // namespace aorsf
