
//  ForestClassification.h

#ifndef FORESTCLASSIFICATION_H
#define FORESTCLASSIFICATION_H

#include "Data.h"
#include "globals.h"
#include "Forest.h"

namespace aorsf {

class ForestClassification: public Forest {

public:

 ForestClassification();

 virtual ~ForestClassification() override = default;

 ForestClassification(const ForestClassification&) = delete;
 ForestClassification& operator=(const ForestClassification&) = delete;

 void resize_pred_mat_internal(arma::mat& p) override;

 void compute_prediction_accuracy_internal(
   arma::mat& y,
   arma::vec& w,
   arma::mat& predictions,
   arma::uword row_fill
 ) override;

};

}



#endif /* ForestClassification_H */
