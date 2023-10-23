
//  Forest.h

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

};

}



#endif /* Forest_H */
