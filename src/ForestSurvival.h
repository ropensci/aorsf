
//  Forest.h

#ifndef ForestSurvival_H
#define ForestSurvival_H

#include "Data.h"
#include "globals.h"
#include "Forest.h"

namespace aorsf {

class ForestSurvival: public Forest {

public:

 ForestSurvival();

 ForestSurvival(const ForestSurvival&) = delete;
 ForestSurvival& operator=(const ForestSurvival&) = delete;


};

}



#endif /* Forest_H */
