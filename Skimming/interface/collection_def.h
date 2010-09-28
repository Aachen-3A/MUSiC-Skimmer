#ifndef collection_def_H
#define collection_def_H

#include <string>
#include "FWCore/Utilities/interface/InputTag.h"

//holds everything used to analyze a collection
struct collection_def{
   std::string   name;
   edm::InputTag MCLabel;
   edm::InputTag RecoLabel;
};

#endif
