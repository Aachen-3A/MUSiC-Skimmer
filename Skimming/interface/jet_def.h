#ifndef jet_def_H
#define jet_def_H

#include <string>
#include "FWCore/ParameterSet/interface/InputTag.h"

//holds everything used to analyze a jet collection
struct jet_def{
   std::string   name;
   edm::InputTag MCLabel;
   edm::InputTag RecoLabel;
};

#endif
