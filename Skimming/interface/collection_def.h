// -*- C++ -*-
// Copyright [2015] <RWTH Aachen, III. Phys. Inst. A>

#ifndef SKIMMING_INTERFACE_COLLECTION_DEF_H_
#define SKIMMING_INTERFACE_COLLECTION_DEF_H_

namespace pat {
class Jet;
}

template< class T > class Selector;

#include <string>
#include <vector>
#include <utility>

#include "FWCore/Utilities/interface/InputTag.h"

struct jet_def{
    std::string   name;
    edm::InputTag MCLabel;
    edm::InputTag RecoLabel;
    bool          isPF;
    std::vector< std::pair< std::string, Selector<pat::Jet>* > >   IDs;
};

#endif  // SKIMMING_INTERFACE_COLLECTION_DEF_H_
