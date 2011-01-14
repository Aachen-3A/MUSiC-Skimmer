#ifndef collection_def_H
#define collection_def_H

#include <string>
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"

//holds everything used to analyze a collection
struct collection_def{
   std::string   name;
   edm::InputTag MCLabel;
   edm::InputTag RecoLabel;
};

// typedef std::vector< std::pair< std::string, Selector<pat::Jet>* > > jet_id_list;
struct jet_def{
   std::string   name;
   edm::InputTag MCLabel;
   edm::InputTag RecoLabel;
   bool          isPF;
   std::vector< std::pair< std::string, Selector<pat::Jet>* > >   IDs;
};
#endif
