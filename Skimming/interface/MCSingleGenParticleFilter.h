#ifndef MCSINGLEGENPARTICLEFILTER
#define MCSINGLEGENPARTICLEFILTER

//
// Filter based on the CMSSW "MCSingleParticleFilter".
//
// This module filters events based on Pythia particleID, the pt, eta. In
// contrast to MCSingleParticleFilter it takes the reco::GenParticle
// information rather than the HepMC::GenEvent information.
//
// Inherits from EDFilter
//

#include <memory>

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//

class MCSingleGenParticleFilter : public edm::EDFilter {
public:
   explicit MCSingleGenParticleFilter(const edm::ParameterSet&);
   ~MCSingleGenParticleFilter();


   virtual bool filter(edm::Event&, const edm::EventSetup&);
private:
   // ----------member data ---------------------------
    edm::InputTag label_;
    std::vector< int > particleID;
    std::vector< double > ptMin;
    std::vector< double > etaMin;
    std::vector< double > etaMax;
    std::vector< int > status;
};

#endif /*MCSINGLEGENPARTICLEFILTER*/
