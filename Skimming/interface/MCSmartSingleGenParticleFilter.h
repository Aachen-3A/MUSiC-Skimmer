#ifndef MCSMARTSINGLEGENPARTICLEFILTER
#define MCSMARTSINGLEGENPARTICLEFILTER

//
// Filter based on the CMSSW "MCSmartSingleParticleFilter".
//
// This module filters events based on Pythia particleID, the pt, eta,
// production vertex. In contrast to MCSmartSingleParticleFilter it takes the
// reco::GenParticle information rather than the HepMC::GenEvent information.
//
// Inherits from EDFilter
//

#include <memory>

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class declaration
//

class MCSmartSingleGenParticleFilter : public edm::EDFilter {
public:
   explicit MCSmartSingleGenParticleFilter( const edm::ParameterSet& );
   ~MCSmartSingleGenParticleFilter();

   virtual bool filter( edm::Event&, const edm::EventSetup& );

private:
   // ----------member data ---------------------------
    edm::InputTag label_;
    std::vector< int > particleID;
    std::vector< double > pMin;
    std::vector< double > ptMin;
    std::vector< double > etaMin;
    std::vector< double > etaMax;
    std::vector< int > status;
    std::vector< double > decayRadiusMin;
    std::vector< double > decayRadiusMax;
    std::vector< double > decayZMin;
    std::vector< double > decayZMax;
};

#endif /*MCSMARTSINGLEGENPARTICLEFILTER*/

