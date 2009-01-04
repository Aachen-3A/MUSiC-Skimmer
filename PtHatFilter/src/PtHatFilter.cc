// -*- C++ -*-
//
// Package:    PtHatFilter
// Class:      PtHatFilter
// 
/**\class PtHatFilter PtHatFilter.cc ePaxDemo/PtHatFilter/src/PtHatFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Carsten Hof
//         Created:  Wed Dec  3 20:07:09 CET 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class PtHatFilter : public edm::EDFilter {
   public:
      explicit PtHatFilter(const edm::ParameterSet&);
      ~PtHatFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      double pt_hat_lower_bound;
      double pt_hat_upper_bound;
      
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PtHatFilter::PtHatFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   // get upper and lower bound from Cfg
   pt_hat_lower_bound = iConfig.getParameter<double>("pt_hat_lower_bound"); 
   pt_hat_upper_bound = iConfig.getParameter<double>("pt_hat_upper_bound"); 
}


PtHatFilter::~PtHatFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------

bool PtHatFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // get and store EventScale aka pt_hat 
   edm::Handle<double> genEventScale;
   iEvent.getByLabel("genEventScale", genEventScale);
   double pthat = *genEventScale;  // pt_hat
  
   if (pthat < pt_hat_lower_bound || pthat >= pt_hat_upper_bound) {
      //std::cout << "pthat = " << pthat << " event rejected!" << std::endl;
      return false; 
   }
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
PtHatFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PtHatFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PtHatFilter);
