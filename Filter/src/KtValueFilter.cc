// -*- C++ -*-
//
// Package:    KtValueFilter
// Class:      KtValueFilter
// 
/**\class KtValueFilter KtValueFilter.cc ePaxDemo/KtValueFilter/src/KtValueFilter.cc

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

class KtValueFilter : public edm::EDFilter {
   public:
      explicit KtValueFilter(const edm::ParameterSet&);
      ~KtValueFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      double kt_value_lower_bound;
      double kt_value_upper_bound;
      
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
KtValueFilter::KtValueFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   // get upper and lower bound from Cfg
   kt_value_lower_bound = iConfig.getParameter<double>("kt_value_lower_bound"); 
   kt_value_upper_bound = iConfig.getParameter<double>("kt_value_upper_bound"); 
}


KtValueFilter::~KtValueFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------

bool KtValueFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // get ktValue
   edm::Handle<double> genktValue; 
   iEvent.getByLabel("genEventKTValue", genktValue);
//   cout << "the KT = " << *ktValue << std::endl;
   double ktValue = *genktValue;  // pt_hat
  
   if (ktValue < kt_value_lower_bound || ktValue >= kt_value_upper_bound) {
      std::cout << "kt_Value = " << ktValue << " event rejected!" << std::endl;
      return false; 
   }
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
KtValueFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KtValueFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(KtValueFilter);
