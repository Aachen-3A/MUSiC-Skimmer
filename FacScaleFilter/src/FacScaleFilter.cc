// -*- C++ -*-
//
// Package:    FacScaleFilter
// Class:      FacScaleFilter
// 
/**\class FacScaleFilter FacScaleFilter.cc ePaxDemo/FacScaleFilter/src/FacScaleFilter.cc

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
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"

//
// class declaration
//

class FacScaleFilter : public edm::EDFilter {
   public:
      explicit FacScaleFilter(const edm::ParameterSet&);
      ~FacScaleFilter();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      double fac_scale_lower_bound;
      double fac_scale_upper_bound;
      
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
FacScaleFilter::FacScaleFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   // get upper and lower bound from Cfg
   fac_scale_lower_bound = iConfig.getParameter<double>("fac_scale_lower_bound"); 
   fac_scale_upper_bound = iConfig.getParameter<double>("fac_scale_upper_bound"); 
}


FacScaleFilter::~FacScaleFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------

bool FacScaleFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // get factorization scale 
   edm::Handle<reco::PdfInfo> pdfstuff;
   iEvent.getByLabel("genEventPdfInfo", pdfstuff);
   double FacScale = pdfstuff->scalePDF;  // fac_scale
  
   if (FacScale < fac_scale_lower_bound || FacScale >= fac_scale_upper_bound) {
      //std::cout << "FacScale = " << FacScale << " event rejected!" << std::endl;
      return false; 
   }
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
FacScaleFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FacScaleFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(FacScaleFilter);
