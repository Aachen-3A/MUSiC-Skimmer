// -*- C++ -*-
//
// Package:    ePaxAnalyzer
// Class:      ePaxAnalyzer
// 
/**\class ePaxAnalyzer ePaxAnalyzer.cc PaxDemo/ePaxAnalyzer/src/ePaxAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  
//         Created:  Mo Okt 30 12:03:52 CET 2006
// $Id$
//
//


// system include files
#include <memory>
// include Message Logger for Debug Output
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Access to Generator stuff
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Pax stuff
#include "ePaxPxl/ePax/interface/ePAX.h"

using namespace edm;
using namespace std;
using namespace ePax;

iotl__declareSpyTypeExplicit(reco::Muon, "reco::Muon", reco__Muon, ptr, {}, {})
typedef const reco::Muon CmsTestClass;

// possible user typedef:
typedef ptl::SpyObject<CmsTestClass> myCmsTestClass;
typedef ptl::SpyWkPtr<CmsTestClass> myCmsTestClassWkPtr;
typedef myCmsTestClass& myCmsTestClassRef;

//
// class declaration ==> should be in the Header File ...
//

class ePaxAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ePaxAnalyzer(const edm::ParameterSet&);
      ~ePaxAnalyzer();


   private:
      virtual void beginJob(const EventSetup&) ;
      virtual void analyze(const Event&, const EventSetup&);
      virtual void endJob() ;
      
      // split analysis in logical parts
      
      virtual void analyzeGenInfo(const Event&, const EventSetup&);
      virtual void analyzeGenJets(const Event&, const EventSetup&);
      virtual void analyzeGenMET(const Event&, const EventSetup&);
      
      virtual void analyzeSimInfo(const Event&, const EventSetup&);
      
      virtual void analyzeRecMuons(const Event&, const EventSetup&);
      virtual void analyzeRecElectrons(const Event&, const EventSetup&);
      virtual void analyzeRecJets(const Event&, const EventSetup&);
      virtual void analyzeRecMET(const Event&, const EventSetup&);
      virtual void analyzeRecGammas(const Event&, const EventSetup&);
      virtual void analyzeRecbJets(const Event&, const EventSetup&);
      
      // ----------member data ---------------------------
      string fMuonRecoLabel;
      iotl::oDiskFile output;

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
ePaxAnalyzer::ePaxAnalyzer(const ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   
   output.open("ePaxAnalyzer.pxlio");
}


ePaxAnalyzer::~ePaxAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
   output.close();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void ePaxAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup) {
   
   // own code
   analyzeGenInfo(iEvent, iSetup);

   ePaxEventView evRec;
   ePaxEventView evGen;
   
   ePaxVertexRef vx = evRec.set().createIndexed<ePaxVertex>("primary vertex");


   edm::Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(fMuonRecoLabel, muons);

 

   cout << "We got " << muons->size() << " Muons" << endl;
   for(reco::MuonCollection::const_iterator  muon = muons->begin(); 
         muon != muons->end(); ++muon ) {
	 
//    ptl::SpyObject<CmsTestClass>& mu = ev.set().create<ptl::SpyObject<CmsTestClass> >(&(*muon));
//    mu.linkMother(vx);
	 
      ePaxParticleRef pa = evRec.set().create<ePaxParticle>();
      pa.linkMother(vx);
      pa.set().setName("muon");
      pa.set().setCharge(muon->charge());
      pa.set().vector(ePax::set).setPx(muon->px());
      pa.set().vector(ePax::set).setPy(muon->py());
      pa.set().vector(ePax::set).setPz(muon->pz());
      pa.set().vector(ePax::set).setE(muon->energy());
      pa.set().setUserRecord("isolation", 0.9);
      
      //pa.get().findUserRecord("isolation"); // kracht"s wenns ihn nicht gibt
      //pa.get().findUserRecord("isolation", 0.); // default wenns ihn nicht gibt
      

//       cout << " Found a Rec Muon: \n" 
//            << "    pt : " << muon->pt() << endl
// 	   << "    eta: " << muon->eta() << endl
// 	   << "    q  : " << muon->charge() << endl;
   }
  
   vx.printDecayTree();
   
   output.storeObject(evRec);
   output.storeObject(evGen);
   output.writeEvent();

}
// ------------ reading the Generator Stuff ------------

//void ePaxAnalyzer::analyzeGenInfo(const Event& iEvent, const EventSetup& iSetup, ePaxEventViewRef ev) {
void ePaxAnalyzer::analyzeGenInfo(const Event& iEvent, const EventSetup& iSetup) {

   Handle< HepMCProduct > HepMC_Handle ;
   //iEvent.getByLabel( "source", EvtHandle ) ;
   // same as source but with a vertex smeared inside the CMS Detector interaction zone
   iEvent.getByLabel( "VtxSmeared", HepMC_Handle );

   const  HepMC::GenEvent* myGenEvent = HepMC_Handle->GetEvent();
   
   //myGenEvent->print();
   
   // Iterate over all particles in the Event
   for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p ) {
      if ( abs((*p)->pdg_id()) == 13) {
         cout << "Found a MC Muon: " << endl;
         (*p)->print();
      }
   }
}


// ------------ reading Reconstructed Muons ------------

void ePaxAnalyzer::analyzeRecMuons(const Event& iEvent, const EventSetup& iSetup) {

   edm::Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(fMuonRecoLabel, muons);

   for(reco::MuonCollection::const_iterator  muon = muons->begin(); 
         muon != muons->end(); ++muon ) {
      cout << " Found a Rec Muon: \n" 
           << "    pt : " << muon->pt() << endl
	   << "    eta: " << muon->eta() << endl
	   << "    q  : " << muon->charge() << endl;
      /*fMuonRec_Px[fNumMuonRec] = muon->px();
      fMuonRec_Py[fNumMuonRec] = muon->py();
      fMuonRec_Pz[fNumMuonRec] = muon->pz();
      fMuonRec_Pt[fNumMuonRec] = muon->pt();
      fMuonRec_E[fNumMuonRec] = muon->energy();
      fMuonRec_Phi[fNumMuonRec] = muon->phi();
      fMuonRec_Theta[fNumMuonRec] = muon->theta();
      fMuonRec_Eta[fNumMuonRec] = muon->eta();
      fMuonRec_Charge[fNumMuonRec] = muon->charge();
      fMuonRec_Vtx_X[fNumMuonRec] = muon->vx();
      fMuonRec_Vtx_Y[fNumMuonRec] = muon->vy();
      fMuonRec_Vtx_Z[fNumMuonRec] = muon->vz();
      fNumMuonRec++;*/
   }
}



// ------------ reading Reconstructed Electrons ------------
// 
// in ORCA there were only candidates, that requiered some probabilistic discrimination 
// what about CMSSW?
//

void ePaxAnalyzer::analyzeRecElectrons(const Event& iEvent, const EventSetup& iSetup) {


}

// ------------ reading Reconstructed Jets ------------
//
// Which kind of Jets?
//

void ePaxAnalyzer::analyzeRecJets(const Event& iEvent, const EventSetup& iSetup) {


}

// ------------ reading Reconstructed MET ------------

void ePaxAnalyzer::analyzeRecMET(const Event& iEvent, const EventSetup& iSetup) {


}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const Event& iEvent, const EventSetup& iSetup) {


}

// ------------ reading Reconstructed b-Jets ------------

void ePaxAnalyzer::analyzeRecbJets(const Event& iEvent, const EventSetup& iSetup) {


}
// ------------ reading the Generator Jets ------------

void ePaxAnalyzer::analyzeGenJets(const Event& iEvent, const EventSetup& iSetup) {

//++recoGenJets "GenJetIcone5" ""
//++recoGenJets "GenJetKt" ""
//++recoGenJets "GenJetMcone5" ""
//++recoGenJets "GenJetMcone7" ""

/*   Handle< GenJetCollection > GenJetIcone5_Handle ;
   Handle< GenJetCollection > GenJetKt_Handle ;
   Handle< GenJetCollection > GenJetMcone5_Handle ;
   Handle< GenJetCollection > GenJetMcone7_Handle ;
   iEvent.getByLabel( "GenJetIcone5", GenJetIcone5_Handle );
   iEvent.getByLabel( "GenJetKt", GenJetKt_Handle );
   iEvent.getByLabel( "GenJetMcone5", GenJetMcone5_Handle );
   iEvent.getByLabel( "GenJetMcone7", GenJetMcone7_Handle );
   
   for( GenJetCollection::const_iterator jetItr = genJets->begin() ; jetItr != genJets->end() ;	 ++jetItr ) {
   
   }*/


}

// ------------ reading the Generator MET ------------

void ePaxAnalyzer::analyzeGenMET(const Event& iEvent, const EventSetup& iSetup) {


}

// ------------ reading the GEANT Stuff ------------

void ePaxAnalyzer::analyzeSimInfo(const Event& iEvent, const EventSetup& iSetup) {

}

// ------------ method called once each job just before starting event loop  ------------

void ePaxAnalyzer::beginJob(const EventSetup& iSetup) {
}

// ------------ method called once each job just after ending the event loop  ------------

void ePaxAnalyzer::endJob() {

}


//define this as a plug-in
DEFINE_FWK_MODULE(ePaxAnalyzer)
