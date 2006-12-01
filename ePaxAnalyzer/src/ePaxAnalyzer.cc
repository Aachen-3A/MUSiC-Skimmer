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
// own header file
#include "ePaxDemo/ePaxAnalyzer/interface/ePaxAnalyzer.h"

// system include files
#include <memory>
// include Message Logger for Debug fePaxFile
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// necessary objects:
#include "FWCore/Framework/interface/ESHandle.h"

using namespace std;

//
// constructors and destructor
//
ePaxAnalyzer::ePaxAnalyzer(const edm::ParameterSet& iConfig) {
   //now do what ever initialization is needed
   // The labels used in cfg-file 
   fHepMCLabel = iConfig.getUntrackedParameter<string>("HepMCLabel");
   fJetMCLabel = iConfig.getUntrackedParameter<string>("JetMCLabel");
   fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   fJetRecoLabel = iConfig.getUntrackedParameter<string>("JetRecoLabel");
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");

   fePaxFile.open("ePaxAnalyzer.pxlio");
}

// ------------ MIS Destructor  ------------

ePaxAnalyzer::~ePaxAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
   fePaxFile.close();

}

// ------------ method called to for each event  ------------

void ePaxAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   // set event counter   
   fNumEvt++;    // set this as setUserRecord ?

   // create two ePaxEventViews for Generator/Reconstructed Objects
   ePax::ePaxEventView GenEvtView;
   ePax::ePaxEventView RecEvtView;

   // probably suits better in a method where the Rec/Gen Particles are connected to the Vertices
   ePax::ePaxVertexRef GenVtx = GenEvtView.set().createIndexed<ePax::ePaxVertex>("primary gen vertex");
   ePax::ePaxVertexRef RecVtx = RecEvtView.set().createIndexed<ePax::ePaxVertex>("primary rec vertex");

   // Generator stuff
   analyzeGenInfo(iEvent, GenEvtView);
   analyzeGenJets(iEvent, GenEvtView);
   analyzeGenMET(iEvent, GenEvtView);
   // Reconstructed stuff
   analyzeRecMuons(iEvent, RecEvtView);
   analyzeRecElectrons(iEvent, RecEvtView);
   analyzeRecJets(iEvent, RecEvtView);
   analyzeRecMET(iEvent, RecEvtView);
   analyzeRecGammas(iEvent, RecEvtView);

   GenVtx.printDecayTree();
   RecVtx.printDecayTree();    
   
   fePaxFile.storeObject(GenEvtView);
   fePaxFile.storeObject(RecEvtView);
   fePaxFile.writeEvent();

}
// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   edm::Handle<edm::HepMCProduct> HepMC_Handle ;
   iEvent.getByLabel( fHepMCLabel, HepMC_Handle );

   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;

   const  HepMC::GenEvent* myGenEvent = HepMC_Handle->GetEvent();

   for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
         p != myGenEvent->particles_end(); ++p ) {
      // fill Gen Muons passing some basic cuts
      if ( abs((*p)->pdg_id()) == 13) {
         if ( MuonMC_cuts(p) ) { 
            ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
            part.set().setName("GenMuon");
            part.set().setCharge(((*p)->pdg_id() > 0) ? -1 : 1);
            part.set().vector(ePax::set).setPx((*p)->momentum().x());
            part.set().vector(ePax::set).setPy((*p)->momentum().y());
            part.set().vector(ePax::set).setPz((*p)->momentum().z());
            part.set().vector(ePax::set).setMass((*p)->Mass());
	    part.set().setUserRecord("Vtx_X", (*p)->creationVertex().x());
	    part.set().setUserRecord("Vtx_Y", (*p)->creationVertex().y());
	    part.set().setUserRecord("Vtx_Z", (*p)->creationVertex().z());
	    numMuonMC++; 
         }
      }
      // fill Gen Electrons passing some basic cuts
      if ( abs((*p)->pdg_id()) == 11) {
         if ( EleMC_cuts(p) ) { 
            ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
            part.set().setName("GenEle");
            part.set().setCharge(((*p)->pdg_id() > 0) ? -1 : 1);
            part.set().vector(ePax::set).setPx((*p)->momentum().x());
            part.set().vector(ePax::set).setPy((*p)->momentum().y());
            part.set().vector(ePax::set).setPz((*p)->momentum().z());
            part.set().vector(ePax::set).setMass((*p)->Mass());
	    part.set().setUserRecord("Vtx_X", (*p)->creationVertex().x());
	    part.set().setUserRecord("Vtx_Y", (*p)->creationVertex().y());
	    part.set().setUserRecord("Vtx_Z", (*p)->creationVertex().z());
	    numEleMC++; 
         }
      }
      // fill Gen Gammas passing some basic cuts
      if ( abs((*p)->pdg_id()) == 22) {
         if ( GammaMC_cuts(p) ) { 
            ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
            part.set().setName("GenGamma");
            part.set().setCharge(0);
            part.set().vector(ePax::set).setPx((*p)->momentum().x());
            part.set().vector(ePax::set).setPy((*p)->momentum().y());
            part.set().vector(ePax::set).setPz((*p)->momentum().z());
            part.set().vector(ePax::set).setMass((*p)->Mass());
	    numGammaMC++;
         }
      }
   } //end of loop over generated-particles
   EvtView.set().setUserRecord("NumMuonMC", numMuonMC);
   EvtView.set().setUserRecord("NumEleMC", numEleMC);
   EvtView.set().setUserRecord("NumGammaMC", numGammaMC);
}

// ------------ reading the Generator Jets ------------

void ePaxAnalyzer::analyzeGenJets(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   //Get the GenJet collections
   edm::Handle<reco::GenJetCollection> Kt_GenJets;
//   edm::Handle<reco::GenJetCollection> Icone5_GenJets;
//   edm::Handle<reco::GenJetCollection> Mcone5_GenJets;
//   edm::Handle<reco::GenJetCollection> Mcone7_GenJets;
   iEvent.getByLabel(fJetMCLabel, Kt_GenJets);
//   iEvent.getByLabel(theJetMCLabel, Icone5_GenJets);
//  iEvent.getByLabel(theJetMCLabel, Mcone5_GenJets);
//   iEvent.getByLabel(theJetMCLabel, Mcone7_GenJets);

   int numKtJetsMC = 0;
   
   //Loop over GenJets
   for( reco::GenJetCollection::const_iterator genJet = Kt_GenJets->begin(); 
         genJet != Kt_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) { 
         ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
         part.set().setName("GenKtJets");
         part.set().vector(ePax::set).setPx(genJet->px());
         part.set().vector(ePax::set).setPy(genJet->py());
         part.set().vector(ePax::set).setPz(genJet->pz());
         part.set().vector(ePax::set).setE((genJet->hadEnergy() + genJet->emEnergy()));
	 numKtJetsMC++;
      }
   }   
   EvtView.set().setUserRecord("NumGenKtJets", numKtJetsMC);
}

// ------------ reading the Generator MET ------------

void ePaxAnalyzer::analyzeGenMET(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   edm::Handle<reco::GenMETCollection> GenMet;
   iEvent.getByLabel(fMETMCLabel, GenMet);
   const GenMETCollection *genmetcol = GenMet.product();
   const GenMET genmet = genmetcol->front();  // MET exists only once!
 
   int numMETMC = 0; //means no MET in event

   //there is always MET in event, just decide if cuts passed
   if (METMC_cuts(genmet)) { 
      ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
      part.set().setName("GenMET");
      part.set().vector(ePax::set).setPx(genmet.px());
      part.set().vector(ePax::set).setPy(genmet.py());
      part.set().vector(ePax::set).setPz(0.);
      part.set().vector(ePax::set).setMass(0.);
      part.set().setUserRecord("sumEt", genmet.sumEt());
      part.set().setUserRecord("mEtSig", genmet.mEtSig());
      numMETMC++; 
   }
   EvtView.set().setUserRecord("NumMETMC", numMETMC);
}


// ------------ reading Reconstructed Muons ------------

void ePaxAnalyzer::analyzeRecMuons(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   edm::Handle<reco::MuonCollection> muons;
   iEvent.getByLabel(fMuonRecoLabel, muons);
 
   int numMuonRec = 0;

   for(reco::MuonCollection::const_iterator  muon = muons->begin(); 
         muon != muons->end(); ++muon ) {
      if (Muon_cuts(muon)) { 
         cout << " Found a Rec Muon: \n" 
              << "    pt : " << muon->pt() << endl
	      << "    eta: " << muon->eta() << endl
	      << "    q  : " << muon->charge() << endl;
         ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
         part.set().setName("RecMuon");
         part.set().setCharge(muon->charge());
         part.set().vector(ePax::set).setPx(muon->px());
         part.set().vector(ePax::set).setPy(muon->py());
         part.set().vector(ePax::set).setPz(muon->pz());
         part.set().vector(ePax::set).setE(muon->energy());
         part.set().setUserRecord("Vtx_X", muon->vx());
         part.set().setUserRecord("Vtx_Y", muon->vy());
         part.set().setUserRecord("Vtx_Z", muon->vz()); 
	 // get isolation ;0 
         // part.set().setUserRecord("isolation", 0.9);
         numMuonRec++;
      }
   }
   EvtView.set().setUserRecord("NumMuonRec", numMuonRec);
}

// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   edm::Handle<reco::ElectronCollection> electrons;
   iEvent.getByLabel(fElectronRecoLabel, electrons);

   int numEleRec = 0;
   
   for ( reco::ElectronCollection::const_iterator ele = electrons->begin(); 
            ele != electrons->end(); ++ele ) {
      if (Ele_cuts(ele)) { 
         ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
         part.set().setName("RecElectron");
         part.set().setCharge(ele->charge());
         part.set().vector(ePax::set).setPx(ele->px());
         part.set().vector(ePax::set).setPy(ele->py());
         part.set().vector(ePax::set).setPz(ele->pz());
         part.set().vector(ePax::set).setE(ele->energy());         
         part.set().setUserRecord("Vtx_X", ele->vx());
         part.set().setUserRecord("Vtx_Y", ele->vy());
         part.set().setUserRecord("Vtx_Z", ele->vz());          
	 numEleRec++;
      }
   }
   EvtView.set().setUserRecord("NumEleRec", numEleRec);
}

// ------------ reading Reconstructed Jets ------------
//
// Which kind of Jets?
//

void ePaxAnalyzer::analyzeRecJets(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   edm::Handle<reco::CaloJetCollection> jets;
   iEvent.getByLabel(fJetRecoLabel, jets);
   
   int numJetRec = 0;

   for(reco::CaloJetCollection::const_iterator jet = jets->begin(); 
           jet != jets->end(); ++jet ) {
      if (Jet_cuts(jet)) { 
         ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
         part.set().setName("RecJet");
         part.set().vector(ePax::set).setPx(jet->px());
         part.set().vector(ePax::set).setPy(jet->py());
         part.set().vector(ePax::set).setPz(jet->pz());
         part.set().vector(ePax::set).setE(jet->energy());         
	 numJetRec++;
      }
   }
   EvtView.set().setUserRecord("NumJetRec", numJetRec);
}

// ------------ reading Reconstructed MET ------------

void ePaxAnalyzer::analyzeRecMET(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   //MET produces segmentation violation?!!!
   edm::Handle<reco::CaloMETCollection> CaloMet;
   iEvent.getByLabel(fMETRecoLabel, CaloMet);
   const CaloMETCollection *calometcol = CaloMet.product();
   const CaloMET calomet = calometcol->front();  // MET exists only once!
 
   int numMETRec = 0;
   //there is always MET in event, just decide if cuts passed
   //WHAT IS SIZE OF MET-COLLECTION???
   if (MET_cuts(calomet)) { 
      ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
      part.set().setName("RecMET");
      part.set().vector(ePax::set).setPx(calomet.px());
      part.set().vector(ePax::set).setPy(calomet.py());
      part.set().vector(ePax::set).setPz(0.);
      part.set().vector(ePax::set).setMass(0.);
      part.set().setUserRecord("sumEt", calomet.sumEt());
      part.set().setUserRecord("mEtSig", calomet.mEtSig());
      numMETRec++; 
   }
   EvtView.set().setUserRecord("NumMETRec", numMETRec);
}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, ePax::ePaxEventViewRef EvtView) {

   edm::Handle<reco::PhotonCollection> Photons;
   iEvent.getByLabel(fGammaRecoLabel, Photons);
   
   int numGammaRec = 0;

   for (reco::PhotonCollection::const_iterator photon = Photons->begin(); 
	photon != Photons->end(); ++photon) {  
      if ( Gamma_cuts(photon) ) { 
         ePax::ePaxParticleRef part = EvtView.set().create<ePax::ePaxParticle>();
         part.set().setName("RecGamma");
         part.set().setCharge(0);
         part.set().vector(ePax::set).setPx(photon->px());
         part.set().vector(ePax::set).setPy(photon->py());
         part.set().vector(ePax::set).setPz(photon->pz());
         part.set().vector(ePax::set).setE(photon->energy());
	 numGammaRec++;
      }	 
   }
   EvtView.set().setUserRecord("NumGammaRec", numGammaRec);
}

// ------------ method called once each job just before starting event loop  ------------

void ePaxAnalyzer::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------

void ePaxAnalyzer::endJob() {

   cout << "++++++++++++++++++++++++++++++++++++++" << endl;
   cout << "analyzed " << fNumEvt << " events " << endl;
   
}

// ------------ method to define MC-MUON-cuts

bool ePaxAnalyzer::MuonMC_cuts(HepMC::GenEvent::particle_const_iterator MCmuon) const {
   //
   return ((*MCmuon)->momentum().perp() > 10) ? 1 : 0;
}

// ------------ method to define MC-Electron-cuts

bool ePaxAnalyzer::EleMC_cuts(HepMC::GenEvent::particle_const_iterator MCele) const {
  //
  return ((*MCele)->momentum().perp() > 10) ? 1 : 0;
}

// ------------ method to define MC-Gamma-cuts

bool ePaxAnalyzer::GammaMC_cuts(HepMC::GenEvent::particle_const_iterator MCgamma) const {
  //
  return ((*MCgamma)->momentum().perp() > 10) ? 1 : 0;
}

// ------------ method to define MC-Jet-cuts

bool ePaxAnalyzer::JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const {
  // 
  return (MCjet->pt() > 10) ? 1 : 0;
}

// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::METMC_cuts(const reco::GenMET MCmet) const {
  // 
  return (MCmet.pt() > 30) ? 1 : 0;
}

// ------------ method to define MUON-cuts

bool ePaxAnalyzer::Muon_cuts(reco::MuonCollection::const_iterator muon) const {
   //
   return (muon->pt() > 10) ? 1 : 0;
}

// ------------ method to define ELECTRON-cuts
bool ePaxAnalyzer::Ele_cuts(reco::ElectronCollection::const_iterator ele) const {
   //
   return (ele->pt() > 10) ? 1 : 0;
}

// ------------ method to define JET-cuts
bool ePaxAnalyzer::Jet_cuts(reco::CaloJetCollection::const_iterator jet) const {
   //
   return (jet->pt() > 10) ? 1 : 0;
}

// ------------ method to define GAMMA-cuts
bool ePaxAnalyzer::Gamma_cuts(reco::PhotonCollection::const_iterator photon) const {
   //
   return (photon->energy() > 10) ? 1 : 0;
}

// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::MET_cuts(const reco::MET met) const {
  // 
  return (met.pt() > 30) ? 1 : 0;
}
//define this as a plug-in
DEFINE_FWK_MODULE(ePaxAnalyzer)

/*   edm::Handle<reco::MuonCollection> muons;
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
*/  
