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
#include <sstream>
// include Message Logger for Debug
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
   // Get Filename from cfg File
   fFileName = iConfig.getUntrackedParameter<string>("FileName");
   // Debugging
   fDebug = iConfig.getUntrackedParameter<int>("debug");
   // The labels used in cfg-file 
   fHepMCLabel = iConfig.getUntrackedParameter<string>("HepMCLabel");
   fKtJetMCLabel = iConfig.getUntrackedParameter<string>("KtJetMCLabel");
   fItCone5JetMCLabel = iConfig.getUntrackedParameter<string>("ItCone5JetMCLabel");
   fMidCone5JetMCLabel = iConfig.getUntrackedParameter<string>("MidCone5JetMCLabel");
   fMidCone7JetMCLabel = iConfig.getUntrackedParameter<string>("MidCone7JetMCLabel");  
   fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fSAMuonRecoLabel = iConfig.getUntrackedParameter<string>("SAMuonRecoLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   fPixelMatchElectronRecoLabel = iConfig.getUntrackedParameter<string>("PixelMatchElectronRecoLabel");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   fKtJetRecoLabel = iConfig.getUntrackedParameter<string>("KtJetRecoLabel");
   fItCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("ItCone5JetRecoLabel");
   fMidCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("MidCone5JetRecoLabel");
   fMidCone7JetRecoLabel = iConfig.getUntrackedParameter<string>("MidCone7JetRecoLabel");
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");
   
   fNumEvt=0;
   
   fePaxFile.open(fFileName);
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
   ePaxEventView GenEvtView;
   ePaxEventView RecEvtView;

   // probably suits better in a method where the Rec/Gen Particles are connected to the Vertices

   // Do we need that vertex here? Maybe in the "real" analyzer 

//   ePaxVertexRef GenVtx = GenEvtView.set().createIndexed<ePaxVertex>("primary gen vertex");
//   ePaxVertexRef RecVtx = RecEvtView.set().createIndexed<ePaxVertex>("primary rec vertex");

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
   
//   cout << "GenVtx Decay Tree" << endl;
//   GenVtx.printDecayTree();
//   cout << "RecVtx Decay Tree" << endl;
//   RecVtx.printDecayTree();    
   
   fePaxFile.storeObject(GenEvtView);
   fePaxFile.storeObject(RecEvtView);
   if (fDebug > 0) {
      cout << "UserRecord  " <<  "   e   " << "  mu   " << " Gamma " << " KtJet " << "  MET  " << endl;
      cout << "Found (MC): " << setw(4) << GenEvtView.get().findUserRecord("NumEleMC") 
           << setw(7) << GenEvtView.get().findUserRecord("NumMuonMC")
           << setw(7) << GenEvtView.get().findUserRecord("NumGammaMC") 
           << setw(4) << GenEvtView.get().findUserRecord("NumKtJetMC") << "/" 
           << GenEvtView.get().findUserRecord("NumItCone5JetMC") << "/" 
           << GenEvtView.get().findUserRecord("NumMidCone5JetMC") << "/" 
           << GenEvtView.get().findUserRecord("NumMidCone7JetMC")  
           << setw(7) << GenEvtView.get().findUserRecord("NumMETMC") << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView.get().findUserRecord("NumEleRec")    
           << setw(7) << RecEvtView.get().findUserRecord("NumMuonRec") 
           << setw(7) << RecEvtView.get().findUserRecord("NumGammaRec")
           << setw(4) << RecEvtView.get().findUserRecord("NumKtJetRec") << "/"
           << RecEvtView.get().findUserRecord("NumItCone5JetRec") << "/"
           << RecEvtView.get().findUserRecord("NumMidCone5JetRec") << "/"
           << RecEvtView.get().findUserRecord("NumMidCone7JetRec") 
           << setw(7) << RecEvtView.get().findUserRecord("NumMETRec", 0.) << endl;
   }
   ostringstream EventType;
   EventType << GenEvtView.get().findUserRecord("NumEleMC") <<  "e"
             <<  GenEvtView.get().findUserRecord("NumMuonMC") << "mu"
             << GenEvtView.get().findUserRecord("NumGammaMC") << "gam"
             << GenEvtView.get().findUserRecord("NumKtJetMC") << "kt"
             << GenEvtView.get().findUserRecord("NumMETMC") << "met";
   if (fDebug > 0) cout << "Event Type: " << EventType.str() << endl;
   fePaxFile.writeEvent(EventType.str());

}
// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   edm::Handle<edm::HepMCProduct> HepMC_Handle ;
   iEvent.getByLabel( fHepMCLabel, HepMC_Handle );

   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;

   const  HepMC::GenEvent* myGenEvent = HepMC_Handle->GetEvent();

   for ( HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin();
         p != myGenEvent->particles_end(); ++p ) {
      // fill Gen Muons passing some basic cuts
      if ( abs((*p)->pdg_id()) == 13 && (*p)->status() == 1 && abs((*p)->momentum().eta()) < 2.1) {
         if ( MuonMC_cuts(p) ) { 
            ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
            part.set().setName("GenMuon");
            part.set().setCharge(((*p)->pdg_id() > 0) ? -1 : 1);
            part.set().vector(ePaxSet).setPx((*p)->momentum().x());
            part.set().vector(ePaxSet).setPy((*p)->momentum().y());
            part.set().vector(ePaxSet).setPz((*p)->momentum().z());
            part.set().vector(ePaxSet).setMass((*p)->Mass());
	    part.set().setUserRecord("Vtx_X", (*p)->creationVertex().x());
	    part.set().setUserRecord("Vtx_Y", (*p)->creationVertex().y());
	    part.set().setUserRecord("Vtx_Z", (*p)->creationVertex().z());
	    numMuonMC++; 
         }
      }
      // fill Gen Electrons passing some basic cuts
      if ( abs((*p)->pdg_id()) == 11 && (*p)->status() == 1 && abs((*p)->momentum().eta()) < 2.5) {
         if ( EleMC_cuts(p) ) { 
            ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
            part.set().setName("GenEle");
            part.set().setCharge(((*p)->pdg_id() > 0) ? -1 : 1);
            part.set().vector(ePaxSet).setPx((*p)->momentum().x());
            part.set().vector(ePaxSet).setPy((*p)->momentum().y());
            part.set().vector(ePaxSet).setPz((*p)->momentum().z());
            part.set().vector(ePaxSet).setMass((*p)->Mass());
	    part.set().setUserRecord("Vtx_X", (*p)->creationVertex().x());
	    part.set().setUserRecord("Vtx_Y", (*p)->creationVertex().y());
	    part.set().setUserRecord("Vtx_Z", (*p)->creationVertex().z());
	    numEleMC++; 
         }
      }
      // fill Gen Gammas passing some basic cuts
      if ( abs((*p)->pdg_id()) == 22 && (*p)->status() == 1 && abs((*p)->momentum().eta()) < 2.5) {
         if ( GammaMC_cuts(p) ) { 
            ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
            part.set().setName("GenGamma");
            part.set().setCharge(0);
            part.set().vector(ePaxSet).setPx((*p)->momentum().x());
            part.set().vector(ePaxSet).setPy((*p)->momentum().y());
            part.set().vector(ePaxSet).setPz((*p)->momentum().z());
            part.set().vector(ePaxSet).setMass((*p)->Mass());
	    numGammaMC++;
         }
      }
   } //end of loop over generated-particles
   if (fDebug > 1)  cout << "MC Found:  " << numMuonMC <<  " mu  " << numEleMC << " e  " << numGammaMC << " gamma" << endl;
   EvtView.set().setUserRecord("NumMuonMC", numMuonMC);
   EvtView.set().setUserRecord("NumEleMC", numEleMC);
   EvtView.set().setUserRecord("NumGammaMC", numGammaMC);
}

// ------------ reading the Generator Jets ------------

void ePaxAnalyzer::analyzeGenJets(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   //Get the GenJet collections
   edm::Handle<reco::GenJetCollection> Kt_GenJets;
   edm::Handle<reco::GenJetCollection> ItCone5_GenJets;
   edm::Handle<reco::GenJetCollection> MidCone5_GenJets;
   edm::Handle<reco::GenJetCollection> MidCone7_GenJets;
   iEvent.getByLabel(fKtJetMCLabel, Kt_GenJets);
   iEvent.getByLabel(fItCone5JetMCLabel, ItCone5_GenJets);
   iEvent.getByLabel(fMidCone5JetMCLabel, MidCone5_GenJets);
   iEvent.getByLabel(fMidCone7JetMCLabel, MidCone7_GenJets);

   int numKtJetMC = 0;
   int numItCone5JetMC = 0;
   int numMidCone5JetMC = 0;
   int numMidCone7JetMC = 0;
   
   //Loop over KtGenJets
   for( reco::GenJetCollection::const_iterator genJet = Kt_GenJets->begin(); 
         genJet != Kt_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) { 
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("GenKtJets");
         part.set().vector(ePaxSet).setPx(genJet->px());
         part.set().vector(ePaxSet).setPy(genJet->py());
         part.set().vector(ePaxSet).setPz(genJet->pz());
         part.set().vector(ePaxSet).setE((genJet->hadEnergy() + genJet->emEnergy()));
	 numKtJetMC++;
      }
   }   
   EvtView.set().setUserRecord("NumKtJetMC", numKtJetMC);

   //Loop over ItCone5GenJets
   for( reco::GenJetCollection::const_iterator genJet = ItCone5_GenJets->begin();
         genJet != ItCone5_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("GenItCone5Jets");
         part.set().vector(ePaxSet).setPx(genJet->px());
         part.set().vector(ePaxSet).setPy(genJet->py());
         part.set().vector(ePaxSet).setPz(genJet->pz());
         part.set().vector(ePaxSet).setE((genJet->hadEnergy() + genJet->emEnergy()));
         numItCone5JetMC++;
      }
   }
   EvtView.set().setUserRecord("NumItCone5JetMC", numItCone5JetMC);

   //Loop over MidCone5GenJets
   for( reco::GenJetCollection::const_iterator genJet = MidCone5_GenJets->begin();
         genJet != MidCone5_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("GenMidCone5Jets");
         part.set().vector(ePaxSet).setPx(genJet->px());
         part.set().vector(ePaxSet).setPy(genJet->py());
         part.set().vector(ePaxSet).setPz(genJet->pz());
         part.set().vector(ePaxSet).setE((genJet->hadEnergy() + genJet->emEnergy()));
         numMidCone7JetMC++;
      }
   }
   EvtView.set().setUserRecord("NumMidCone5JetMC", numMidCone5JetMC);

   //Loop over MidCone7GenJets
   for( reco::GenJetCollection::const_iterator genJet = MidCone7_GenJets->begin();
         genJet != MidCone7_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("GenMidCone7Jets");
         part.set().vector(ePaxSet).setPx(genJet->px());
         part.set().vector(ePaxSet).setPy(genJet->py());
         part.set().vector(ePaxSet).setPz(genJet->pz());
         part.set().vector(ePaxSet).setE((genJet->hadEnergy() + genJet->emEnergy()));
         numMidCone7JetMC++;
      }
   }
   EvtView.set().setUserRecord("NumMidCone7JetMC", numMidCone5JetMC);
  
   if (fDebug > 1) cout << "Found MC Jets:  " << numKtJetMC << " Kt  " << numItCone5JetMC << " It5  " 
        << numMidCone5JetMC << " Mid5 " << numMidCone7JetMC << " Mid7 " << endl;
   
}

// ------------ reading the Generator MET ------------

void ePaxAnalyzer::analyzeGenMET(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   edm::Handle<reco::GenMETCollection> GenMet;
   iEvent.getByLabel(fMETMCLabel, GenMet);
   const GenMETCollection *genmetcol = GenMet.product();
   const GenMET genmet = genmetcol->front();  // MET exists only once!
 
   int numMETMC = 0; //means no MET in event

   ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
   part.set().setName("GenMET");
   part.set().vector(ePaxSet).setPx(genmet.px());
   part.set().vector(ePaxSet).setPy(genmet.py());
   part.set().vector(ePaxSet).setPz(0.);
   part.set().vector(ePaxSet).setMass(0.);
   part.set().setUserRecord("sumEt", genmet.sumEt());
   part.set().setUserRecord("mEtSig", genmet.mEtSig());

   if (fDebug > 1) cout << "GenMET before muon corr: Px = " << genmet.px() << "   Py = " << genmet.py() << "   Pt = " << part.get().vector().getPt() << endl;
   // Perform Muon Corrections!
   // loop over muons and subtract them
   if (EvtView.get().findUserRecord("NumMuonMC") > 0) { 
      for (pxl::Objects::TypeIterator<pxl::Particle> iter(EvtView().getObjects()); !iter.isDone(); iter.next()) { 
         if (iter.object().get().getName() == "GenMuon")  {
            if (fDebug > 1) cout << "Correcting with " << iter.object().get().getName() << " px = " << iter.object().get().vector().getPx() 
                                 << " Py = " << iter.object().get().vector().getPy() << endl;
            part.set() -= iter.object().get();
         }
      }
   } 
   if (fDebug) cout << "GenMET after muon corr: Px = " << part.get().vector().getPx() << "   Py = " << part.get().vector().getPy() << "   Pt = " << part.get().vector().getPt() << endl;     
   //there is always MET in event, just decide if cuts passed
   if (part.get().vector().getPt() > 20) {
      numMETMC++; 
   }
   EvtView.set().setUserRecord("NumMETMC", numMETMC);
   if (numMETMC && fDebug > 1) cout << "Event contains MET" << endl;
}


// ------------ reading Reconstructed Muons ------------

void ePaxAnalyzer::analyzeRecMuons(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   edm::Handle<reco::MuonCollection> muons;
//   edm::Handle<reco::MuonCollection> SAmuons;
   iEvent.getByLabel(fMuonRecoLabel, muons);
//   iEvent.getByLabel(fSAMuonRecoLabel, SAmuons);
 
   int numMuonRec = 0;
//   int numSAMuonRec = 0;

   for(reco::MuonCollection::const_iterator  muon = muons->begin(); 
         muon != muons->end(); ++muon ) {
      if (Muon_cuts(muon)) { 
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecMuon");
         part.set().setCharge(muon->charge());
         part.set().vector(ePaxSet).setPx(muon->px());
         part.set().vector(ePaxSet).setPy(muon->py());
         part.set().vector(ePaxSet).setPz(muon->pz());
         part.set().vector(ePaxSet).setE(muon->energy());
         part.set().setUserRecord("Vtx_X", muon->vx());
         part.set().setUserRecord("Vtx_Y", muon->vy());
         part.set().setUserRecord("Vtx_Z", muon->vz()); 
	 // get isolation ;0 
         // part.set().setUserRecord("isolation", 0.9);
         numMuonRec++;
      }
   }
   EvtView.set().setUserRecord("NumMuonRec", numMuonRec);
  
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 
  
/*   for(reco::MuonCollection::const_iterator  muon = SAmuons->begin();
         muon != SAmuons->end(); ++muon ) {
      if (Muon_cuts(muon)) {
         cout << " Found a Rec Muon: \n"
              << "    pt : " << muon->pt() << endl
              << "    eta: " << muon->eta() << endl
              << "    q  : " << muon->charge() << endl;
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecSAMuon");
         part.set().setCharge(muon->charge());
         part.set().vector(ePaxSet).setPx(muon->px());
         part.set().vector(ePaxSet).setPy(muon->py());
         part.set().vector(ePaxSet).setPz(muon->pz());
         part.set().vector(ePaxSet).setE(muon->energy());
         part.set().setUserRecord("Vtx_X", muon->vx());
         part.set().setUserRecord("Vtx_Y", muon->vy());
         part.set().setUserRecord("Vtx_Z", muon->vz());
         // get isolation ;0 
         // part.set().setUserRecord("isolation", 0.9);
         numSAMuonRec++;
      }
   }
   EvtView.set().setUserRecord("NumSAMuonRec", numSAMuonRec);
*/
}

// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   edm::Handle<reco::ElectronCollection> electrons;
//   edm::Handle<reco::ElectronCollection> pixelelectrons;
   iEvent.getByLabel(fElectronRecoLabel, electrons);
//   iEvent.getByLabel(fPixelMatchElectronRecoLabel, pixelelectrons);

   int numEleRec = 0;
//   int numPixelEleRec = 0;   

   for ( reco::ElectronCollection::const_iterator ele = electrons->begin(); 
            ele != electrons->end(); ++ele ) {
      if (Ele_cuts(ele)) { 
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecElectron");
         part.set().setCharge(ele->charge());
         part.set().vector(ePaxSet).setPx(ele->px());
         part.set().vector(ePaxSet).setPy(ele->py());
         part.set().vector(ePaxSet).setPz(ele->pz());
         part.set().vector(ePaxSet).setE(ele->energy());         
         part.set().setUserRecord("Vtx_X", ele->vx());
         part.set().setUserRecord("Vtx_Y", ele->vy());
         part.set().setUserRecord("Vtx_Z", ele->vz());          
	 numEleRec++;
      }
   }
   EvtView.set().setUserRecord("NumEleRec", numEleRec);

   if (fDebug > 1) cout << "RecEle:  " << numEleRec << endl;

/*   for ( reco::ElectronCollection::const_iterator ele = pixelelectrons->begin();
            ele != pixelelectrons->end(); ++ele ) {
      if (Ele_cuts(ele)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecPixelElectron");
         part.set().setCharge(ele->charge());
         part.set().vector(ePaxSet).setPx(ele->px());
         part.set().vector(ePaxSet).setPy(ele->py());
         part.set().vector(ePaxSet).setPz(ele->pz());
         part.set().vector(ePaxSet).setE(ele->energy());
         part.set().setUserRecord("Vtx_X", ele->vx());
         part.set().setUserRecord("Vtx_Y", ele->vy());
         part.set().setUserRecord("Vtx_Z", ele->vz());
         numPixelEleRec++;
      }
   }
   EvtView.set().setUserRecord("NumPixelEleRec", numPixelEleRec);
*/
}

// ------------ reading Reconstructed Jets ------------
//
// Which kind of Jets?
//

void ePaxAnalyzer::analyzeRecJets(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   edm::Handle<reco::CaloJetCollection> Ktjets;
   edm::Handle<reco::CaloJetCollection> ItCone5jets;
   edm::Handle<reco::CaloJetCollection> MidCone5jets;
   edm::Handle<reco::CaloJetCollection> MidCone7jets;

   iEvent.getByLabel(fKtJetRecoLabel, Ktjets);
   iEvent.getByLabel(fItCone5JetRecoLabel, ItCone5jets);
   iEvent.getByLabel(fMidCone5JetRecoLabel, MidCone5jets);
   iEvent.getByLabel(fMidCone7JetRecoLabel, MidCone7jets);
   
   int numKtJetRec = 0;
   int numItCone5JetRec = 0;
   int numMidCone5JetRec = 0;
   int numMidCone7JetRec = 0;

   for(reco::CaloJetCollection::const_iterator jet = Ktjets->begin(); 
           jet != Ktjets->end(); ++jet ) {
      if (Jet_cuts(jet)) { 
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecKtJet");
         part.set().vector(ePaxSet).setPx(jet->px());
         part.set().vector(ePaxSet).setPy(jet->py());
         part.set().vector(ePaxSet).setPz(jet->pz());
         part.set().vector(ePaxSet).setE(jet->energy());         
	 numKtJetRec++;
      }
   }
   EvtView.set().setUserRecord("NumKtJetRec", numKtJetRec);

   for(reco::CaloJetCollection::const_iterator jet = ItCone5jets->begin();
           jet != ItCone5jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecItConeJet");
         part.set().vector(ePaxSet).setPx(jet->px());
         part.set().vector(ePaxSet).setPy(jet->py());
         part.set().vector(ePaxSet).setPz(jet->pz());
         part.set().vector(ePaxSet).setE(jet->energy());
         numItCone5JetRec++;
      }
   }
   EvtView.set().setUserRecord("NumItCone5JetRec", numItCone5JetRec);

   for(reco::CaloJetCollection::const_iterator jet = MidCone5jets->begin();
           jet != MidCone5jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecMidCone5Jet");
         part.set().vector(ePaxSet).setPx(jet->px());
         part.set().vector(ePaxSet).setPy(jet->py());
         part.set().vector(ePaxSet).setPz(jet->pz());
         part.set().vector(ePaxSet).setE(jet->energy());
         numMidCone5JetRec++;
      }
   }
   EvtView.set().setUserRecord("NumMidCone5JetRec", numMidCone5JetRec);

   for(reco::CaloJetCollection::const_iterator jet = MidCone7jets->begin();
           jet != MidCone7jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecMidCone7Jet");
         part.set().vector(ePaxSet).setPx(jet->px());
         part.set().vector(ePaxSet).setPy(jet->py());
         part.set().vector(ePaxSet).setPz(jet->pz());
         part.set().vector(ePaxSet).setE(jet->energy());
         numMidCone7JetRec++;
      }
   }
   EvtView.set().setUserRecord("NumMidCone7JetRec", numMidCone7JetRec);

   if (fDebug > 1) cout << "Found Rec Jets:  " << numKtJetRec << " Kt  " << numItCone5JetRec << " It5  "
                        << numMidCone5JetRec << " Mid5 " << numMidCone7JetRec << " Mid7 " << endl;

}

// ------------ reading Reconstructed MET ------------

void ePaxAnalyzer::analyzeRecMET(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   //MET produces segmentation violation?!!!
   edm::Handle<reco::CaloMETCollection> CaloMet;
   iEvent.getByLabel(fMETRecoLabel, CaloMet);
   const CaloMETCollection *calometcol = CaloMet.product();
   const CaloMET calomet = calometcol->front();  // MET exists only once!
 
   int numMETRec = 0;
   ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
   part.set().setName("RecMET");
   part.set().vector(ePaxSet).setPx(calomet.px());
   part.set().vector(ePaxSet).setPy(calomet.py());
   part.set().vector(ePaxSet).setPz(0.);
   part.set().vector(ePaxSet).setMass(0.);
   part.set().setUserRecord("sumEt", calomet.sumEt());
   part.set().setUserRecord("mEtSig", calomet.mEtSig());
   
   if (fDebug > 1) cout << "RecMET before muon corr: Px = " << calomet.px() << "   Py = " << calomet.py() << "   Pt = " << part.get().vector().getPt() << endl;   
   cout << " MET (uncorr): ( " << part.get().vector().getPx() << ", " << part.get().vector().getPy() << ", "
                       << part.get().vector().getPz() <<  " )    E: " << part.get().vector().getE() << "   Et = "
                       << part.get().vector().getEt() << " Theta: " << part.get().vector().getTheta()
                       << " Mass: " << part.get().vector().getMass() << endl;
   // Perform Muon Corrections!   
   // loop over muons and subtract them   
   if (EvtView.get().findUserRecord("NumMuonRec") > 0) {      
   for (pxl::Objects::TypeIterator<pxl::Particle> iter(EvtView().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().get().getName() == "RecMuon")  {
         if (fDebug > 1) cout << "Correcting with " << iter.object().get().getName() << " px = " << iter.object().get().vector().getPx() 
                              << " Py = " << iter.object().get().vector().getPy() << endl; 
            part.set() -= iter.object().get(); 
         }  
      }   
   } 
   part.set().vector(ePaxSet).setPz(0.);  
   part.set().vector(ePaxSet).setMass(0.);
   cout << " MET (corr): ( " << part.get().vector().getPx() << ", " << part.get().vector().getPy() << ", "
                       << part.get().vector().getPz() <<  " )    E: " << part.get().vector().getE() << "   Et = "
                       << part.get().vector().getEt() << " Theta: " << part.get().vector().getTheta()
                       << " Mass: " << part.get().vector().getMass() << endl;
   if (fDebug) cout << "RecMET after muon corr: Px = " << part.get().vector().getPx() << "   Py = " << part.get().vector().getPy() 
                    << "   Pt = " << part.get().vector().getPt() << endl;   
   
   if (part.get().vector().getPt() > 20) {
      numMETRec++;
   }

   EvtView.set().setUserRecord("NumMETRec", numMETRec);
   if (numMETRec && fDebug > 1) cout << "Found RecMET" << endl;
}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, ePaxEventViewRef EvtView) {

   edm::Handle<reco::PhotonCollection> Photons;
   iEvent.getByLabel(fGammaRecoLabel, Photons);
   
   int numGammaRec = 0;

   for (reco::PhotonCollection::const_iterator photon = Photons->begin(); 
	photon != Photons->end(); ++photon) {  
      if ( Gamma_cuts(photon) ) { 
         ePaxParticleRef part = EvtView.set().create<ePaxParticle>();
         part.set().setName("RecGamma");
         part.set().setCharge(0);
         part.set().vector(ePaxSet).setPx(photon->px());
         part.set().vector(ePaxSet).setPy(photon->py());
         part.set().vector(ePaxSet).setPz(photon->pz());
         part.set().vector(ePaxSet).setE(photon->energy());
	 numGammaRec++;
      }	 
   }
   EvtView.set().setUserRecord("NumGammaRec", numGammaRec);
   if (fDebug > 1) cout << "Rec Gamma: " << numGammaRec << endl;
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
  return ((*MCgamma)->momentum().perp() > 30) ? 1 : 0;
}

// ------------ method to define MC-Jet-cuts

bool ePaxAnalyzer::JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const {
  // 
  return (MCjet->pt() > 30) ? 1 : 0;
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
   return (jet->pt() > 30) ? 1 : 0;
}

// ------------ method to define GAMMA-cuts
bool ePaxAnalyzer::Gamma_cuts(reco::PhotonCollection::const_iterator photon) const {
   //
   return (photon->energy() > 30) ? 1 : 0;
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
      pa.set().vector(ePaxSet).setPx(muon->px());
      pa.set().vector(ePaxSet).setPy(muon->py());
      pa.set().vector(ePaxSet).setPz(muon->pz());
      pa.set().vector(ePaxSet).setE(muon->energy());
      pa.set().setUserRecord("isolation", 0.9);
      
      //pa.get().findUserRecord("isolation"); // kracht"s wenns ihn nicht gibt
      //pa.get().findUserRecord("isolation", 0.); // default wenns ihn nicht gibt
      

//       cout << " Found a Rec Muon: \n" 
//            << "    pt : " << muon->pt() << endl
// 	   << "    eta: " << muon->eta() << endl
// 	   << "    q  : " << muon->charge() << endl;
   }
*/  
