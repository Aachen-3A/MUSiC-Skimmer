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
// for electron shapes:
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
// for ECAL enumerator
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
// for HCAL navigation used in HadOverEm calculation
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

//for GenPartciles
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

//for isolation
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/TrackReco/interface/Track.h"

using namespace std;

//
// constructors and destructor
//
ePaxAnalyzer::ePaxAnalyzer(const edm::ParameterSet& iConfig) {
   //now do what ever initialization is needed
   // Get Filename from cfg File
   fFileName = iConfig.getUntrackedParameter<string>("FileName");
   // Get Physics process
   fProcess = iConfig.getUntrackedParameter<string>("Process");
   // Gen-Only or also Rec-information
   fGenOnly = iConfig.getUntrackedParameter<bool>("GenOnly");
   // Debugging
   fDebug = iConfig.getUntrackedParameter<int>("debug");
   // The labels used in cfg-file 
   fHepMCLabel = iConfig.getUntrackedParameter<string>("HepMCLabel");
   fgenParticleCandidatesLabel  = iConfig.getUntrackedParameter<string>("genParticleCandidatesLabel");
   fKtJetMCLabel = iConfig.getUntrackedParameter<string>("KtJetMCLabel");
   fItCone5JetMCLabel = iConfig.getUntrackedParameter<string>("ItCone5JetMCLabel");
   fMidCone5JetMCLabel = iConfig.getUntrackedParameter<string>("MidCone5JetMCLabel");
   fMidCone7JetMCLabel = iConfig.getUntrackedParameter<string>("MidCone7JetMCLabel");  
   fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fVertexRecoLabel = iConfig.getUntrackedParameter<string>("VertexRecoLabel");
   fSAMuonRecoLabel = iConfig.getUntrackedParameter<string>("SAMuonRecoLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   fBarrelClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("barrelClusterShapeAssociation");
   fEndcapClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("endcapClusterShapeAssociation");   
   fPixelMatchElectronRecoLabel = iConfig.getUntrackedParameter<string>("PixelMatchElectronRecoLabel");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   fKtJetRecoLabel = iConfig.getUntrackedParameter<string>("KtJetRecoLabel");
   fItCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("ItCone5JetRecoLabel");
   fMidCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("MidCone5JetRecoLabel");
   fMidCone7JetRecoLabel = iConfig.getUntrackedParameter<string>("MidCone7JetRecoLabel");
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");
   fHBHELabel = iConfig.getUntrackedParameter<string>("fHBHELabel");
   fHBHEInstanceName = iConfig.getUntrackedParameter<string>("fHBHEInstanceName");
   
   
   fNumEvt=0;
   Matcher = new ParticleMatcher();
   
   fePaxFile.open(fFileName);
}

// ------------ MIS Destructor  ------------

ePaxAnalyzer::~ePaxAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
   fePaxFile.close();
   delete Matcher;
}

// ------------ method called to for each event  ------------

void ePaxAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //cout<<"Event Number: "<<fNumEvt<<endl;
   // set event counter   
   fNumEvt++;    // set this as setUserRecord ?
   // get the calorimeter geometry in order to navigate through it for HadOverEm calculation
   iSetup.get<IdealGeometryRecord>().get(theCaloGeom);

   // create two ePaxEventViews for Generator/Reconstructed Objects
   pxl::EventView GenEvtView;
   pxl::EventView RecEvtView;
   GenEvtView.set().setUserRecord<string>("Type", "Gen");
   RecEvtView.set().setUserRecord<string>("Type", "Rec");
   // set physics process   
   GenEvtView.set().setUserRecord<string>("Process", fProcess);
   RecEvtView.set().setUserRecord<string>("Process", fProcess);

   // Generator stuff
   analyzeGenInfo(iEvent, GenEvtView);
   analyzeGenJets(iEvent, GenEvtView);
   if( fGenOnly == false ){ //only if info is in event
     analyzeGenMET(iEvent, GenEvtView);
     // Reconstructed stuff
     analyzeRecVertices(iEvent, RecEvtView);
     analyzeRecMuons(iEvent, RecEvtView);
     analyzeRecElectrons(iEvent, RecEvtView);
     analyzeRecJets(iEvent, RecEvtView);
     analyzeRecMET(iEvent, RecEvtView);
     analyzeRecGammas(iEvent, RecEvtView);
   }

   Matcher->matchObjects(GenEvtView, RecEvtView);
   //matchObjects(GenEvtView, RecEvtView);

   // set event class strings
   GenEvtView.set().setUserRecord<string>("EventClass", getEventClass(GenEvtView));
   RecEvtView.set().setUserRecord<string>("EventClass", getEventClass(RecEvtView));
   
   fePaxFile.storeObject(GenEvtView);
   fePaxFile.storeObject(RecEvtView);
   if (fDebug > 0) {
      cout << "UserRecord  " <<  "   e   " << "  mu   " << " Gamma " << " KtJet " << "  MET  " << endl;
      cout << "Found (MC): " << setw(4) << GenEvtView.get().findUserRecord<int>("NumEle", 0) 
           << setw(7) << GenEvtView.get().findUserRecord<int>("NumMuon", 0)
           << setw(7) << GenEvtView.get().findUserRecord<int>("NumGamma", 0) 
           << setw(4) << GenEvtView.get().findUserRecord<int>("NumKtJet", 0) << "/" 
           << GenEvtView.get().findUserRecord<int>("NumItCone5Jet", 0) << "/" 
           << GenEvtView.get().findUserRecord<int>("NumMidCone5Jet", 0) << "/" 
           << GenEvtView.get().findUserRecord<int>("NumMidCone7Jet", 0)  
           << setw(7) << GenEvtView.get().findUserRecord<int>("NumMET", 0) << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView.get().findUserRecord<int>("NumEle", 0)    
           << setw(7) << RecEvtView.get().findUserRecord<int>("NumMuon", 0) 
           << setw(7) << RecEvtView.get().findUserRecord<int>("NumGamma", 0)
           << setw(4) << RecEvtView.get().findUserRecord<int>("NumKtJet", 0) << "/"
           << RecEvtView.get().findUserRecord<int>("NumItCone5Jet", 0) << "/"
           << RecEvtView.get().findUserRecord<int>("NumMidCone5Jet", 0) << "/"
           << RecEvtView.get().findUserRecord<int>("NumMidCone7Jet", 0) 
           << setw(7) << RecEvtView.get().findUserRecord<int>("NumMET", 0) << endl;
   }

   if (fDebug > 0) { 
      cout << "Gen Event Type: " << GenEvtView.get().findUserRecord<string>("EventClass") << endl;
      cout << "Rec Event Type: " << RecEvtView.get().findUserRecord<string>("EventClass") << endl;
   }   
   fePaxFile.writeEvent(RecEvtView.get().findUserRecord<string>("EventClass"));
}

// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

  //gen particles
  edm::Handle<reco::CandidateCollection> genParticleHandel;
  iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );

   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;

   //save mother of stable particle
   const Candidate* p_mother; 
   int mother = 0;

   //FIXME: get primary vertex also from candidates
   edm::Handle<edm::HepMCProduct> HepMC_Handle ;
   iEvent.getByLabel( fHepMCLabel, HepMC_Handle );
   const  HepMC::GenEvent* myGenEvent = HepMC_Handle->GetEvent();
   
   // Store primary Vertex:
   pxl::VertexRef GenVtx = EvtView.set().create<pxl::Vertex>();
   GenVtx.set().setName("PrimaryVertex");
   GenVtx.set().vector(pxl::set).setX((*(myGenEvent->vertices_begin()))->position().x());
   GenVtx.set().vector(pxl::set).setY((*(myGenEvent->vertices_begin()))->position().y());
   GenVtx.set().vector(pxl::set).setZ((*(myGenEvent->vertices_begin()))->position().z());

   // loop over all particles
   for( reco::CandidateCollection::const_iterator pa = genParticleHandel->begin(); 
	pa != genParticleHandel->end(); ++ pa ) {

     //cast iterator into GenParticleCandidate
     const GenParticleCandidate* p = (const GenParticleCandidate*) &(*pa);

     // fill Gen Muons passing some basic cuts
     if ( abs((p)->pdgId()) == 13 && (p)->status() == 1) {
         if ( MuonMC_cuts(p) ) { 
            pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
            part.set().setName("Muon");
            part.set().setCharge( p->charge() );
            part.set().vector(pxl::set).setPx(p->px());
            part.set().vector(pxl::set).setPy(p->py());
            part.set().vector(pxl::set).setPz(p->pz());
            part.set().vector(pxl::set).setMass(p->mass());
	    part.set().setUserRecord<double>("Vtx_X", p->vx());
	    part.set().setUserRecord<double>("Vtx_Y", p->vy());
	    part.set().setUserRecord<double>("Vtx_Z", p->vz());

	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, p->pt(), p->eta(), p->phi(), 0.2, 0.1);
	    part.set().setUserRecord<double>("GenIso", GenIso);

	    //save mother of stable muon
	    p_mother =p->mother(); 
	    mother = p_mother->pdgId();
	    //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	    while( fabs(mother) == 13){
	      p_mother = p_mother->mother();
	      mother = p_mother->pdgId();
	    }
	    part.set().setUserRecord<int>("mother_id", mother);

	    numMuonMC++; 
	   	   
         }
      }
      // fill Gen Electrons passing some basic cuts
      if ( abs(p->pdgId()) == 11 && p->status() == 1) {
         if ( EleMC_cuts(p) ) { 
            pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
            part.set().setName("Ele");
            part.set().setCharge(p->charge());
            part.set().vector(pxl::set).setPx(p->px());
            part.set().vector(pxl::set).setPy(p->py());
            part.set().vector(pxl::set).setPz(p->pz());
            part.set().vector(pxl::set).setMass(p->mass());
	    part.set().setUserRecord<double>("Vtx_X", p->vx());
	    part.set().setUserRecord<double>("Vtx_Y", p->vy());
	    part.set().setUserRecord<double>("Vtx_Z", p->vz());

	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, p->pt(), p->eta(), p->phi(), 0.2, 0.1);
	    part.set().setUserRecord<double>("GenIso", GenIso);
	    
	    //save mother of stable electron
	    p_mother =p->mother(); 
	    mother = p_mother->pdgId();
	    //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	    while( fabs(mother) == 13){
	      p_mother = p_mother->mother();
	      mother = p_mother->pdgId();
	    }
	    part.set().setUserRecord<int>("mother_id", mother);

	    numEleMC++; 

         }
      }
      // fill Gen Gammas passing some basic cuts
      if ( abs(p->pdgId()) == 22 && p->status() == 1) {
         if ( GammaMC_cuts(p) ) { 
            pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
            part.set().setName("Gamma");
            part.set().setCharge(0);
            part.set().vector(pxl::set).setPx(p->px());
            part.set().vector(pxl::set).setPy(p->py());
            part.set().vector(pxl::set).setPz(p->pz());
            part.set().vector(pxl::set).setMass(p->mass());

	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, 0., p->eta(), p->phi(), 0.2, 0.1); //0. since gamma not charged!
	    part.set().setUserRecord<double>("GenIso", GenIso);
	    
	    //save mother of stable gamma
	    p_mother =p->mother(); 
	    mother = p_mother->pdgId();
	    //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	    while( fabs(mother) == 13){
	      p_mother = p_mother->mother();
	      mother = p_mother->pdgId();
	    }
	    part.set().setUserRecord<int>("mother_id", mother);

	    numGammaMC++;

         }
      }
   } //end of loop over generated-particles
   if (fDebug > 1)  cout << "MC Found:  " << numMuonMC <<  " mu  " << numEleMC << " e  " << numGammaMC << " gamma" << endl;
   EvtView.set().setUserRecord<int>("NumMuon", numMuonMC);
   EvtView.set().setUserRecord<int>("NumEle", numEleMC);
   EvtView.set().setUserRecord<int>("NumGamma", numGammaMC);
}

// ------------ reading the Generator Jets ------------

void ePaxAnalyzer::analyzeGenJets(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

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
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("KtJet");
         part.set().vector(pxl::set).setPx(genJet->px());
         part.set().vector(pxl::set).setPy(genJet->py());
         part.set().vector(pxl::set).setPz(genJet->pz());
         part.set().vector(pxl::set).setE((genJet->hadEnergy() + genJet->emEnergy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
	 numKtJetMC++;
      }
   }   
   EvtView.set().setUserRecord<int>("NumKtJet", numKtJetMC);

   //Loop over ItCone5GenJets
   for( reco::GenJetCollection::const_iterator genJet = ItCone5_GenJets->begin();
         genJet != ItCone5_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("ItCone5Jet");
         part.set().vector(pxl::set).setPx(genJet->px());
         part.set().vector(pxl::set).setPy(genJet->py());
         part.set().vector(pxl::set).setPz(genJet->pz());
         part.set().vector(pxl::set).setE((genJet->hadEnergy() + genJet->emEnergy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
         numItCone5JetMC++;
      }
   }
   EvtView.set().setUserRecord<int>("NumItCone5Jet", numItCone5JetMC);

   //Loop over MidCone5GenJets
   for( reco::GenJetCollection::const_iterator genJet = MidCone5_GenJets->begin();
         genJet != MidCone5_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("MidCone5Jet");
         part.set().vector(pxl::set).setPx(genJet->px());
         part.set().vector(pxl::set).setPy(genJet->py());
         part.set().vector(pxl::set).setPz(genJet->pz());
         part.set().vector(pxl::set).setE((genJet->hadEnergy() + genJet->emEnergy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
         numMidCone7JetMC++;
      }
   }
   EvtView.set().setUserRecord<int>("NumMidCone5Jet", numMidCone5JetMC);

   //Loop over MidCone7GenJets
   for( reco::GenJetCollection::const_iterator genJet = MidCone7_GenJets->begin();
         genJet != MidCone7_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("MidCone7Jet");
         part.set().vector(pxl::set).setPx(genJet->px());
         part.set().vector(pxl::set).setPy(genJet->py());
         part.set().vector(pxl::set).setPz(genJet->pz());
         part.set().vector(pxl::set).setE((genJet->hadEnergy() + genJet->emEnergy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
         numMidCone7JetMC++;
      }
   }
   EvtView.set().setUserRecord<int>("NumMidCone7Jet", numMidCone5JetMC);
  
   if (fDebug > 1) cout << "Found MC Jets:  " << numKtJetMC << " Kt  " << numItCone5JetMC << " It5  " 
        << numMidCone5JetMC << " Mid5 " << numMidCone7JetMC << " Mid7 " << endl;
   
}

// ------------ reading the Generator MET ------------

void ePaxAnalyzer::analyzeGenMET(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   edm::Handle<reco::GenMETCollection> GenMet;
   iEvent.getByLabel(fMETMCLabel, GenMet);
   const GenMETCollection *genmetcol = GenMet.product();
   const GenMET genmet = genmetcol->front();  // MET exists only once!
 
   int numMETMC = 0; //means no MET in event

   pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
   part.set().setName("MET");
   part.set().vector(pxl::set).setPx(genmet.px());
   part.set().vector(pxl::set).setPy(genmet.py());
   part.set().vector(pxl::set).setPz(0.);
   part.set().vector(pxl::set).setMass(0.);
   part.set().setUserRecord<double>("sumEt", genmet.sumEt());
   part.set().setUserRecord<double>("mEtSig", genmet.mEtSig());
   //fill additional jet-related infos
   part.set().setUserRecord<double>("EmE", genmet.emEnergy());
   part.set().setUserRecord<double>("HadE", genmet.hadEnergy());
   part.set().setUserRecord<double>("InvE", genmet.invisibleEnergy());

   if (fDebug > 1) cout << "GenMET before muon corr: Px = " << genmet.px() << "   Py = " << genmet.py() << "   Pt = " << part.get().vector().getPt() << endl;
   // Perform Muon Corrections!
   // loop over muons and subtract them
   // FIXME: Really only correct for selected muons? Better take GenMET-collection without muons?
   if (EvtView.get().findUserRecord<int>("NumMuon") > 0) { 
      for (pxl::Objects::TypeIterator<pxl::Particle> iter(EvtView().getObjects()); !iter.isDone(); iter.next()) { 
         if (iter.object().get().getName() == "Muon")  {
            if (fDebug > 1) cout << "Correcting with " << iter.object().get().getName() << " px = " << iter.object().get().vector().getPx() 
                                 << " Py = " << iter.object().get().vector().getPy() << endl;
            part.set() -= iter.object().get();
         }
      }
   } 
   if (fDebug > 1) cout << "GenMET after muon corr: Px = " << part.get().vector().getPx() << "   Py = " << part.get().vector().getPy() << "   Pt = " << part.get().vector().getPt() << endl;  

   //reset eta-info after muon corrections
   part.set().vector(pxl::set).setPz(0.);  
   part.set().vector(pxl::set).setMass(0.);   
   //there is always MET in event, just decide if cuts passed (do this after muon corrections!)
   if (METMC_cuts(part)) { 
     numMETMC++; 
   }
   EvtView.set().setUserRecord<int>("NumMET", numMETMC);
   if (numMETMC && fDebug > 1) cout << "Event contains MET" << endl;
}

// ------------ reading Reconstructed Primary Vertices ------------

void ePaxAnalyzer::analyzeRecVertices(const edm::Event& iEvent, pxl::EventViewRef EvtView) {
  
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(fVertexRecoLabel, vertices);
   
   int numVertices = 0;
   
   for(reco::VertexCollection::const_iterator  vertex = vertices->begin(); vertex != vertices->end(); ++vertex ) {
      //only fill primary vertex if cuts passed
      if (Vertex_cuts(vertex)) { 
         pxl::VertexRef vtx = EvtView.set().create<pxl::Vertex>();
         vtx.set().setName("PrimaryVertex");
         vtx.set().vector(pxl::set).setX(vertex->x());
         vtx.set().vector(pxl::set).setY(vertex->y());
         vtx.set().vector(pxl::set).setZ(vertex->z());
	 // chi2 of vertex-fit
         vtx.set().setUserRecord<double>("NormChi2", vertex->normalizedChi2() );        
         // number of tracks with origin in that vertex???
	 vtx.set().setUserRecord<int>("NumTracks", vertex->tracksSize());
         numVertices++;
      }
   }
   EvtView.set().setUserRecord<int>("NumVertices", numVertices); 
}


// ------------ reading Reconstructed Muons ------------

void ePaxAnalyzer::analyzeRecMuons(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   edm::Handle<reco::MuonCollection> muons;
//   edm::Handle<reco::MuonCollection> SAmuons;
   iEvent.getByLabel(fMuonRecoLabel, muons);
//   iEvent.getByLabel(fSAMuonRecoLabel, SAmuons);
 
   int numMuonRec = 0;
//   int numSAMuonRec = 0;

   for(reco::MuonCollection::const_iterator  muon = muons->begin(); 
         muon != muons->end(); ++muon ) {
      if (Muon_cuts(muon)) { 
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("Muon");
         part.set().setCharge(muon->charge());
         part.set().vector(pxl::set).setPx(muon->px());
         part.set().vector(pxl::set).setPy(muon->py());
         part.set().vector(pxl::set).setPz(muon->pz());
         part.set().vector(pxl::set).setE(muon->energy());
         part.set().setUserRecord<double>("Vtx_X", muon->vx());
         part.set().setUserRecord<double>("Vtx_Y", muon->vy());
         part.set().setUserRecord<double>("Vtx_Z", muon->vz()); 
/*	 if (fDebug > 1) {
	    const HitPattern& Pattern = muon->combinedMuon()->hitPattern();
	    cout << "Pattern gives: Valid/Lost: " << endl
	    	 << "	     Total: " << Pattern.numberOfValidHits() << " / " << Pattern.numberOfLostHits() << endl
	    	 << "	      Muon: " << Pattern.numberOfValidMuonHits() << " / " << Pattern.numberOfLostMuonHits() << endl
	    	 << "	   Tracker: " << Pattern.numberOfValidTrackerHits() << " / " << Pattern.numberOfLostTrackerHits() << endl
	    	 << "	     Pixel: " << Pattern.numberOfValidPixelHits() << " / " << Pattern.numberOfLostPixelHits() << endl;
	    
		 cout << "Muon gives: Valid/Lost: " << muon->combinedMuon()->numberOfValidHits() << " / " << muon->combinedMuon()->numberOfLostHits() << endl;
         } */
         part.set().setUserRecord<double>("NormChi2", muon->combinedMuon()->normalizedChi2());
         part.set().setUserRecord<int>("ValidHits", muon->combinedMuon()->numberOfValidHits());
         part.set().setUserRecord<int>("LostHits", muon->combinedMuon()->numberOfLostHits()); 
	 
	 
	 // TEMPORARY: calculate isolation ourselves
	 double CaloIso = IsoCalSum(iEvent, 0., muon->track()->outerEta(), muon->track()->outerPhi(), 0.2, 0.1);
	 double TrkIso = IsoTrkSum(iEvent,  muon->track()->pt(), muon->track()->eta(), muon->track()->phi(), 0.2, 0.1);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);

         numMuonRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumMuon", numMuonRec);
  
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 
  
/*   for(reco::MuonCollection::const_iterator  muon = SAmuons->begin();
         muon != SAmuons->end(); ++muon ) {
      if (Muon_cuts(muon)) {
         cout << " Found a Rec Muon: \n"
              << "    pt : " << muon->pt() << endl
              << "    eta: " << muon->eta() << endl
              << "    q  : " << muon->charge() << endl;
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("SAMuon");
         part.set().setCharge(muon->charge());
         part.set().vector(pxl::set).setPx(muon->px());
         part.set().vector(pxl::set).setPy(muon->py());
         part.set().vector(pxl::set).setPz(muon->pz());
         part.set().vector(pxl::set).setE(muon->energy());
         part.set().setUserRecord<double>("Vtx_X", muon->vx());
         part.set().setUserRecord<double>("Vtx_Y", muon->vy());
         part.set().setUserRecord<double>("Vtx_Z", muon->vz());
         // get isolation ;0 
         // part.set().setUserRecord<double>("isolation", 0.9);
         numSAMuonRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumSAMuon", numSAMuonRec);
*/
}

// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

//   edm::Handle<SiStripElectronCollection> electrons;
   edm::Handle<PixelMatchGsfElectronCollection> pixelelectrons;
//   iEvent.getByLabel(fElectronRecoLabel, electrons);
   iEvent.getByLabel(fPixelMatchElectronRecoLabel, pixelelectrons);
   // Get association maps linking BasicClusters to ClusterShape
   edm::Handle<reco::BasicClusterShapeAssociationCollection> barrelClShpHandle;
   iEvent.getByLabel(fBarrelClusterShapeAssocProducer, barrelClShpHandle);
   edm::Handle<reco::BasicClusterShapeAssociationCollection> endcapClShpHandle;
   iEvent.getByLabel(fEndcapClusterShapeAssocProducer, endcapClShpHandle);

//   int numEleRec = 0;
   int numPixelEleRec = 0;   
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;

/*   for ( SiStripElectronCollection::const_iterator ele = electrons->begin(); 
            ele != electrons->end(); ++ele ) {
      if (Ele_cuts(ele)) { 
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("SiEle");
         part.set().setCharge(ele->charge());
         part.set().vector(pxl::set).setPx(ele->px());
         part.set().vector(pxl::set).setPy(ele->py());
         part.set().vector(pxl::set).setPz(ele->pz());
         part.set().vector(pxl::set).setE(ele->energy());         
         part.set().setUserRecord<double>("Vtx_X", ele->vx());
         part.set().setUserRecord<double>("Vtx_Y", ele->vy());
         part.set().setUserRecord<double>("Vtx_Z", ele->vz());
	 
	 part.set().setUserRecord<double>("SClusterE", ele->superCluster()->energy());
	
	 //DOES NOT WORK IN CMSSW_1_2_0, try with later version //fEleRec_TrackerP->push_back(ele->track()->p());
	 //DOES NOT WORK IN CMSSW_1_2_0, try with later version //fEleRec_TrackerNormChi2->push_back(ele->track()->normalizedChi2());
	 //DOES NOT WORK IN CMSSW_1_2_0, try with later version //fEleRec_TrackerNhitsValid->push_back(ele->track()->numberOfValidHits());
	 //hopefully soon more info on deposited HCAL-energy and shower shapes, as well as isolation (->PixelMatch(Gsf)Electron) !!!
           
	 numEleRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumEle", numEleRec);

   if (fDebug > 1) cout << "RecEle:  " << numEleRec << endl;
*/

   for ( PixelMatchGsfElectronCollection::const_iterator ele = pixelelectrons->begin();
            ele != pixelelectrons->end(); ++ele ) {
      if (Ele_cuts(ele)) {
         if (fDebug > 1) {
	    cout << "Electron Energy scale corrected: " << ele->isEnergyScaleCorrected() 
	         << "  Momentum corrected: " << ele->isMomentumCorrected() << endl;
         }
	 pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("Ele");
         part.set().setCharge(ele->charge());
         part.set().vector(pxl::set).setPx(ele->px());
         part.set().vector(pxl::set).setPy(ele->py());
         part.set().vector(pxl::set).setPz(ele->pz());
         part.set().vector(pxl::set).setE(ele->energy());
         part.set().setUserRecord<double>("Vtx_X", ele->vx());
         part.set().setUserRecord<double>("Vtx_Y", ele->vy());
         part.set().setUserRecord<double>("Vtx_Z", ele->vz());
	 part.set().setUserRecord<float>("EoP", ele->eSuperClusterOverP());
         part.set().setUserRecord<float>("HoE", ele->hadronicOverEm());
	 part.set().setUserRecord<float>("DEtaSCVtx", ele->deltaEtaSuperClusterTrackAtVtx());
	 part.set().setUserRecord<float>("DPhiSCVtx", ele->deltaPhiSuperClusterTrackAtVtx());
	 //! the seed cluster eta - track eta at calo from outermost state
	 part.set().setUserRecord<float>("DEtaSeedTrk", ele->deltaEtaSeedClusterTrackAtCalo());
	 part.set().setUserRecord<float>("DPhiSeedTrk", ele->deltaPhiSeedClusterTrackAtCalo());
	 //part.set().setUserRecord<double>("SCE", ele->superCluster()->energy());
         // the super cluster energy corrected by EnergyScaleFactor
	 part.set().setUserRecord<float>("SCE", ele->caloEnergy());
         // the errors on the supercluster energy and track momentum
         part.set().setUserRecord<float>("SCEErr", ele->caloEnergyError());
	 // the errors on the track momentum
         part.set().setUserRecord<float>("PErr", ele->trackMomentumError());	 
         part.set().setUserRecord<double>("TrackerP", ele->gsfTrack()->p());
	 //the seed cluster energy / track momentum at calo from outermost state
	 part.set().setUserRecord<float>("ESCSeedPout", ele->eSeedClusterOverPout());
         //part.set().setUserRecord<double>("NormChi2", ele->gsfTrack()->normalizedChi2());
         part.set().setUserRecord<int>("ValidHits", ele->gsfTrack()->numberOfValidHits());
	 part.set().setUserRecord<int>("Class", ele->classification());
         
	 // Get the supercluster (ref) of the Electron
	 // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
	 // a SuperClusterCollection is a std::vector<SuperCluster>
	 // although we get a vector of SuperClusters an electron is only made out of ONE SC
	 // therefore only the first element of the vector should be available!
	 const SuperClusterRef SCRef = ele->superCluster();
	 //const BasicClusterRef& SCSeed = SCRef.seed(); 
         // Find the entry in the map corresponding to the seed BasicCluster of the SuperCluster
         DetId id = SCRef->seed()->getHitsByDetId()[0];
         if (id.subdetId() == EcalBarrel) {
            seedShpItr = barrelClShpHandle->find(SCRef->seed());
         } else {
            seedShpItr = endcapClShpHandle->find(SCRef->seed());
         }
	 
         const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;
	 part.set().setUserRecord<double>("e3x3", seedShapeRef->e3x3());
	 part.set().setUserRecord<double>("e5x5", seedShapeRef->e5x5());
	 part.set().setUserRecord<double>("EtaEta", seedShapeRef->covEtaEta());
	 part.set().setUserRecord<double>("EtaPhi", seedShapeRef->covEtaPhi());
	 part.set().setUserRecord<double>("PhiPhi", seedShapeRef->covPhiPhi());

	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection)
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
	 /*cout<<"seedShapeRef->eMaxId().rawId(): "<<seedShapeRef->eMaxId().rawId()<<endl;
	 cout<<"SCRef->seed()->phi(): "<<SCRef->seed()->phi()<<endl;
	 cout<<"SCRef->seed()->eta(): "<<SCRef->seed()->eta()<<endl;
	 cout<<"ele->charge(): "<<ele->charge()<<endl;
	 cout<<"ele->pt(): "<<ele->pt()<<endl;
	 cout<<"ele->phi(): "<<ele->phi()<<endl;
	 cout<<"ele->eta(): "<<ele->eta()<<endl;
	 cout<<"ele->energy(): "<<ele->energy()<<endl;
	 cout<<"ele->caloEnergy(): "<<ele->caloEnergy()<<endl;
	 cout<<"ele->superCluster()->energy(): "<<ele->superCluster()->energy()<<endl;
	 cout<<"ele->superCluster()->rawEnergy(): "<<ele->superCluster()->rawEnergy()<<endl;*/
	 //if(ele->eSuperClusterOverP() - ele->caloEnergy()/ele->trackMomentumAtVtx().R() >0.01){
	 
	 /*cout<<""<<endl;
	 cout<<"ele->eSuperClusterOverP()                                                 : "<<ele->eSuperClusterOverP()<<endl;
	 cout<<"EOverP by hand (ele->caloEnergy()/ele->trackMomentumAtVtx().R())          : "<< ele->caloEnergy()/ele->trackMomentumAtVtx().R()<<endl;
	 cout<<"ele->trackMomentumAtVtx().R(): "<<ele->trackMomentumAtVtx().R()<<endl;
	 cout<<"ele->TrackPostionAtVtx(): "<<ele->TrackPositionAtVtx()<<endl;
	 cout<<"ele->gsfTrack()->p(): "<<ele->gsfTrack()->p()<<endl;
	 cout<<"ele->gsfTrack()->innerMomentum().R(): "<<ele->gsfTrack()->innerMomentum().R()<<endl;
	 cout<<"ele->gsfTrack()->innerPosition(): "<<ele->gsfTrack()->innerPosition()<<endl;
	 cout<<""<<endl;*/
	 
	 //}
	 /*cout<<"ele->gsfTrack()->innerMomentum().R(): "<<ele->gsfTrack()->innerMomentum().R()<<endl;
	 cout<<"ele->trackMomentumAtVtx().R(): "<<ele->trackMomentumAtVtx().R()<<endl;
	 cout<<"ele->gsfTrack()->p(): "<<ele->gsfTrack()->p()<<endl;*/

	 // TEMPORARY: calculate isolation ourselves
	 double CaloPt = ( ele->superCluster()->rawEnergy() + ele->superCluster()->rawEnergy()*ele->hadronicOverEm() ) / cosh(ele->trackMomentumAtCalo().eta());
	 double CaloIso = IsoCalSum(iEvent, CaloPt, ele->trackMomentumAtCalo().eta(), ele->trackMomentumAtCalo().phi(), 0.2, 0.1);
	 double TrkIso = IsoTrkSum(iEvent, ele->gsfTrack()->pt(), ele->gsfTrack()->eta(), ele->gsfTrack()->phi(), 0.2, 0.1);
	 //double TrkIso = IsoTrkSum(iEvent, ele->trackMomentumAtVtx().Rho(), ele->trackMomentumAtVtx().Eta(), ele->trackMomentumAtVtx().Phi(), 0.2, 0.1);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);
	 	 
         numPixelEleRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumEle", numPixelEleRec);
}

// ------------ reading Reconstructed Jets ------------
//
// Which kind of Jets?
//

void ePaxAnalyzer::analyzeRecJets(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

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
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("KtJet");
         part.set().vector(pxl::set).setPx(jet->px());
         part.set().vector(pxl::set).setPy(jet->py());
         part.set().vector(pxl::set).setPz(jet->pz());
         part.set().vector(pxl::set).setE(jet->energy());         
	 part.set().setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part.set().setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part.set().setUserRecord<double>("N90", jet->n90());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
	 numKtJetRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumKtJet", numKtJetRec);

   for(reco::CaloJetCollection::const_iterator jet = ItCone5jets->begin();
           jet != ItCone5jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("ItCone5Jet");
         part.set().vector(pxl::set).setPx(jet->px());
         part.set().vector(pxl::set).setPy(jet->py());
         part.set().vector(pxl::set).setPz(jet->pz());
         part.set().vector(pxl::set).setE(jet->energy());
	 part.set().setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part.set().setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part.set().setUserRecord<double>("N90", jet->n90());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
         numItCone5JetRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumItCone5Jet", numItCone5JetRec);

   for(reco::CaloJetCollection::const_iterator jet = MidCone5jets->begin();
           jet != MidCone5jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("MidCone5Jet");
         part.set().vector(pxl::set).setPx(jet->px());
         part.set().vector(pxl::set).setPy(jet->py());
         part.set().vector(pxl::set).setPz(jet->pz());
         part.set().vector(pxl::set).setE(jet->energy());
	 part.set().setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part.set().setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part.set().setUserRecord<double>("N90", jet->n90());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
         numMidCone5JetRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumMidCone5Jet", numMidCone5JetRec);

   for(reco::CaloJetCollection::const_iterator jet = MidCone7jets->begin();
           jet != MidCone7jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("MidCone7Jet");
         part.set().vector(pxl::set).setPx(jet->px());
         part.set().vector(pxl::set).setPy(jet->py());
         part.set().vector(pxl::set).setPz(jet->pz());
         part.set().vector(pxl::set).setE(jet->energy());
	 part.set().setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part.set().setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part.set().setUserRecord<double>("N90", jet->n90());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
         numMidCone7JetRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumMidCone7Jet", numMidCone7JetRec);

   if (fDebug > 1) cout << "Found Rec Jets:  " << numKtJetRec << " Kt  " << numItCone5JetRec << " It5  "
                        << numMidCone5JetRec << " Mid5 " << numMidCone7JetRec << " Mid7 " << endl;

}

// ------------ reading Reconstructed MET ------------

void ePaxAnalyzer::analyzeRecMET(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   edm::Handle<reco::CaloMETCollection> CaloMet;
   iEvent.getByLabel(fMETRecoLabel, CaloMet);
   const CaloMETCollection *calometcol = CaloMet.product();
   const CaloMET calomet = calometcol->front();  // MET exists only once!
 
   int numMETRec = 0;
   pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
   part.set().setName("MET");
   part.set().vector(pxl::set).setPx(calomet.px());
   part.set().vector(pxl::set).setPy(calomet.py());
   part.set().vector(pxl::set).setPz(0.);
   part.set().vector(pxl::set).setMass(0.);
   part.set().setUserRecord<double>("sumEt", calomet.sumEt());
   part.set().setUserRecord<double>("mEtSig", calomet.mEtSig());
   part.set().setUserRecord<double>("EmEt", calomet.emEtFraction());          //not muon corrected
   part.set().setUserRecord<double>("HadEt", calomet.etFractionHadronic());   //not muon corrected
   part.set().setUserRecord<double>("MaxEtEm", calomet.maxEtInEmTowers());   //not muon corrected
   part.set().setUserRecord<double>("MaxEtHad", calomet.maxEtInHadTowers()); //not muon corrected
   
   
   if (fDebug > 1) {
      cout << "RecMET before muon corr: Px = " << calomet.px() << "   Py = " << calomet.py() << "   Pt = " << part.get().vector().getPt() << endl;   
      cout << " MET (uncorr): ( " << part.get().vector().getPx() << ", " << part.get().vector().getPy() << ", "
           << part.get().vector().getPz() <<  " )    E: " << part.get().vector().getE() << "   Et = "
           << part.get().vector().getEt() << " Theta: " << part.get().vector().getTheta()
           << " Mass: " << part.get().vector().getMass() << endl;
   }
   // Perform Muon Corrections!   
   // loop over muons and subtract them   
   // FIXME: Really only correct for selected muons? Is there official muon correction out yet?
   if (EvtView.get().findUserRecord<int>("NumMuon") > 0) {      
   for (pxl::Objects::TypeIterator<pxl::Particle> iter(EvtView().getObjects()); !iter.isDone(); iter.next()) {
      if (iter.object().get().getName() == "Muon")  {
         if (fDebug > 1) cout << "Correcting with " << iter.object().get().getName() << " px = " << iter.object().get().vector().getPx() 
                              << " Py = " << iter.object().get().vector().getPy() << endl; 
            part.set() -= iter.object().get(); 
         }  
      }   
   } 
   //reset eta-info after muon corrections
   part.set().vector(pxl::set).setPz(0.);  
   part.set().vector(pxl::set).setMass(0.);
   if (fDebug > 1) {
      cout << " MET (corr): ( " << part.get().vector().getPx() << ", " << part.get().vector().getPy() << ", "
           << part.get().vector().getPz() <<  " )    E: " << part.get().vector().getE() << "   Et = "
           << part.get().vector().getEt() << " Theta: " << part.get().vector().getTheta()
           << " Mass: " << part.get().vector().getMass() << endl;
      cout << "RecMET after muon corr: Px = " << part.get().vector().getPx() << "   Py = " << part.get().vector().getPy() 
           << "   Pt = " << part.get().vector().getPt() << endl;   
   } 
   //there is always MET in event, just decide if cuts passed (do this after muon corrections!)
   if (MET_cuts(part)) { 
     numMETRec++;
   }

   EvtView.set().setUserRecord<int>("NumMET", numMETRec);
   if (numMETRec && fDebug > 1) cout << "Found RecMET" << endl;
}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventViewRef EvtView) {
   
   // get Photon Collection
   edm::Handle<reco::PhotonCollection> Photons;
   iEvent.getByLabel(fGammaRecoLabel, Photons);
   // get ECAL Cluster shapes
   edm::Handle<reco::BasicClusterShapeAssociationCollection> barrelClShpHandle;
   iEvent.getByLabel(fBarrelClusterShapeAssocProducer, barrelClShpHandle);
   edm::Handle<reco::BasicClusterShapeAssociationCollection> endcapClShpHandle;
   iEvent.getByLabel(fEndcapClusterShapeAssocProducer, endcapClShpHandle);
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
   // get Handle to HCAL RecHits for HadOverEm calculation:
   edm::Handle<HBHERecHitCollection> hbhe;
   //HBHERecHitMetaCollection *HBHE_RecHits = 0;
   iEvent.getByLabel(fHBHELabel, fHBHEInstanceName, hbhe);  
   HBHERecHitMetaCollection HBHE_RecHits(*hbhe); 
   
   int numGammaRec = 0;
   double HoE = 0.;

   for (reco::PhotonCollection::const_iterator photon = Photons->begin(); 
	photon != Photons->end(); ++photon) {  
      if ( Gamma_cuts(photon) ) { 
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("Gamma");
         part.set().setCharge(0);
         part.set().vector(pxl::set).setPx(photon->px());
         part.set().vector(pxl::set).setPy(photon->py());
         part.set().vector(pxl::set).setPz(photon->pz());
         part.set().vector(pxl::set).setE(photon->energy());
	 // ratio of E(3x3)/ESC
	 part.set().setUserRecord<double>("r9", photon->r9());
	 // ratio of Emax/E(3x3)
	 part.set().setUserRecord<double>("r19", photon->r19());
	 // 5x5 energy
	 part.set().setUserRecord<double>("e5x5", photon->e5x5());
	 /// Whether or not the SuperCluster has a matched pixel seed
	 part.set().setUserRecord<bool>("HasSeed", photon->hasPixelSeed());
	 
	 const SuperClusterRef SCRef = photon->superCluster();
	 //const BasicClusterRef& SCSeed = SCRef.seed(); 
         // Find the entry in the map corresponding to the seed BasicCluster of the SuperCluster
         DetId id = SCRef->seed()->getHitsByDetId()[0];
         if (id.subdetId() == EcalBarrel) {
            seedShpItr = barrelClShpHandle->find(SCRef->seed());
         } else {
            seedShpItr = endcapClShpHandle->find(SCRef->seed());
         }
	 
         const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;
	 part.set().setUserRecord<double>("e3x3", seedShapeRef->e3x3());
	 //part.set().setUserRecord<double>("e5x5", seedShapeRef->e5x5());
	 //cout << "Comparing E5x5 from Photon: " << photon->e5x5() 
	 //     << " and from ClusterShape: " << seedShapeRef->e5x5() << endl; // they are identical! 
  	 part.set().setUserRecord<double>("EtaEta", seedShapeRef->covEtaEta());
	 part.set().setUserRecord<double>("EtaPhi", seedShapeRef->covEtaPhi());
	 part.set().setUserRecord<double>("PhiPhi", seedShapeRef->covPhiPhi());

	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection)
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
         
	 // calculate HadoverEm - now its getting tough ...
	 CaloConeSelector sel(0.1, theCaloGeom.product(), DetId::Hcal);
         GlobalPoint pclu(SCRef->x(),SCRef->y(),SCRef->z());
	 double hcalEnergy = 0.;
         std::auto_ptr<CaloRecHitMetaCollectionV> chosen=sel.select(pclu,HBHE_RecHits);
         for (CaloRecHitMetaCollectionV::const_iterator i=chosen->begin(); i!=chosen->end(); i++) {
            //std::cout << HcalDetId(i->detid()) << " : " << (*i) << std::endl;
            hcalEnergy += i->energy();
         }
         HoE = hcalEnergy/photon->energy();
         //cout << "H/E : " << HoE << endl;
         part.set().setUserRecord<float>("HoE", HoE);

	 // TEMPORARY: calculate isolation ourselves
	 double CaloPt = ( hcalEnergy + photon->superCluster()->rawEnergy() ) / cosh(photon->eta());
	 double CaloIso = IsoCalSum(iEvent, CaloPt, photon->eta(), photon->phi(), 0.2, 0.1);
	 double TrkIso = IsoTrkSum(iEvent, 0., photon->eta(), photon->phi(), 0.2, 0.1);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);

	 numGammaRec++;
      }	 
   }
   EvtView.set().setUserRecord<int>("NumGamma", numGammaRec);
   if (fDebug > 1) cout << "Rec Gamma: " << numGammaRec << endl;
}

// ------------ method returning the EventClassType  ------------

std::string ePaxAnalyzer::getEventClass(pxl::EventViewRef EvtView) {

   ostringstream EventType;
   //set default values to 0 for Gen-only mode
   EventType << EvtView.get().findUserRecord<int>("NumEle", 0) <<  "e"
             << EvtView.get().findUserRecord<int>("NumMuon", 0) << "mu"
             << EvtView.get().findUserRecord<int>("NumGamma", 0) << "gam"
             << EvtView.get().findUserRecord<int>("NumKtJet", 0) << "kt"
             << EvtView.get().findUserRecord<int>("NumMET", 0) << "met";
   EventType.flush();
   return EventType.str();
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

bool ePaxAnalyzer::MuonMC_cuts(const GenParticleCandidate* MCmuon) const {
   //
   if (MCmuon->pt() < 10) return false;
   if (fabs(MCmuon->eta()) > 2.4) return false;
   return true;
}

// ------------ method to define MC-Electron-cuts

bool ePaxAnalyzer::EleMC_cuts(const GenParticleCandidate* MCele) const {
   //
   if (MCele->pt() < 10) return false;
   if (fabs(MCele->eta()) > 5.) return false;
   return true;
}

// ------------ method to define MC-Gamma-cuts

bool ePaxAnalyzer::GammaMC_cuts(const GenParticleCandidate* MCgamma) const {
   //
   if (MCgamma->pt() < 30) return false;
   if (fabs(MCgamma->eta()) > 5.0) return false;
   return true;
}

// ------------ method to define MC-Jet-cuts

bool ePaxAnalyzer::JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const {
   // 
   if (MCjet->pt() < 30) return false;
   if (fabs(MCjet->eta()) > 5.) return false;
   return true;
}

// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::METMC_cuts(const pxl::ParticleRef MCmet) const {
   // 
   if (MCmet.get().vector().getPt() < 30) return false;
   return true; 
}

// ------------ method to define RecVertex-cuts

bool ePaxAnalyzer::Vertex_cuts(reco::VertexCollection::const_iterator vertex) const {
  //
  //check compatibility of vertex with beam spot
  double zV = vertex->z();
  double rV = sqrt( vertex->x()*vertex->x() + vertex->y()*vertex->y() );
  if (fabs(zV)>20. || rV>1. ) return false;
  return true;
}

// ------------ method to define MUON-cuts

bool ePaxAnalyzer::Muon_cuts(reco::MuonCollection::const_iterator muon) const {
   //
   if (muon->pt() < 10)  return false;
   if (fabs(muon->eta()) > 2.4) return false;
   return true;
}

// ------------ method to define ELECTRON-cuts

bool ePaxAnalyzer::Ele_cuts(SiStripElectronCollection::const_iterator ele) const {
   //
   if (ele->pt() < 10) return false;
   if (fabs(ele->eta()) > 5) return false;
   return true;
}

// ------------ method to define ELECTRON-cuts

bool ePaxAnalyzer::Ele_cuts(PixelMatchGsfElectronCollection::const_iterator ele) const {
   //
   if (ele->pt() < 10) return false;
   if (fabs(ele->eta()) > 5) return false;
   return true;
}

// ------------ method to define JET-cuts

bool ePaxAnalyzer::Jet_cuts(reco::CaloJetCollection::const_iterator jet) const {
   //
   if (jet->pt() < 30) return false;
   if (fabs(jet->eta()) > 5) return false;
   return true;
}

// ------------ method to define GAMMA-cuts

bool ePaxAnalyzer::Gamma_cuts(reco::PhotonCollection::const_iterator photon) const {
   //
   if (photon->pt() < 30) return false;
   if (fabs(photon->eta()) > 5.0) return false;
   return true;
}

// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::MET_cuts(const pxl::ParticleRef met) const {
   // 
   if (met.get().vector().getPt() < 30) return false;
   return true;
}


// TEMPORARY STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//------------------------------------------------------------------------------

double ePaxAnalyzer::IsoCalSum (const edm::Event& iEvent, double ParticleCalPt, double ParticleCalEta, double ParticleCalPhi, double iso_DR, double iso_Seed){
// Computes the sum of Pt inside a cone of R=iso_DR
// using 4-vectors stored in CaloTower objects

  float sum = 0.;

  edm::Handle<CaloTowerCollection> CaloTowerData ;
  iEvent.getByLabel( "towerMaker", CaloTowerData );

 for( CaloTowerCollection::const_iterator tower = CaloTowerData->begin(); 
         tower != CaloTowerData->end(); ++tower ) {
   if (tower->energy() > iso_Seed){
      float eta = tower->eta();
      float phi = tower->phi();
      float DR = GetDeltaR(ParticleCalEta, eta, ParticleCalPhi, phi);
      if (DR <= 0.) {DR = 0.001;}
      if (DR < iso_DR){
          float pt = tower->energy() / cosh(eta);
          sum += pt;
      }
    }
  }

  sum -= ParticleCalPt;
  
  return sum;

}

//------------------------------------------------------------------------------

double ePaxAnalyzer::IsoTrkSum (const edm::Event& iEvent, double ParticleTrkPt, double ParticleTrkEta, double ParticleTrkPhi, double iso_DR, double iso_Seed){
// Computes the sum of Pt inside a cone of R=iso_DR
// using 4-vectors stored in Track objects

  float sum = 0.;

  edm::Handle<TrackCollection> TrackData ;
  iEvent.getByLabel( "ctfWithMaterialTracks", TrackData );

  for( reco::TrackCollection::const_iterator track = TrackData->begin(); 
         track != TrackData->end(); ++track ) {
    if (track->p() > iso_Seed){
      float eta = track->eta();
      float phi = track->phi();
      float DR = GetDeltaR(ParticleTrkEta, eta, ParticleTrkPhi, phi);
      if (DR <= 0.) {DR = 0.001;}
      if (DR < iso_DR){
          sum += track->pt();
      }
    }
  }
  
  sum -= ParticleTrkPt;

  return sum;

}

//------------------------------------------------------------------------------

double ePaxAnalyzer::IsoGenSum (const edm::Event& iEvent, double ParticleGenPt, double ParticleGenEta, double ParticleGenPhi, double iso_DR, double iso_Seed){
// Computes the sum of Pt inside a cone of R=iso_DR
// using 4-vectors stored in GenParticle objects

  float sum = 0.;

  //gen particles
  edm::Handle<reco::CandidateCollection> genParticleHandel;
  iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );

  // loop over all particles
  for( reco::CandidateCollection::const_iterator pa = genParticleHandel->begin(); 
       pa != genParticleHandel->end(); ++ pa ) {

    //cast iterator into GenParticleCandidate
    const GenParticleCandidate* p = (const GenParticleCandidate*) &(*pa);

    // only consider stable partciles and charged particles in order to be more comparable with track-isolation
    if ( p->status() == 1 && p->charge() != 0 ) {

      if (p->energy() > iso_Seed){
	float eta = p->eta();
	float phi = p->phi();
	float DR = GetDeltaR(ParticleGenEta, eta, ParticleGenPhi, phi);
	if (DR <= 0.) {DR = 0.001;}
	if (DR < iso_DR){
          sum += p->pt();
	}
      }
    }
  }
  
  sum -= ParticleGenPt;

  return sum;

}

//------------------------------------------------------------------------------

double ePaxAnalyzer::DeltaPhi(double v1, double v2)
{ // Computes the correctly normalized phi difference
  // v1, v2 = phi of object 1 and 2
 
 double pi    = 3.141592654;
 double twopi = 6.283185307;
 
 double diff = fabs(v2 - v1);
 double corr = twopi - diff;
 if (diff < pi){ return diff;} else { return corr;} 
 
}

double ePaxAnalyzer::GetDeltaR(double eta1, double eta2, double phi1, double phi2)
{ // Computes the DeltaR of two objects from their eta and phi values

 return sqrt( (eta1-eta2)*(eta1-eta2) 
            + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
}



//define this as a plug-in
DEFINE_FWK_MODULE(ePaxAnalyzer);
