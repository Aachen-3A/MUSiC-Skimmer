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
#include <iostream>
// include Message Logger for Debug
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h" 
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//ok

// necessary objects:
#include "FWCore/Framework/interface/ESHandle.h"
//ok

// for electron shapes:
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
//ok

// for ECAL enumerator
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
//ok

// for HCAL navigation used in HadOverEm calculation
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
//ok

//for GenPartciles
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

//for TrackingVertex
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

//for muon-isolation
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/TrackReco/interface/Track.h"

//for Electron ID
#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"

//Photon pi0 rejection
#include "DataFormats/EgammaCandidates/interface/PhotonPi0DiscriminatorAssociation.h"

//no longer available in CMSSW 2.0.7
//#include "DataFormats/EgammaCandidates/interface/ConvertedPhoton.h"

//for electron-isolation
#include "DataFormats/Candidate/src/classes.h"
//no longer there in CMSSW 2.0.7
//#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoNumCollection.h"
//#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoCollection.h"

//for Trigger Bits
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
//ok

//For L1 and Hlt objects
#include "DataFormats/Common/interface/RefToBase.h"

//no longer there in CMSSW 2.0.7
//#include "DataFormats/HLTReco/interface/HLTFilterObject.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//ok

//math stuff from Physics tools
#include "DataFormats/Math/interface/deltaR.h"

//PAT related stuff
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


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
   //fGenOnly = iConfig.getUntrackedParameter<bool>("GenOnly");
   // Debugging
   fDebug = iConfig.getUntrackedParameter<int>("debug");
   // The labels used in cfg-file 
   //fTruthVertexLabel = iConfig.getUntrackedParameter<string>("TruthVertexLabel");
   //fgenParticleCandidatesLabel  = iConfig.getUntrackedParameter<string>("genParticleCandidatesLabel");
   //fKtJetMCLabel = iConfig.getUntrackedParameter<string>("KtJetMCLabel");
   //fItCone5JetMCLabel = iConfig.getUntrackedParameter<string>("ItCone5JetMCLabel");  
   //fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fVertexRecoLabel = iConfig.getUntrackedParameter<string>("VertexRecoLabel");
   //fSAMuonRecoLabel = iConfig.getUntrackedParameter<string>("SAMuonRecoLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   //fBarrelClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("barrelClusterShapeAssociation");
   //fEndcapClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("endcapClusterShapeAssociation");   
   //fPixelMatchElectronRecoLabel = iConfig.getUntrackedParameter<string>("PixelMatchElectronRecoLabel");
   //fElectronIDAssocProducer = iConfig.getParameter<std::string>("ElectronIDAssocProducer");
   //fElectronHcalIsolationProducer = iConfig.getParameter<std::string>("ElectronHcalIsolationProducer");
   //fElectronEcalIsolationProducer = iConfig.getParameter<std::string>("ElectronEcalIsolationProducer");
   //fElectronTrackIsolationProducer = iConfig.getParameter<std::string>("ElectronTrackIsolationProducer");
   //fElectronTrackNumProducer = iConfig.getParameter<std::string>("ElectronTrackNumProducer");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   //fGammaHcalIsolationProducer = iConfig.getParameter<std::string>("GammaHcalIsolationProducer");
   //fGammaEcalIsolationProducer = iConfig.getParameter<std::string>("GammaEcalIsolationProducer");
   //fGammaTrackIsolationProducer = iConfig.getParameter<std::string>("GammaTrackIsolationProducer");
   //fGammaTrackNumProducer = iConfig.getParameter<std::string>("GammaTrackNumProducer");
   //fKtJetRecoLabel = iConfig.getUntrackedParameter<string>("KtJetRecoLabel");
   //fItCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("ItCone5JetRecoLabel");
   //fL2L3JESic5JetRecoLabel = iConfig.getUntrackedParameter<string>("L2L3JESic5JetRecoLabel");
   //fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");
   //fMETCorrRecoLabel = iConfig.getUntrackedParameter<string>("METCorrRecoLabel");
   //fHBHELabel = iConfig.getUntrackedParameter<string>("fHBHELabel");
   //fHBHEInstanceName = iConfig.getUntrackedParameter<string>("fHBHEInstanceName");
   
   
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
  //cout<<"Event Number: "<<fNumEvt<<endl;
   // set event counter   
   fNumEvt++;    // set this as setUserRecord ?
   
   // get the calorimeter geometry in order to navigate through it for HadOverEm calculation
   //probably no longer needed and not available in this way in CMSSW 2.0.7 anyway
   //iSetup.get<IdealGeometryRecord>().get(theCaloGeom);

   // create two ePaxEventViews for Generator/Reconstructed Objects
   pxl::EventView GenEvtView;
   pxl::EventView RecEvtView;
   GenEvtView.set().setUserRecord<string>("Type", "Gen");
   RecEvtView.set().setUserRecord<string>("Type", "Rec");
   
   // Store Run and Event ID
   GenEvtView.set().setUserRecord<int>("Run", iEvent.id().run());
   GenEvtView.set().setUserRecord<int>("ID", iEvent.id().event());
   RecEvtView.set().setUserRecord<int>("Run", iEvent.id().run());
   RecEvtView.set().setUserRecord<int>("ID", iEvent.id().event());
   //cout << "Run " << iEvent.id().run() << "   EventID = " << iEvent.id().event() << endl;


   double processID = 0.;
   double pthat = 0.;

   
   // else use stuff inside AOD : can it be done in PAT?
  
      //edm::Handle< int > genProcessID;
      //iEvent.getByLabel( "genEventProcID", genProcessID );
      //processID = *genProcessID;
 
      //edm::Handle< double > genEventScale;
      //iEvent.getByLabel( "genEventScale", genEventScale );
      //pthat = *genEventScale;
    

   //set process name
   GenEvtView.set().setUserRecord<string>("Process", fProcess);
   RecEvtView.set().setUserRecord<string>("Process", fProcess);
   
 
   
   // store the ID in both event views? 
   GenEvtView.set().setUserRecord<double>("pthat", pthat);
   RecEvtView.set().setUserRecord<double>("pthat", pthat); 

   
   double weight = 1.;

   //cout << "Weight = " << weight << endl;
   //cout << "Process ID: " << processID << " and Event Scale (pthat): " << pthat << endl;

   GenEvtView.set().setUserRecord<double>("Weight", weight);
   RecEvtView.set().setUserRecord<double>("Weight", weight);
   
   GenEvtView.set().setUserRecord<double>("ProcID", processID);
   RecEvtView.set().setUserRecord<double>("ProcID", processID);

   // Generator stuff
   //analyzeGenInfo(iEvent, GenEvtView);
   //analyzeGenJets(iEvent, GenEvtView);
   //if( fGenOnly == false ){ //only if info is in event
     //analyzeGenMET(iEvent, GenEvtView);

     //Trigger bits
     analyzeTrigger(iEvent, RecEvtView);
   
   // Reconstructed stuff
     analyzeRecVertices(iEvent, RecEvtView);
     analyzeRecMuons(iEvent, RecEvtView);
     analyzeRecElectrons(iEvent, RecEvtView);
     //analyzeRecJets(iEvent, RecEvtView);
     //analyzeRecMET(iEvent, RecEvtView);
     analyzeRecGammas(iEvent, RecEvtView);


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
           << GenEvtView.get().findUserRecord<int>("NumMidCone5Jet", 0)
           << setw(7) << GenEvtView.get().findUserRecord<int>("NumMET", 0) << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView.get().findUserRecord<int>("NumEle", 0)    
           << setw(7) << RecEvtView.get().findUserRecord<int>("NumMuon", 0) 
           << setw(7) << RecEvtView.get().findUserRecord<int>("NumGamma", 0)
           << setw(4) << RecEvtView.get().findUserRecord<int>("NumKtJet", 0) << "/"
           << RecEvtView.get().findUserRecord<int>("NumItCone5Jet", 0) 
           << setw(7) << RecEvtView.get().findUserRecord<int>("NumMET", 0) << endl;
   }

   if (fDebug > 0) { 
      cout << "Gen Event Type: " << GenEvtView.get().findUserRecord<string>("EventClass") << endl;
      cout << "Rec Event Type: " << RecEvtView.get().findUserRecord<string>("EventClass") << endl;
   }   
   fePaxFile.writeEvent(RecEvtView.get().findUserRecord<string>("EventClass"));
}


/*

// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   //gen particles
   edm::Handle<reco::CandidateCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );
   
   //TrackingVertexCollection ONLY available in RECO, so need ugly catch workaround...
   try{
     //get primary vertex from TrackingVertex. First vertex should be equal to HepMC-info, but HepMC is depreciated in AOD...
     edm::Handle<TrackingVertexCollection> TruthVertexContainer;
     iEvent.getByLabel(fTruthVertexLabel, TruthVertexContainer);
     const TrackingVertexCollection* tVC = TruthVertexContainer.product();
     
     // take only first vertex as this should correspond to generated primary vertex. Seems to be carbon-copy of HepMC-genVertex
     // Question: What about pile-up/min-bias? Will there be additional primary vertices??? How do we find them (BX-id, event-id) and should we save them?
     TrackingVertexCollection::const_iterator EventVertices = tVC->begin(); 
     
     // Store primary Vertex:
     pxl::VertexRef GenVtx = EvtView.set().create<pxl::Vertex>();
     GenVtx.set().setName("PV");
     GenVtx.set().vector(pxl::set).setX(EventVertices->position().x());
     GenVtx.set().vector(pxl::set).setY(EventVertices->position().y());
     GenVtx.set().vector(pxl::set).setZ(EventVertices->position().z());
     // we only have a single PV at Generator Level. Due to EventView Consistency .i.e. GenEvtView should look identical to RecEvtView
     // this variable is explicetly set
     EvtView.set().setUserRecord<int>("NumVertices", 1);               
     // do we need this BX/event identification???
     //GenVtx.set().setUserRecord<int>("Vtx_BX", EventVertices->eventId().bunchCrossing());
     //GenVtx.set().setUserRecord<int>("Vtx_event", EventVertices->eventId().event());

   }catch(...){ //for AOD use daughter of first GenParticle, this should be equal to generated primary vertex...  

     const GenParticleCandidate* p = (const GenParticleCandidate*) &(*genParticleHandel->begin()); //this is the incoming proton
     pxl::VertexRef GenVtx = EvtView.set().create<pxl::Vertex>();
     GenVtx.set().setName("PV");
     if(p->daughter(0)){ //check validity, otherwise sometimes segmentation violation - WHY ???
       GenVtx.set().vector(pxl::set).setX(p->daughter(0)->vx()); //need daughter since first particle (proton) has position zero
       GenVtx.set().vector(pxl::set).setY(p->daughter(0)->vy());
       GenVtx.set().vector(pxl::set).setZ(p->daughter(0)->vz());
     }else{ //set dummy values to catch segmentation violation
       GenVtx.set().vector(pxl::set).setX(0.); //need daughter since first particle (proton) has position zero
       GenVtx.set().vector(pxl::set).setY(0.);
       GenVtx.set().vector(pxl::set).setZ(0.);
     }
     EvtView.set().setUserRecord<int>("NumVertices", 1);       
   }

  
   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;

   //save mother of stable particle
   const Candidate* p_mother; 
   int mother = 0;
//   int numOfProton = 0;
//   vector<const Candidate*> cand_prot1;
//   vector<const Candidate*> cand_prot2;
   
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
	    double GenIso = IsoGenSum(iEvent, p->pt(), p->eta(), p->phi(), 0.3, 1.5);
	    part.set().setUserRecord<double>("GenIso", GenIso);

	    //save mother of stable muon
	    p_mother =p->mother(); 
	    if (p_mother != 0) {
	       mother = p_mother->pdgId();
	       //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	       while( abs(mother) == 13 ){
	         p_mother = p_mother->mother();
	         mother = p_mother->pdgId();
	       }	       
	       part.set().setUserRecord<int>("mother_id", mother);
            } else {
	       part.set().setUserRecord<int>("mother_id", -1);
	    }          
	    numMuonMC++; 
         }
      }
      // fill Gen Electrons passing some basic cuts
      


      // fill Gen Gammas passing some basic cuts
      
   } //end of loop over generated-particles


   
   if (fDebug > 1)  cout << "MC Found:  " << numMuonMC <<  " mu  " << numEleMC << " e  " << numGammaMC << " gamma" << endl;
   EvtView.set().setUserRecord<int>("NumMuon", numMuonMC);
   //EvtView.set().setUserRecord<int>("NumEle", numEleMC);
   //EvtView.set().setUserRecord<int>("NumGamma", numGammaMC);
}

*/

// ------------ reading the Generator Jets ------------


// ------------ reading the Generator MET ------------



// ------------ reading HLT and L1 Trigger Bits ------------

void ePaxAnalyzer::analyzeTrigger(const edm::Event& iEvent, pxl::EventViewRef EvtView) {
  
  //HLT trigger bits
  string errMsg("");
  edm::Handle<edm::TriggerResults> hltresults;
  //HLT producer is called several times within production steps, thus need Input tag with label and process name here
  edm::InputTag hlt = edm::InputTag("TriggerResults","","HLT");
  try {iEvent.getByLabel(hlt, hltresults);} catch (...) { errMsg=errMsg + "  -- No HLTRESULTS";}
  // trigger names
  edm::TriggerNames triggerNames_;

  //HLT: set to false as default
  EvtView.set().setUserRecord<bool>("HLT1Electron", false);
  EvtView.set().setUserRecord<bool>("HLT1ElectronRelaxed", false);
  EvtView.set().setUserRecord<bool>("HLT2Electron", false);
  EvtView.set().setUserRecord<bool>("HLT2ElectronRelaxed", false);
  EvtView.set().setUserRecord<bool>("HLT1EMHighEt", false);
  EvtView.set().setUserRecord<bool>("HLT1EMVeryHighEt", false);
  EvtView.set().setUserRecord<bool>("HLT1MuonIso", false);
  EvtView.set().setUserRecord<bool>("HLT1MuonNonIso", false);
  EvtView.set().setUserRecord<bool>("CandHLT2MuonIso", false);
  EvtView.set().setUserRecord<bool>("HLT2MuonNonIso", false);
  EvtView.set().setUserRecord<bool>("HLTNMuonNonIso", false);
  EvtView.set().setUserRecord<bool>("HLTXElectronMuon", false);
  EvtView.set().setUserRecord<bool>("HLTXElectronMuonRelaxed", false);
  EvtView.set().setUserRecord<bool>("HLTXElectron1Jet", false);
  EvtView.set().setUserRecord<bool>("HLTXElectron2Jet", false);
  EvtView.set().setUserRecord<bool>("HLTXElectron3Jet", false);
  EvtView.set().setUserRecord<bool>("HLTXElectron4Jet", false);
  EvtView.set().setUserRecord<bool>("HLTXMuonJets", false);
  //L1: set to false as default (have all prescale = 1)
  EvtView.set().setUserRecord<bool>("L1_SingleMu7", false);
  EvtView.set().setUserRecord<bool>("L1_DoubleMu3", false);
  EvtView.set().setUserRecord<bool>("L1_SingleIsoEG12", false);
  EvtView.set().setUserRecord<bool>("L1_SingleEG15", false);
  EvtView.set().setUserRecord<bool>("L1_DoubleIsoEG8", false);
  EvtView.set().setUserRecord<bool>("L1_DoubleEG10", false);

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
         vtx.set().setName("PV");
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

   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(fMuonRecoLabel, muonHandle);
   std::vector<pat::Muon> muons = *muonHandle;
   int numMuonRec = 0;


   for( std::vector<pat::Muon>::const_iterator muon = muons.begin(); 
         muon != muons.end(); ++muon ) {
      if (Muon_cuts(muon)) { 
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("Muon");
         part.set().setCharge(muon->charge());
         part.set().vector(pxl::set).setPx(muon->px()); //!!!!!
         part.set().vector(pxl::set).setPy(muon->py());
         part.set().vector(pxl::set).setPz(muon->pz());
         part.set().vector(pxl::set).setE(muon->energy());
         part.set().setUserRecord<double>("Vtx_X", muon->vx());
         part.set().setUserRecord<double>("Vtx_Y", muon->vy());
         part.set().setUserRecord<double>("Vtx_Z", muon->vz()); 

	 //mind that from CMSSW 2.1 on the information will be stored in outerTrack, innerTrack, globalTrack
	 //combinedMuon and standAloneMuon and track will then be deprecated !

	 //save info about quality of track-fit for combined muon (muon system + tracker)
         part.set().setUserRecord<double>("NormChi2", muon->combinedMuon()->normalizedChi2()); //ok
         part.set().setUserRecord<int>("ValidHits", muon->combinedMuon()->numberOfValidHits()); //ok
         part.set().setUserRecord<int>("LostHits", muon->combinedMuon()->numberOfLostHits()); //ok
	 //error info also used in muon-Met corrections, thus store variable to save info for later re-corrections
	 part.set().setUserRecord<double>("dPtRelTrack", muon->combinedMuon()->error(0)/(muon->combinedMuon()->qoverp()) ); //ok
	 part.set().setUserRecord<double>("dPtRelTrack_off", muon->combinedMuon()->ptError()/muon->combinedMuon()->pt() ); //ok
	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part.set().setUserRecord<double>("DSZ", muon->combinedMuon()->dsz()); //ok

	 //save info about quality of track-fit for standAloneMuon (muon system only)
	 part.set().setUserRecord<double>("NormChi2_SA", muon->standAloneMuon()->normalizedChi2()); //ok
         part.set().setUserRecord<int>("ValidHits_SA", muon->standAloneMuon()->numberOfValidHits()); //ok
         part.set().setUserRecord<int>("LostHits_SA", muon->standAloneMuon()->numberOfLostHits()); //ok

	 //save info about quality of track-fit for track (tracker only)
	 part.set().setUserRecord<double>("NormChi2_TM", muon->track()->normalizedChi2()); //ok
         part.set().setUserRecord<int>("ValidHits_TM", muon->track()->numberOfValidHits()); //ok
         part.set().setUserRecord<int>("LostHits_TM", muon->track()->numberOfLostHits()); //ok
	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part.set().setUserRecord<double>("DSZ", muon->combinedMuon()->dsz()); //ok
	 
       
	 //official CaloIso and TrkIso
	 //Def:  aMuon.setCaloIso(aMuon.isolationR03().emEt + aMuon.isolationR03().hadEt + aMuon.isolationR03().hoEt);

	 double CaloIso = muon->caloIso();
	 double TrkIso = muon->trackIso();

	 // TEMPORARY: calculate isolation ourselves
	 //double CaloIso = IsoCalSum(iEvent, 0., muon->track()->eta(), muon->track()->phi(), 0.3, 1.5);
	 //double TrkIso = IsoTrkSum(iEvent,  muon->track()->pt(), muon->track()->eta(), muon->track()->phi(), 0.3, 1.5);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);
 
	 //save offical isolation information
	 const MuonIsolation& muonIsoR03 = muon->isolationR03();
	 part.set().setUserRecord<double>("IsoR03SumPt", muonIsoR03.sumPt);
	 part.set().setUserRecord<double>("IsoR03EmEt", muonIsoR03.emEt);
	 part.set().setUserRecord<double>("IsoR03HadEt", muonIsoR03.hadEt);
	 part.set().setUserRecord<double>("IsoR03HoEt", muonIsoR03.hoEt);
	 part.set().setUserRecord<int>("IsoR03NTracks", muonIsoR03.nTracks);
	 part.set().setUserRecord<int>("IsoR03NJets", muonIsoR03.nJets);
 
	 const MuonIsolation& muonIsoR05 = muon->isolationR05();
	 part.set().setUserRecord<double>("IsoR05SumPt", muonIsoR05.sumPt);
	 part.set().setUserRecord<double>("IsoR05EmEt", muonIsoR05.emEt);
	 part.set().setUserRecord<double>("IsoR05HadEt", muonIsoR05.hadEt);
	 part.set().setUserRecord<double>("IsoR05HoEt", muonIsoR05.hoEt);
	 part.set().setUserRecord<int>("IsoR05NTracks", muonIsoR05.nTracks);
	 part.set().setUserRecord<int>("IsoR05NJets", muonIsoR05.nJets);

	 //save some stuff related to Muon-ID (Calo-info etc.)
	 //Stefan: not found in PAT so far
	 part.set().setUserRecord<double>("CaloCompatibility", muon->caloCompatibility());
	 part.set().setUserRecord<int>("NumberOfChambers", muon->numberOfChambers());
	 part.set().setUserRecord<int>("NumberOfMatches", muon->numberOfMatches());
	 part.set().setUserRecord<double>("MuonDepositEM", muon->calEnergy().em);
	 part.set().setUserRecord<double>("MuonDepositHCAL", muon->calEnergy().had);
	 
         numMuonRec++;
	cout << "numMuonRec=" << numMuonRec <<endl;
      }
   }
   EvtView.set().setUserRecord<int>("NumMuon", numMuonRec);
  
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 

}


// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   int numEleRec = 0;   
   int numEleAll = 0;   // for matching
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;

   edm::Handle<std::vector<pat::Electron> > electronHandle;
   iEvent.getByLabel(fElectronRecoLabel, electronHandle);
   const std::vector<pat::Electron> &electrons = *electronHandle;


//   edm::Handle<SiStripElectronCollection> electrons;
//   edm::Handle<PixelMatchGsfElectronCollection> pixelelectrons;
//   iEvent.getByLabel(fElectronRecoLabel, electrons);
//   iEvent.getByLabel(fPixelMatchElectronRecoLabel, pixelelectrons);

   // Get association maps linking BasicClusters to ClusterShape
   edm::Handle<reco::BasicClusterShapeAssociationCollection> barrelClShpHandle;
   iEvent.getByLabel(fBarrelClusterShapeAssocProducer, barrelClShpHandle);
   edm::Handle<reco::BasicClusterShapeAssociationCollection> endcapClShpHandle;
   iEvent.getByLabel(fEndcapClusterShapeAssocProducer, endcapClShpHandle);

//   int numEleRec = 0;
//   int numPixelEleRec = 0;  
//   int numPixelEleAll = 0; //need counter for EleID 
//   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;

   for ( std::vector<pat::Electron>::const_iterator ele = electrons.begin();
            ele != electrons.end(); ++ele ) {
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
	 part.set().setUserRecord<float>("EoP", ele->eSuperClusterOverP()); //used for CutBasedElectronID
         part.set().setUserRecord<float>("HoEm", ele->hadronicOverEm()); //used for CutBasedElectronID
	 part.set().setUserRecord<float>("DEtaSCVtx", ele->deltaEtaSuperClusterTrackAtVtx()); //used for CutBasedElectronID
	 part.set().setUserRecord<float>("DPhiSCVtx", ele->deltaPhiSuperClusterTrackAtVtx()); //used for CutBasedElectronID
	 //! the seed cluster eta - track eta at calo from outermost state
	 part.set().setUserRecord<float>("DEtaSeedTrk", ele->deltaEtaSeedClusterTrackAtCalo()); //ok
	 part.set().setUserRecord<float>("DPhiSeedTrk", ele->deltaPhiSeedClusterTrackAtCalo()); //ok
	 //part.set().setUserRecord<double>("SCE", ele->superCluster()->energy());
         // the super cluster energy corrected by EnergyScaleFactor
	 part.set().setUserRecord<float>("SCE", ele->caloEnergy()); //ok
         // the errors on the supercluster energy and track momentum
         part.set().setUserRecord<float>("SCEErr", ele->caloEnergyError()); //ok
	 // the errors on the track momentum
         part.set().setUserRecord<float>("PErr", ele->trackMomentumError());	 
         part.set().setUserRecord<double>("TrackerP", ele->gsfTrack()->p()); //ok
	 //the seed cluster energy / track momentum at calo from outermost state
	 part.set().setUserRecord<float>("ESCSeedPout", ele->eSeedClusterOverPout()); //ok
         //part.set().setUserRecord<double>("NormChi2", ele->gsfTrack()->normalizedChi2()); //why was this left out?
         part.set().setUserRecord<int>("ValidHits", ele->gsfTrack()->numberOfValidHits()); //ok
	 part.set().setUserRecord<int>("Class", ele->classification()); //ok

	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part.set().setUserRecord<double>("DSZ", ele->gsfTrack()->dsz()); //ok
	 //part.set().setUserRecord<double>("DZ", ele->gsfTrack()->dz()); //not needed since vz()==dz()
         

	 /*
	 // Get the supercluster (ref) of the Electron
	 // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
	 // a SuperClusterCollection is a std::vector<SuperCluster>
	 // although we get a vector of SuperClusters an electron is only made out of ONE SC
	 // therefore only the first element of the vector should be available!
	 const SuperClusterRef SCRef = ele->superCluster();
	 //const BasicClusterRef& SCSeed = SCRef->seed(); 
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
	 part.set().setUserRecord<double>("EtaEta", seedShapeRef->covEtaEta()); //used for CutBasedElectronID
	 part.set().setUserRecord<double>("EtaPhi", seedShapeRef->covEtaPhi());
	 part.set().setUserRecord<double>("PhiPhi", seedShapeRef->covPhiPhi());

	 //save eta/phi and DetId info from seed-cluster to prevent duplication of Electron/Photon-Candidates (in final selection)
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());

	 //additional data used for cutBasedElectronId in CMSSW 2_0_7
	 part.set().setUserRecord<double>("eSeed", SCRef->seed()->energy()); //used for CutBasedElectronID
	 part.set().setUserRecord<double>("pin", ele->trackMomentumAtVtx().R() ); //used for CutBasedElectronID	 
	 part.set().setUserRecord<double>( "pout", ele->trackMomentumOut().R() ); //used for CutBasedElectronID	
	 
	 */

	 //store ID information
	 //the default leptonID is cut-based for CMSSW_2_0_7 
	 //I have not yet found if "robust" or "tight" cuts are applied probably "robust"(problem is solved in 2_1_0)
	 part.set().setUserRecord<float>("CutBasedID", ele->leptonID());
	 part.set().setUserRecord<float>("CutBasedIDRobust", ele->electronIDRobust());
	 //

	 // TEMPORARY: calculate isolation ourselves
	 //double CaloPt = ( ele->superCluster()->rawEnergy() + ele->superCluster()->rawEnergy()*ele->hadronicOverEm() ) / cosh(ele->trackMomentumAtCalo().eta());
	 //double CaloIso = IsoCalSum(iEvent, CaloPt, ele->trackMomentumAtCalo().eta(), ele->trackMomentumAtCalo().phi(), 0.3, 1.5);
	 //double TrkIso = IsoTrkSum(iEvent, ele->gsfTrack()->pt(), ele->gsfTrack()->eta(), ele->gsfTrack()->phi(), 0.3, 1.5);
	 //double TrkIso = IsoTrkSum(iEvent, ele->trackMomentumAtVtx().Rho(), ele->trackMomentumAtVtx().Eta(), ele->trackMomentumAtVtx().Phi(), 0.2, 0.1);
	 
         //save official isolation information

	 double CaloIso = ele->caloIso();
	 double TrkIso = ele->trackIso();

	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);

	 

	//save additional isolation information(not working for the PAT version corresponding to 2_0_9 I checked out) 
	//not working for > 2_0_7 anyway
	//in CMSSW_2_0_9 the data can be found in lepton.h
	//FIXME after moving to CMSSW_2_0_9 
	//part.set().setUserRecord<float>("egammaTkIso", ele->egammaTkIso());
	//part.set().setUserRecord<int>("egammaTkNumIso", ele->egammaTkNumIso());
	//part.set().setUserRecord<float>("egammaEcalIso", ele->egammaEcalIso());
	//part.set().setUserRecord<float>("egammaHcalIso", ele->egammaHcalIso());

	//FIXME this should somehow be accessible after moving to CMSSW_2_0_9
 	// Get the association vector for track number
	//edm::Handle<reco::PMGsfElectronIsoNumCollection> trackNumHandle;
	//iEvent.getByLabel(fElectronTrackNumProducer, trackNumHandle);

	//direct access (by object). We have candidate already from HCAL-isolation
	//int numVal = (*trackNumHandle)[ electronIsoRef ];
	//part.set().setUserRecord<int>("TrackNum", numVal);


	 	 
         numEleRec++;
      }
      numEleAll++;
   }
  
   EvtView.set().setUserRecord<int>("NumEle", numEleRec);
}






// ------------ reading Reconstructed Jets ------------
//
// Which kind of Jets?
//


// ------------ reading Reconstructed MET ------------


// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventViewRef EvtView) {
   
   // get Photon Collection     
   edm::Handle<std::vector<pat::Photon> > photonHandle;
   iEvent.getByLabel(fGammaRecoLabel, photonHandle);
   std::vector<pat::Photon> photons = *photonHandle;


    
   // get ECAL Cluster shapes
   edm::Handle<reco::BasicClusterShapeAssociationCollection> barrelClShpHandle;
   iEvent.getByLabel(fBarrelClusterShapeAssocProducer, barrelClShpHandle);
   edm::Handle<reco::BasicClusterShapeAssociationCollection> endcapClShpHandle;
   iEvent.getByLabel(fEndcapClusterShapeAssocProducer, endcapClShpHandle);
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
   
   
   int numGammaRec = 0;
   int numGammaAll = 0; //need counter for ID ???


   for (std::vector<pat::Photon>::const_iterator photon = photons.begin(); 
	photon != photons.end(); ++photon) {  
      if ( Gamma_cuts(photon) ) { 
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("Gamma");
         part.set().setCharge(0);
         part.set().vector(pxl::set).setPx(photon->px());
         part.set().vector(pxl::set).setPy(photon->py());
         part.set().vector(pxl::set).setPz(photon->pz());
         part.set().vector(pxl::set).setE(photon->energy());


	 /*	 	
	 // ratio of E(3x3)/ESC
	 part.set().setUserRecord<double>("r9", photon->r9());
	 // ratio of Emax/E(3x3)
	 part.set().setUserRecord<double>("r19", photon->r19());
	 // 5x5 energy
	 part.set().setUserRecord<double>("e5x5", photon->e5x5());
	 /// Whether or not the SuperCluster has a matched pixel seed
	 part.set().setUserRecord<bool>("HasSeed", photon->hasPixelSeed());
	 */
	 /*
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
	 part.set().setUserRecord<double>("e5x5", seedShapeRef->e5x5());
	 //there was a remark that this value equals that from photon->e5x5() but I guess it makes sense to safe this anyway
  	 part.set().setUserRecord<double>("EtaEta", seedShapeRef->covEtaEta());
	 part.set().setUserRecord<double>("EtaPhi", seedShapeRef->covEtaPhi());
	 part.set().setUserRecord<double>("PhiPhi", seedShapeRef->covPhiPhi());


	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection) adn to reject converted photons
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
	 */         

	 
	 //set hadronic over electromagnetic energy fraction
	 part.set().setUserRecord<float>("HoEm", photon->hadronicOverEm());
	 
 	 //save official isolation information

	 part.set().setUserRecord<double>("HCALIso", photon->hcalIso());
	 part.set().setUserRecord<double>("ECALIso", photon->ecalIso());
	 part.set().setUserRecord<double>("TrackIso", photon->trackIso());
	 part.set().setUserRecord<double>("CALIso", photon->caloIso());
	 

	 //FIXME This part does not work. Should this not be possible in an easy way with PAT?
	 //edm::Ref<reco::PhotonCollection> gammaIsoRef(photons, numGammaAll);
	 //reco::CandidateBaseRef candRef(gammaIsoRef);
	 // Get the association vector for track number
	 //edm::Handle<reco::CandViewIntAssociations> trackNumHandle;
	 //iEvent.getByLabel(fGammaTrackNumProducer, trackNumHandle);
	 //direct access (by object). We have candidate already from HCAL-isolation
	 //int numVal = (*trackNumHandle)[ candRef ];
	 //part.set().setUserRecord<int>("TrackNum", numVal);

	 //store information about converted state
	 part.set().setUserRecord<bool>("IsConverted", photon->isConvertedPhoton());  
	 //no way to get this?
	 //part.set().setUserRecord<int>("ConvertedNtrk", Ntrk_conv);  

	 //store photon-Id
	 part.set().setUserRecord<float>("photonID", photon->photonID());

	 /*
	 // store Gamma info corrected for primary vertex (this changes direction but leaves energy of SC unchanged 
	 //get primary vertex (hopefully correct one) for physics eta THIS NEEDS TO BE CHECKED !!!
         math::XYZPoint vtx(0., 0., 0.);
         if (EvtView().findUserRecord<int>("NumVertices") > 0) {
	    pxl::Objects::TypeIterator<pxl::Vertex> iter(EvtView().getObjects()); 
	    pxl::VertexWkPtr pxlvtx = iter.object();
	    if(pxlvtx.valid()) vtx = math::XYZPoint(pxlvtx.get().vector().getX(), pxlvtx.get().vector().getY(), pxlvtx.get().vector().getZ());
         }
	 /////  Set event vertex
	 pat::Photon localPho(*photon);
	 localPho.setVertex(vtx);
	 part.set().setUserRecord<double>("PhysicsEta", localPho.p4().eta());
	 part.set().setUserRecord<double>("PhysicsPhi", localPho.p4().phi());
	*/

	 //FIXME Pi0 stuff still missing
	
	 
         numGammaRec++;
      }	 
      numGammaAll++;
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

/*
// ------------ method to define MC-MUON-cuts

bool ePaxAnalyzer::MuonMC_cuts(const GenParticleCandidate* MCmuon) const {
   //
   if (MCmuon->pt() < 15.) return false;
   if (fabs(MCmuon->eta()) > 3.) return false;
   return true;
}
 


// ------------ method to define MC-Electron-cuts

bool ePaxAnalyzer::EleMC_cuts(const GenParticleCandidate* MCele) const {
   //
   if (MCele->pt() < 15.) return false;
   if (fabs(MCele->eta()) > 3.) return false;
   return true;
}

// ------------ method to define MC-Gamma-cuts

bool ePaxAnalyzer::GammaMC_cuts(const GenParticleCandidate* MCgamma) const {
   //
   if (MCgamma->pt() < 15.) return false;
   if (fabs(MCgamma->eta()) > 3.) return false;
   return true;
}

// ------------ method to define MC-Jet-cuts

bool ePaxAnalyzer::JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const {
   // 
   if (MCjet->pt() < 30.) return false;
   if (fabs(MCjet->eta()) > 3.) return false;
   return true;
}

// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::METMC_cuts(const pxl::ParticleRef MCmet) const {
   // 
   if (MCmet.get().vector().getPt() < 30.) return false;
   return true; 
}

// ------------ method to define RecVertex-cuts

*/
bool ePaxAnalyzer::Vertex_cuts(reco::VertexCollection::const_iterator vertex) const {
  //
  //check compatibility of vertex with beam spot
  double zV = vertex->z();
  double rV = sqrt( vertex->x()*vertex->x() + vertex->y()*vertex->y() );
  if (fabs(zV)>20. || rV>1. ) return false;
  return true;
}

// ------------ method to define MUON-cuts

bool ePaxAnalyzer::Muon_cuts(std::vector<pat::Muon>::const_iterator muon) const {
   if (muon->pt() < 15.)  return false;
   if (fabs(muon->eta()) > 3.) return false;
   return true;
}


// ------------ method to define ELECTRON-cuts

bool ePaxAnalyzer::Ele_cuts(std::vector<pat::Electron>::const_iterator ele) const {
   if (ele->pt() < 15.) return false;
   if (fabs(ele->eta()) > 3.) return false;
   return true;
}

/*
// ------------ method to define ELECTRON-cuts

bool ePaxAnalyzer::Ele_cuts(SiStripElectronCollection::const_iterator ele) const {
   //
   if (ele->pt() < 15.) return false;
   if (fabs(ele->eta()) > 3.) return false;
   return true;
}

// ------------ method to define ELECTRON-cuts

bool ePaxAnalyzer::Ele_cuts(PixelMatchGsfElectronCollection::const_iterator ele) const {
   //
   if (ele->pt() < 15.) return false;
   if (fabs(ele->eta()) > 3.) return false;
   return true;
}

// ------------ method to define JET-cuts

bool ePaxAnalyzer::Jet_cuts(reco::CaloJetCollection::const_iterator jet) const {
   //
   if (jet->pt() < 30.) return false;
   if (fabs(jet->eta()) > 3.) return false;
   return true;
}
*/


// ------------ method to define GAMMA-cuts

bool ePaxAnalyzer::Gamma_cuts(std::vector<pat::Photon>::const_iterator photon) const {
   //
   if (photon->pt() < 15.) return false;
   if (fabs(photon->eta()) > 3.) return false;
   return true;
}

/*
// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::MET_cuts(const pxl::ParticleRef met) const {
   // 
   if (met.get().vector().getPt() < 30.) return false;
   return true;
}


// TEMPORARY STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Stefan: Probably no longer needed because of PAT
//------------------------------------------------------------------------------

double ePaxAnalyzer::IsoCalSum (const edm::Event& iEvent, double ParticleCalPt, double ParticleCalEta, double ParticleCalPhi, double iso_DR, double iso_Seed){
// Computes the sum of Pt inside a cone of R=iso_DR
// using 4-vectors stored in CaloTower objects

  double sum = 0.;

  edm::Handle<CaloTowerCollection> CaloTowerData ;
  iEvent.getByLabel( "towerMaker", CaloTowerData );

  for( CaloTowerCollection::const_iterator tower = CaloTowerData->begin(); 
       tower != CaloTowerData->end(); ++tower ) {
    double eta = tower->eta();
    if ( (tower->energy() / cosh(eta)) > iso_Seed){
      double phi = tower->phi();
      double DR = deltaR(ParticleCalEta, ParticleCalPhi, eta, phi);
      if (DR <= 0.) {DR = 0.001;}
      if (DR < iso_DR){
	double pt = tower->energy() / cosh(eta);
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

  double sum = 0.;

  edm::Handle<TrackCollection> TrackData ;
  iEvent.getByLabel( "ctfWithMaterialTracks", TrackData );

  for( reco::TrackCollection::const_iterator track = TrackData->begin(); 
         track != TrackData->end(); ++track ) {
    if (track->pt() > iso_Seed){
      double eta = track->eta();
      double phi = track->phi();
      double DR = deltaR(ParticleTrkEta, ParticleTrkPhi, eta, phi);
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

  double sum = 0.;

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

      if (p->pt() > iso_Seed){
	double eta = p->eta();
	double phi = p->phi();
	double DR = deltaR(ParticleGenEta, ParticleGenPhi, eta, phi);
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

*/


//define this as a plug-in
DEFINE_FWK_MODULE(ePaxAnalyzer);
