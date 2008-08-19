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

//for GenParticles
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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

//test
#include "DataFormats/Common/interface/Ptr.h"

//test for 2_1_0
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
//#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
 



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
   fgenParticleCandidatesLabel  = iConfig.getUntrackedParameter<string>("genParticleCandidatesLabel");
   //fKtJetMCLabel = iConfig.getUntrackedParameter<string>("KtJetMCLabel");
   fItCone5JetMCLabel = iConfig.getUntrackedParameter<string>("ItCone5JetMCLabel");  
   fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fVertexRecoLabel = iConfig.getUntrackedParameter<string>("VertexRecoLabel");
   //fSAMuonRecoLabel = iConfig.getUntrackedParameter<string>("SAMuonRecoLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   //fBarrelClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("barrelClusterShapeAssociation");
   //fEndcapClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("endcapClusterShapeAssociation"); 
   fbarrelClusterCollection = iConfig.getParameter<edm::InputTag>("barrelClusterCollection");
   fendcapClusterCollection = iConfig.getParameter<edm::InputTag>("endcapClusterCollection");
   freducedBarrelRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
   freducedEndcapRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");
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
   fItCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("ItCone5JetRecoLabel");
   //fL2L3JESic5JetRecoLabel = iConfig.getUntrackedParameter<string>("L2L3JESic5JetRecoLabel");
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");
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


   //create object for EcalClusterLazyTools
   EcalClusterLazyTools lazyTools( iEvent, iSetup, freducedBarrelRecHitCollection, freducedEndcapRecHitCollection);

   // Generator stuff
   analyzeGenInfo(iEvent, GenEvtView);
   analyzeGenJets(iEvent, GenEvtView);
   //if( fGenOnly == false ){ //only if info is in event
     analyzeGenMET(iEvent, GenEvtView);

     //Trigger bits
     analyzeTrigger(iEvent, RecEvtView);
   
   // Reconstructed stuff
     analyzeRecVertices(iEvent, RecEvtView);
     analyzeRecMuons(iEvent, RecEvtView);
     analyzeRecElectrons(iEvent, RecEvtView, lazyTools);
     analyzeRecJets(iEvent, RecEvtView);
     analyzeRecMET(iEvent, RecEvtView);
     analyzeRecGammas(iEvent, RecEvtView, lazyTools);


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




// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   //gen particles
   edm::Handle<reco::GenParticleCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );
   
   //TrackingVertexCollection ONLY available in RECO, so need ugly catch workaround...
/*
   //try{
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
	*/
     const GenParticle* p = (const GenParticle*) &(*genParticleHandel->begin()); //this is the incoming proton
     pxl::VertexRef GenVtx = EvtView.set().create<pxl::Vertex>();
     GenVtx.set().setName("PV");
     //if(p->daughter(0)){ //check validity, otherwise sometimes segmentation violation - WHY ???
       GenVtx.set().vector(pxl::set).setX(p->daughter(0)->vx()); //need daughter since first particle (proton) has position zero
       GenVtx.set().vector(pxl::set).setY(p->daughter(0)->vy());
       GenVtx.set().vector(pxl::set).setZ(p->daughter(0)->vz());
     //}else{ //set dummy values to catch segmentation violation
     //  GenVtx.set().vector(pxl::set).setX(0.); //need daughter since first particle (proton) has position zero
     //  GenVtx.set().vector(pxl::set).setY(0.);
     //  GenVtx.set().vector(pxl::set).setZ(0.);
     //}
     EvtView.set().setUserRecord<int>("NumVertices", 1);  
	     
   //}

  
   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;
   int GenId = 0;

   //save mother of stable particle
   const Candidate* p_mother; 
   int mother = 0;
   
   // loop over all particles
   for( reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); 
	pa != genParticleHandel->end(); ++ pa ) {

     //cast iterator into GenParticleCandidate
     const GenParticle* p = (const GenParticle*) &(*pa);
 
  
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
	    part.set().setUserRecord<int>("GenId", GenId);
	   
	    // TEMPORARY: calculate isolation ourselves 
	    //FIXME: make this at least comparable with pat/lepton isolation
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
	    part.set().setUserRecord<int>("GenId", GenId);
	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, p->pt(), p->eta(), p->phi(), 0.3, 1.5);
	    part.set().setUserRecord<double>("GenIso", GenIso);
	    //save mother of stable electron
	    p_mother =p->mother(); 
	    if (p_mother != 0) {
	       mother = p_mother->pdgId();
	       //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	       while( abs(mother) == 11 ){
	         p_mother = p_mother->mother();
	         mother = p_mother->pdgId();
	       }	       
	       part.set().setUserRecord<int>("mother_id", mother);
            } else {
	       part.set().setUserRecord<int>("mother_id", -1);
	    }          
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
	    part.set().setUserRecord<int>("GenId", GenId);
	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, 0., p->eta(), p->phi(), 0.3, 1.5); //0. since gamma not charged!
	    part.set().setUserRecord<double>("GenIso", GenIso);
	    
	    //save mother of stable gamma
	    p_mother =p->mother(); 
	    if (p_mother != 0) {
	       mother = p_mother->pdgId();
	       //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	       while( abs(mother) == 22 ){
	         p_mother = p_mother->mother();
	         mother = p_mother->pdgId();
	       }	       
	       part.set().setUserRecord<int>("mother_id", mother);
            } else {
	       part.set().setUserRecord<int>("mother_id", -1);
	    }          
	    numGammaMC++;
         }
      }

      GenId++;
   } //end of loop over generated-particles

   if (fDebug > 1)  cout << "MC Found:  " << numMuonMC <<  " mu  " << numEleMC << " e  " << numGammaMC << " gamma" << endl;
   EvtView.set().setUserRecord<int>("NumMuon", numMuonMC);
   EvtView.set().setUserRecord<int>("NumEle", numEleMC);
   EvtView.set().setUserRecord<int>("NumGamma", numGammaMC);
}



// ------------ reading the Generator Jets ------------


void ePaxAnalyzer::analyzeGenJets(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   //Get the GenJet collections
   //edm::Handle<reco::GenJetCollection> Kt_GenJets;
   edm::Handle<reco::GenJetCollection> ItCone5_GenJets;
   //iEvent.getByLabel(fKtJetMCLabel, Kt_GenJets);
   iEvent.getByLabel(fItCone5JetMCLabel, ItCone5_GenJets);
   
   //int numKtJetMC = 0;
   int numItCone5JetMC = 0;
 
   vector <const GenParticle*> genJetConstit; //genJet-constituents
   int numGenJetConstit_withcuts = 0;
   double constit_pT = 5.; //here we have a hardcoded cut, but do we really need cfg-parameter for this?...
   
   

   //Loop over ItCone5GenJets
   for( reco::GenJetCollection::const_iterator genJet = ItCone5_GenJets->begin();
         genJet != ItCone5_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("ItCone5Jet");
         part.set().vector(pxl::set).setPx(genJet->px());
         part.set().vector(pxl::set).setPy(genJet->py());
         part.set().vector(pxl::set).setPz(genJet->pz());
         part.set().vector(pxl::set).setE((genJet->energy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
	 part.set().setUserRecord<double>("AuxE", genJet->auxiliaryEnergy());
         numItCone5JetMC++;

	 //save number of GenJet-constituents fulfilling some cuts
	 numGenJetConstit_withcuts = 0; //reset variable
	 genJetConstit = genJet->getGenConstituents();
	 for( std::vector<const GenParticle*>::iterator constit = genJetConstit.begin(); 
	     constit != genJetConstit.end(); ++constit ) {
	    
	   //raise counter if cut passed
	   if( (*constit)->pt() > constit_pT ) numGenJetConstit_withcuts++; 
	
	 }//end-of loop over all constituents
	
	 part.set().setUserRecord<int>("GenJetConstit", numGenJetConstit_withcuts);

      }

   }
   EvtView.set().setUserRecord<int>("NumItCone5Jet", numItCone5JetMC);
   
  
   if (fDebug > 1) cout << "Found MC Jets:  "  << numItCone5JetMC << " It5  " << endl;


   
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

  
   if (METMC_cuts(part)) { 
     numMETMC++; 
   }
   EvtView.set().setUserRecord<int>("NumMET", numMETMC);
   if (numMETMC && fDebug > 1) cout << "Event contains MET" << endl;

 
}





// ------------ reading the Reconstructed MET ------------
//the stored information should already contain the muon corrections plus several other corrections
//it will have certainly corrections in the future as there are certain functions for uncorrection planned

void ePaxAnalyzer::analyzeRecMET(const edm::Event& iEvent, pxl::EventViewRef EvtView) {

   edm::Handle<std::vector<pat::MET> > METHandle;
   iEvent.getByLabel(fMETRecoLabel, METHandle);
   std::vector<pat::MET> METs = *METHandle;
   std::vector<pat::MET>::const_iterator met = METs.begin();
   

   int numMETRec = 0;
   pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
   part.set().setName("MET");
   part.set().vector(pxl::set).setPx(met->px());
   part.set().vector(pxl::set).setPy(met->py());
   part.set().vector(pxl::set).setPz(0.);
   part.set().vector(pxl::set).setMass(0.);
   part.set().setUserRecord<double>("sumEt", met->sumEt());
   part.set().setUserRecord<double>("mEtSig", met->mEtSig());
   part.set().setUserRecord<double>("EmEt", met->emEtFraction());         
   part.set().setUserRecord<double>("HadEt", met->etFractionHadronic()); 
   part.set().setUserRecord<double>("MaxEtEm", met->maxEtInEmTowers());   
   part.set().setUserRecord<double>("MaxEtHad", met->maxEtInHadTowers()); 


   if (MET_cuts(part)) { 
     numMETRec++;
   }

   EvtView.set().setUserRecord<int>("NumMET", numMETRec);


}


// ------------ reading HLT and L1 Trigger Bits ------------

void ePaxAnalyzer::analyzeTrigger(const edm::Event& iEvent, pxl::EventViewRef EvtView) {
  /*
  //HLT trigger bits
  string errMsg("");
  edm::Handle<edm::TriggerResults> hltresults;
  //HLT producer is called several times within production steps, thus need Input tag with label and process name here
  edm::InputTag hlt = edm::InputTag("TriggerResults","","HLT");
  try {iEvent.getByLabel(hlt, hltresults);} catch (...) { errMsg=errMsg + "  -- No HLTRESULTS";}
  // trigger names
  edm::TriggerNames triggerNames_;
 */

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

//gen particles for PAT-matching
   edm::Handle<reco::GenParticleCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );


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

	//store PAT matching info
	int count = 0; // keeping track of n-th loop
        const reco::Particle* matched = muon->genLepton();
	for( reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); 
	pa != genParticleHandel->end(); ++ pa ) {
		if( &(*pa) == matched ) {
		cout << "MATCH!" << endl;
		double ptrec = muon->pt();
		double ptgen = pa->pt(); 
		cout << "pt des gen-mu: " << ptgen << endl;
		cout << "pt des rec-mu: " << ptrec<< endl;
	part.set().setUserRecord<double>("MatchedGenId", count); 
			count++;
		}
		else{
			part.set().setUserRecord<double>("MatchedGenId", -1);
			count++;
		}
	}
	 


	 //mind that from CMSSW 2.1 on information will be stored in outerTrack, innerTrack, globalTrack
	 //combinedMuon and standAloneMuon and track might then be deprecated sooner or later!

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
	 

	//Variables for outerTrack, innerTrack, globalTrack ...

	//FIXME variables containing information about muon quality ...


       
	 //official CaloIso and TrkIso
	 //Def:  aMuon.setCaloIso(aMuon.isolationR03().emEt + aMuon.isolationR03().hadEt + aMuon.isolationR03().hoEt);

	 double CaloIso = muon->caloIso();
	 double TrkIso = muon->trackIso();

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
	 
	 part.set().setUserRecord<double>("CaloCompatibility", muon->caloCompatibility());
	 part.set().setUserRecord<int>("NumberOfChambers", muon->numberOfChambers());
	 part.set().setUserRecord<int>("NumberOfMatches", muon->numberOfMatches());
	 part.set().setUserRecord<double>("MuonDepositEM", muon->calEnergy().em);
	 part.set().setUserRecord<double>("MuonDepositHCAL", muon->calEnergy().had);
	 
         numMuonRec++;
	//cout << "numMuonRec=" << numMuonRec <<endl;
      }
   }
   EvtView.set().setUserRecord<int>("NumMuon", numMuonRec);
  
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 

}


// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, pxl::EventViewRef EvtView, EcalClusterLazyTools& lazyTools) {

   int numEleRec = 0;   
   int numEleAll = 0;   // for matching
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr; //wegen 2_1_0
BasicClusterRefVector::iterator basicCluster_iterator;


   edm::Handle<std::vector<pat::Electron> > electronHandle;
   iEvent.getByLabel(fElectronRecoLabel, electronHandle);
   const std::vector<pat::Electron> &electrons = *electronHandle;


//gen particles for PAT-matching
   edm::Handle<reco::GenParticleCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );


   // Get association maps linking BasicClusters to ClusterShape FIXME make this work in 2_1_0
   edm::Handle<reco::BasicClusterShapeAssociationCollection> barrelClShpHandle;
   iEvent.getByLabel(fBarrelClusterShapeAssocProducer, barrelClShpHandle);
   edm::Handle<reco::BasicClusterShapeAssociationCollection> endcapClShpHandle;
   iEvent.getByLabel(fEndcapClusterShapeAssocProducer, endcapClShpHandle);


   for ( std::vector<pat::Electron>::const_iterator ele = electrons.begin();
            ele != electrons.end(); ++ele ) {
      if (Ele_cuts(ele)) {
         if (fDebug > 1) {
	    cout << "Electron Energy scale corrected: " << ele->isEnergyScaleCorrected() 
	         << "  Momentum corrected: " << ele->isMomentumCorrected() << endl;
         }


	 pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();



//pxl::ParticleFilter recElectronList(recEvtView().getObjects(), kElectronName);
	//for (pxl::ParticleFilterIterator iter(recElectronList); !iter.isDone(); iter.next()) {
	//	pxl::ParticleWkPtr recPa = iter.wkPtr();
	//	int id = recPa.get().findUserRecord<int>(kCandidateId, -1);
		
		//edm::RefToBase<pat::ElectronType> eleRef = ele->originalObjectRef();
  		//const edm::Ptr<reco::Candidate> & eleRef = ele->originalObjectRef();
		//reco::GenParticleRef genEle = (*genMatch)[eleRef];
		
		//if (genEle.isNonnull() && genEle.isAvailable() ) {
			//cout << "rec pt " << muoeleRef->pt() << " matched gen particle pt= " << genMuon->pt() << " key " << genMuon.key() << endl;
			//part->set().setUserRecord<int>("MatchedKey", genEle.key());
		//}
		
	//}


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
	 part.set().setUserRecord<double>("DEtaSeedTrk", ele->deltaEtaSeedClusterTrackAtCalo()); //ok
	 part.set().setUserRecord<double>("DPhiSeedTrk", ele->deltaPhiSeedClusterTrackAtCalo()); //ok
	 //part.set().setUserRecord<double>("SCE", ele->superCluster()->energy());
         // the super cluster energy corrected by EnergyScaleFactor
	 part.set().setUserRecord<float>("SCE", ele->caloEnergy()); //ok
         // the errors on the supercluster energy and track momentum
         part.set().setUserRecord<float>("SCEErr", ele->caloEnergyError()); //ok
	 // the errors on the track momentum
         part.set().setUserRecord<float>("PErr", ele->trackMomentumError());	 
         part.set().setUserRecord<double>("TrackerP", ele->gsfTrack()->p()); //ok
	 //the seed cluster energy / track momentum at calo from outermost state
	 part.set().setUserRecord<double>("ESCSeedPout", ele->eSeedClusterOverPout()); //ok
         //part.set().setUserRecord<double>("NormChi2", ele->gsfTrack()->normalizedChi2()); //why was this left out?
         part.set().setUserRecord<int>("ValidHits", ele->gsfTrack()->numberOfValidHits()); //ok
	 part.set().setUserRecord<int>("Class", ele->classification()); //ok

	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part.set().setUserRecord<double>("DSZ", ele->gsfTrack()->dsz()); //ok
	 //part.set().setUserRecord<double>("DZ", ele->gsfTrack()->dz()); //not needed since vz()==dz()
         

	//store PAT matching info
	int count = 0; // keeping track of n-th loop
        const reco::Particle* matched = ele->genLepton();
       	for( reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); 
	pa != genParticleHandel->end(); ++ pa ) {
		if( &(*pa) == matched ) {
			cout << "MATCH!" << endl;
			double ptrec = ele->pt();
			double ptgen = pa->pt(); 
			cout << "pt des gen-el: " << ptgen << endl;
			cout << "pt des rec-el: " << ptrec<< endl;
		 	part.set().setUserRecord<double>("MatchedGenId", count); 
			count++;
		}
		else{
			part.set().setUserRecord<double>("MatchedGenId", -1);
			count++;
		}
	}
	 


	 // Get the supercluster (ref) of the Electron
	 // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
	 // a SuperClusterCollection is a std::vector<SuperCluster>
	 // although we get a vector of SuperClusters an electron is only made out of ONE SC
	 // therefore only the first element of the vector should be available!
	 const SuperClusterRef SCRef = ele->superCluster();
	 const BasicClusterRef& SCSeed = SCRef->seed(); 
         // Find the entry in the map corresponding to the seed BasicCluster of the SuperCluster

         //use EcalClusterLazyTools to store ClusterShapeVariables

	 part.set().setUserRecord<double>("e3x3",  lazyTools.e3x3(*SCSeed) );
	 part.set().setUserRecord<double>("e5x5",  lazyTools.e5x5(*SCSeed)  );
         std::vector<float> covariances = lazyTools.covariances(*SCSeed, 4.7 );
	 part.set().setUserRecord<double>("EtaEta", covariances[0] ); //used for CutBasedElectronID
	 part.set().setUserRecord<double>("EtaPhi", covariances[1] );
	 part.set().setUserRecord<double>("PhiPhi", covariances[2] );
	 

	 //save eta/phi and DetId info from seed-cluster to prevent duplication of Electron/Photon-Candidates (in final selection)
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 //part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());

	 
	 //additional data used for cutBasedElectronId in CMSSW 2_0_7
	 part.set().setUserRecord<double>("eSeed", SCRef->seed()->energy()); //used for CutBasedElectronID
	 part.set().setUserRecord<double>("pin", ele->trackMomentumAtVtx().R() ); //used for CutBasedElectronID	 
	 part.set().setUserRecord<double>( "pout", ele->trackMomentumOut().R() ); //used for CutBasedElectronID	
	 
	 

	 //store ID information
	 //Cut based ID is stored as float in 2_1_0 hence it is converted to bool for backwards compability

	float IDfloat =  ele->leptonID("tight");
	bool IDbool = false ;
	if (IDfloat > 0.5) {IDbool = true;}
	part.set().setUserRecord<bool>("CutBasedIDTight", IDbool);
	IDfloat = ele->leptonID("robust");
	IDbool = false;
	if (IDfloat > 0.5) {IDbool = true;}	
	part.set().setUserRecord<bool>("CutBasedIDRobust", IDbool);
	 

	 
         //save official isolation information

	 double CaloIso = ele->caloIso();
	 double TrkIso = ele->trackIso();

	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrackIso", TrkIso);

	 

	//save additional isolation information(not working for the PAT version corresponding to 2_0_9 I checked out) 
	//not working for > 2_0_7 anyway
	//in CMSSW_2_0_9 the data can be found in lepton.h
	//FIXME after moving to CMSSW_2_0_9 
	//part.set().setUserRecord<float>("egammaTkIso", ele->egammaTkIso());
	//part.set().setUserRecord<int>("egammaTkNumIso", ele->egammaTkNumIso());
	part.set().setUserRecord<double>("ECALIso", ele->ecalIso());
	part.set().setUserRecord<double>("HCALIso", ele->hcalIso());

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
// by default for creating the allLayer0Jets the iterativeCone5CaloJets are used
// discuss if we also need informations from other Jet-Collections

void ePaxAnalyzer::analyzeRecJets(const edm::Event& iEvent, pxl::EventViewRef EvtView) 
{
   int numItCone5JetRec = 0;



   //get primary vertex (hopefully correct one) for physics eta
   double VertexZ = 0.;
   if (EvtView().findUserRecord<int>("NumVertices") > 0) {
      pxl::Objects::TypeIterator<pxl::Vertex> iter(EvtView().getObjects()); 
      pxl::VertexWkPtr vtx = iter.object();
      if(vtx.valid()) VertexZ = vtx.get().vector().getZ();
   } 

   
   
   edm::Handle<std::vector<pat::Jet> > jetHandle;
   iEvent.getByLabel(fItCone5JetRecoLabel, jetHandle);
   std::vector<pat::Jet> jets = *jetHandle;
   
//Get the GenJet collections for PAT matching
   //edm::Handle<reco::GenJetCollection> Kt_GenJets;
   edm::Handle<reco::GenJetCollection> ItCone5_GenJets;
   //iEvent.getByLabel(fKtJetMCLabel, Kt_GenJets);
   iEvent.getByLabel(fItCone5JetMCLabel, ItCone5_GenJets);


   for(std::vector<pat::Jet>::const_iterator jet = jets.begin();
           jet != jets.end(); ++jet ) 
   {
   	
      if (Jet_cuts(jet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("ItCone5Jet"); //compare with Clemens; this may need further discussion
         part.set().vector(pxl::set).setPx(jet->px());
         part.set().vector(pxl::set).setPy(jet->py());
         part.set().vector(pxl::set).setPz(jet->pz());
         part.set().vector(pxl::set).setE(jet->energy());
	 part.set().setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part.set().setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part.set().setUserRecord<int>("N90", jet->n90());
	 part.set().setUserRecord<int>("N60", jet->n60());
	 std::vector <CaloTowerPtr> caloRefs = jet->getCaloConstituents();
	 part.set().setUserRecord<int>("NCaloRefs", caloRefs.size());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
 	 part.set().setUserRecord<double>("TowersArea", jet->towersArea());
	 part.set().setUserRecord<double>("PhysicsEta", jet->physicsEta(VertexZ,jet->eta()));


	//store PAT matching info
	//FIXME does not work, because implementation of genJet() differs from genLepton()
	int count = 0;
        const reco::Jet* matched = jet->genJet();
       	for( reco::GenJetCollection::const_iterator pa = ItCone5_GenJets->begin(); 
	pa != ItCone5_GenJets->end(); ++ pa ) {
		if( &(*pa) == matched ) {
		cout << "MATCH!" << endl;
		double ptrec = jet->pt();
		double ptgen = pa->pt(); 
		cout << "pt des gen-jets: " << ptgen << endl;
		cout << "pt des rec-jets: " << ptrec<< endl;
		part.set().setUserRecord<double>("MatchedGenId", count); 
			count++;
		}
		else{
			part.set().setUserRecord<double>("MatchedGenId", -1);
			count++;
		}
	}
	 
      		//if (jet->isCaloJet()) { cout << "isCaloJet!!!!!!!!!!!" << endl;}
	 		//else { cout << "NOT isCaloJet!!!!!!!!!!!" << endl;}

         numItCone5JetRec++;
	 
     }
						
   }
   EvtView.set().setUserRecord<int>("NumItCone5Jet", numItCone5JetRec);

   if (fDebug > 1) cout << "Found Rec Jets:  " << numItCone5JetRec << " It5  "
                        << endl;

}


// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventViewRef EvtView, EcalClusterLazyTools& lazyTools){
   
   // get Photon Collection     
   edm::Handle<std::vector<pat::Photon> > photonHandle;
   iEvent.getByLabel(fGammaRecoLabel, photonHandle);
   std::vector<pat::Photon> photons = *photonHandle;

   
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


	 /// Whether or not the SuperCluster has a matched pixel seed
	 part.set().setUserRecord<bool>("HasSeed", photon->hasPixelSeed());
	 
	 
	 const SuperClusterRef SCRef = photon->superCluster();
         // Find the entry in the map corresponding to the seed BasicCluster of the SuperCluster
         const BasicClusterRef& SCSeed = SCRef->seed();

	 part.set().setUserRecord<double>("rawEnergy",  SCRef->rawEnergy() );
 	 part.set().setUserRecord<double>("preshowerEnergy",  SCRef->preshowerEnergy() );

	


	/*	 
         DetId id = SCRef->seed()->getHitsByDetId()[0];
         if (id.subdetId() == EcalBarrel) {
		//cout << "ECALBarrel TRUE !!!!!" << endl;
            //seedShpItr = barrelClShpHandle->find(SCRef->seed()); //FIXME 2_1_0
         } else {
		//cout << "NOT ECALBarrel TRUE !!!!" << endl;
            //seedShpItr = endcapClShpHandle->find(SCRef->seed()); //FIXME 2_1_0
         }
	*/

         //use EcalClusterLazyTools to store ClusterShapeVariables
	 part.set().setUserRecord<double>("e3x3",  lazyTools.e3x3(*SCSeed) );
	 part.set().setUserRecord<double>("e5x5",  lazyTools.e5x5(*SCSeed)  );
         std::vector<float> covariances = lazyTools.covariances(*SCSeed );
	 part.set().setUserRecord<double>("EtaEta", covariances[0] ); 
	 part.set().setUserRecord<double>("EtaPhi", covariances[1] );
	 part.set().setUserRecord<double>("PhiPhi", covariances[2] );
	 part.set().setUserRecord<double>("Emax",  lazyTools.eMax(*SCSeed)  );

	 part.set().setUserRecord<double>("r9", ( lazyTools.e3x3(*SCSeed) )/( SCRef->rawEnergy() + SCRef->preshowerEnergy() ) );
	 part.set().setUserRecord<double>("r19",  (lazyTools.eMax(*SCSeed) / lazyTools.e3x3(*SCSeed)) );

	 

	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection) adn to reject converted photons
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 //part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
	        
	 
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

	 //FIXME store photon-Id, make sure this is what we want!
	 part.set().setUserRecord<float>("photonIDTight", photon->isTightPhoton());

	 
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
	 //localPho.setVertex(vtx);  //FIXME this line does not work (missing Cluster Shape info)
	 part.set().setUserRecord<double>("PhysicsEta", localPho.p4().eta());
	 part.set().setUserRecord<double>("PhysicsPhi", localPho.p4().phi());
	 part.set().setUserRecord<double>("PhysicsPt", localPho.p4().pt());
	 
	 //FIXME Pi0 stuff still missing
	
	 
         numGammaRec++;
      }	 
      numGammaAll++;
   }

   EvtView.set().setUserRecord<int>("NumGamma", numGammaRec);
   if (fDebug > 1) cout << "Rec Gamma: " << numGammaRec << endl;
  
}






// ------------ method returning the EventClassType ------------

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

bool ePaxAnalyzer::MuonMC_cuts(const GenParticle* MCmuon) const {
   //
   if (MCmuon->pt() < 15.) return false;
   if (fabs(MCmuon->eta()) > 3.) return false;
   return true;
}
 


// ------------ method to define MC-Electron-cuts

bool ePaxAnalyzer::EleMC_cuts(const GenParticle* MCele) const {
   //
   if (MCele->pt() < 15.) return false;
   if (fabs(MCele->eta()) > 3.) return false;
   return true;
}

// ------------ method to define MC-Gamma-cuts

bool ePaxAnalyzer::GammaMC_cuts(const GenParticle* MCgamma) const {
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

// ------------ method to define JET-cuts

bool ePaxAnalyzer::Jet_cuts(std::vector<pat::Jet>::const_iterator jet) const {
   //
   if (jet->pt() < 30.) return false;
   if (fabs(jet->eta()) > 3.) return false;
   return true;
}


// ------------ method to define GAMMA-cuts

bool ePaxAnalyzer::Gamma_cuts(std::vector<pat::Photon>::const_iterator photon) const {
   //
   if (photon->pt() < 15.) return false;
   if (fabs(photon->eta()) > 3.) return false;
   return true;
}


// ------------ method to define MET-cuts

bool ePaxAnalyzer::MET_cuts(const pxl::ParticleRef met) const {
   // 
   if (met.get().vector().getPt() < 30.) return false;
   return true;
}

/*
// TEMPORARY STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Stefan: no longer needed because of PAT
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
*/
//FIXME compare to PAT-isolation 
double ePaxAnalyzer::IsoGenSum (const edm::Event& iEvent, double ParticleGenPt, double ParticleGenEta, double ParticleGenPhi, double iso_DR, double iso_Seed){
// Computes the sum of Pt inside a cone of R=iso_DR
// using 4-vectors stored in GenParticle objects

  double sum = 0.;

  //gen particles
  edm::Handle<reco::GenParticleCollection> genParticleHandel;
  iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );

  // loop over all particles
  for( reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); 
       pa != genParticleHandel->end(); ++ pa ) {

    //cast iterator into GenParticleCandidate
    const GenParticle* p = (const GenParticle*) &(*pa);

    // only consider stable particles and charged particles in order to be more comparable with track-isolation
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



//define this as a plug-in
DEFINE_FWK_MODULE(ePaxAnalyzer);
