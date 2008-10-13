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
#include <map>
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
ePaxAnalyzer::ePaxAnalyzer(const edm::ParameterSet& iConfig) : fFileName(iConfig.getUntrackedParameter<string>("FileName")), fePaxFile(fFileName) {
   //now do what ever initialization is needed
   // Get Physics process
   fProcess = iConfig.getUntrackedParameter<string>("Process");
   // Gen-Only or also Rec-information
   fGenOnly = iConfig.getUntrackedParameter<bool>("GenOnly");
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
   // Owner of all Pxl Objects 
   pxl::Event* event = new pxl::Event();

   // event-specific data
   event->setUserRecord<bool>("Data", false);  //distinguish between data and MC
   event->setUserRecord<int>("Run", iEvent.id().run());
   event->setUserRecord<int>("ID", iEvent.id().event());	
   //cout << "Run " << iEvent.id().run() << "   EventID = " << iEvent.id().event() << endl;  

   // create two ePaxEventViews for Generator/Reconstructed Objects
   pxl::EventView* GenEvtView = event->createIndexed<pxl::EventView>("Gen");
   pxl::EventView* RecEvtView = event->createIndexed<pxl::EventView>("Rec");
   GenEvtView->setUserRecord<std::string>("Type", "Gen");
   RecEvtView->setUserRecord<std::string>("Type", "Rec");
   

   double processID = 0.;
   double pthat = 0.;
   
   // else use stuff inside AOD : can it be done in PAT?
  
      //edm::Handle< int > genProcessID;
      //iEvent.getByLabel( "genEventProcID", genProcessID );
      //processID = *genProcessID;
 
      //edm::Handle< double > genEventScale;
      //iEvent.getByLabel( "genEventScale", genEventScale );
      //pthat = *genEventScale;
    

   //maps for matching
	std::map<const Particle*, pxl::Particle*> genmap;

   //set process name
   GenEvtView->setUserRecord<std::string>("Process", fProcess);
   RecEvtView->setUserRecord<std::string>("Process", fProcess);
   // store the ID in both event views? 
   GenEvtView->setUserRecord<double>("pthat", pthat);
   RecEvtView->setUserRecord<double>("pthat", pthat); 
   
   double weight = 1.;
   //cout << "Weight = " << weight << endl;
   //cout << "Process ID: " << processID << " and Event Scale (pthat): " << pthat << endl;
   GenEvtView->setUserRecord<double>("Weight", weight);
   RecEvtView->setUserRecord<double>("Weight", weight);
   // store Process ID
   GenEvtView->setUserRecord<double>("ProcID", processID);
   RecEvtView->setUserRecord<double>("ProcID", processID);
   //create object for EcalClusterLazyTools
   EcalClusterLazyTools lazyTools( iEvent, iSetup, freducedBarrelRecHitCollection, freducedEndcapRecHitCollection);
   // Generator stuff
   analyzeGenInfo(iEvent, GenEvtView, genmap);
   analyzeGenJets(iEvent, GenEvtView);
   // store Rec Objects only if requested
   if (!fGenOnly) {
      analyzeGenMET(iEvent, GenEvtView);
      //Trigger bits
      analyzeTrigger(iEvent, RecEvtView);
      // Reconstructed stuff
      analyzeRecVertices(iEvent, RecEvtView);
      analyzeRecMuons(iEvent, RecEvtView, GenEvtView, genmap);
      analyzeRecElectrons(iEvent, RecEvtView, GenEvtView, lazyTools, genmap);
      analyzeRecJets(iEvent, RecEvtView);
      analyzeRecMET(iEvent, RecEvtView);
      analyzeRecGammas(iEvent, RecEvtView, GenEvtView, lazyTools, genmap);
   }

   // set event class strings
   GenEvtView->setUserRecord<std::string>("EventClass", getEventClass(GenEvtView));
   RecEvtView->setUserRecord<std::string>("EventClass", getEventClass(RecEvtView));
   
   if (fDebug > 0) {
      cout << "UserRecord  " <<  "   e   " << "  mu   " << " Gamma " << " KtJet " << "  MET  " << endl;
      cout << "Found (MC): " << setw(4) << GenEvtView->findUserRecord<int>("NumEle", 0) 
           << setw(7) << GenEvtView->findUserRecord<int>("NumMuon", 0)
           << setw(7) << GenEvtView->findUserRecord<int>("NumGamma", 0) 
           << setw(4) << GenEvtView->findUserRecord<int>("NumKtJet", 0) << "/" 
           << GenEvtView->findUserRecord<int>("NumItCone5Jet", 0) << "/" 
           << GenEvtView->findUserRecord<int>("NumMidCone5Jet", 0)
           << setw(7) << GenEvtView->findUserRecord<int>("NumMET", 0) << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView->findUserRecord<int>("NumEle", 0)    
           << setw(7) << RecEvtView->findUserRecord<int>("NumMuon", 0) 
           << setw(7) << RecEvtView->findUserRecord<int>("NumGamma", 0)
           << setw(4) << RecEvtView->findUserRecord<int>("NumKtJet", 0) << "/"
           << RecEvtView->findUserRecord<int>("NumItCone5Jet", 0) 
           << setw(7) << RecEvtView->findUserRecord<int>("NumMET", 0) << endl;
   }

   if (fDebug > 0) { 
      cout << "Gen Event Type: " << GenEvtView->findUserRecord<string>("EventClass") << endl;
      cout << "Rec Event Type: " << RecEvtView->findUserRecord<string>("EventClass") << endl;
   }   
   fePaxFile.writeEvent(event, RecEvtView->findUserRecord<string>("EventClass"));
}




// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, pxl::EventView* EvtView, std::map<const Particle*, pxl::Particle*> & genmap ) {

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
     pxl::VertexRef GenVtx = EvtView->create<pxl::Vertex>();
     GenVtx.setName("PV");
	
     GenVtx.set().vector(pxl::set).setX(EventVertices->position().x());
     GenVtx.set().vector(pxl::set).setY(EventVertices->position().y());
     GenVtx.set().vector(pxl::set).setZ(EventVertices->position().z());
     // we only have a single PV at Generator Level. Due to EventView Consistency .i.e. GenEvtView should look identical to RecEvtView
     // this variable is explicetly set
     EvtView->setUserRecord<int>("NumVertices", 1);               
     // do we need this BX/event identification???
     //GenVtx.setUserRecord<int>("Vtx_BX", EventVertices->eventId().bunchCrossing());
     //GenVtx.setUserRecord<int>("Vtx_event", EventVertices->eventId().event());
     
   }catch(...){ //for AOD use daughter of first GenParticle, this should be equal to generated primary vertex...  
	*/
     const GenParticle* p = (const GenParticle*) &(*genParticleHandel->begin()); //this is the incoming proton
     pxl::Vertex* GenVtx = EvtView->create<pxl::Vertex>();
     GenVtx->setName("PV");
     GenVtx->setXYZ(p->daughter(0)->vx(), p->daughter(0)->vy(), p->daughter(0)->vz()); //need daughter since first particle (proton) has position zero
     EvtView->setUserRecord<int>("NumVertices", 1);  
	     
   //}
  
   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;
   int GenId = 0;

   //save mother of stable particle
   const Candidate* p_mother; 
   int mother = 0;
   
   // loop over all particles
   for (reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); pa != genParticleHandel->end(); ++ pa ) {

      //cast iterator into GenParticleCandidate
      const GenParticle* p = (const GenParticle*) &(*pa);
      // fill Gen Muons passing some basic cuts
      if ( abs((p)->pdgId()) == 13 && (p)->status() == 1) {
         if ( MuonMC_cuts(p) ) { 
            pxl::Particle* part = EvtView->create<pxl::Particle>();
	    genmap[p] = part;	//fill genmap
            part->setName("Muon");
            part->setCharge(p->charge());
            part->setP4(p->px(), p->py(), p->pz(), p->energy());
	    part->setUserRecord<float>("Vtx_X", p->vx());
	    part->setUserRecord<float>("Vtx_Y", p->vy());
	    part->setUserRecord<float>("Vtx_Z", p->vz());
	    part->setUserRecord<int>("GenId", GenId);
	   
	    // TEMPORARY: calculate isolation ourselves 
	    //FIXME: make this at least comparable with pat/lepton isolation
	    double GenIso = IsoGenSum(iEvent, p->pt(), p->eta(), p->phi(), 0.3, 1.5);
	    part->setUserRecord<float>("GenIso", GenIso);

	    //save mother of stable muon
	    p_mother =p->mother(); 
	    if (p_mother != 0) {
	       mother = p_mother->pdgId();
	       //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	       while( abs(mother) == 13 ){
	         p_mother = p_mother->mother();
	         mother = p_mother->pdgId();
	       }	       
	       part->setUserRecord<int>("mother_id", mother);
            } else {
	       part->setUserRecord<int>("mother_id", -1);
	    }          
	    numMuonMC++; 
         }
      }

 // fill Gen Electrons passing some basic cuts
      if ( abs(p->pdgId()) == 11 && p->status() == 1) {
         if ( EleMC_cuts(p) ) { 
            pxl::Particle* part = EvtView->create<pxl::Particle>();
	    genmap[p] = part; //fill genmap
            part->setName("Ele");
            part->setCharge(p->charge());
            part->setP4(p->px(), p->py(), p->pz(), p->energy());
	    part->setUserRecord<double>("Vtx_X", p->vx());
	    part->setUserRecord<double>("Vtx_Y", p->vy());
	    part->setUserRecord<double>("Vtx_Z", p->vz());
	    part->setUserRecord<int>("GenId", GenId);
	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, p->pt(), p->eta(), p->phi(), 0.3, 1.5);
	    part->setUserRecord<double>("GenIso", GenIso);
	    //save mother of stable electron
	    p_mother =p->mother(); 
	    if (p_mother != 0) {
	       mother = p_mother->pdgId();
	       //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	       while (abs(mother) == 11) {
	          p_mother = p_mother->mother();
	          mother = p_mother->pdgId();
	       }	       
	       part->setUserRecord<int>("mother_id", mother);
            } else {
	       part->setUserRecord<int>("mother_id", -1);
	    }          
	    numEleMC++; 
         }
      }

	// fill Gen Gammas passing some basic cuts
      if ( abs(p->pdgId()) == 22 && p->status() == 1) {
         if ( GammaMC_cuts(p) ) { 
            pxl::Particle* part = EvtView->create<pxl::Particle>();
	    genmap[p] = part; //fill genmap
            part->setName("Gamma");
            part->setCharge(0);
	    part->setP4(p->px(), p->py(), p->pz(), p->energy());
            part->setUserRecord<int>("GenId", GenId);
	    // TEMPORARY: calculate isolation ourselves
	    double GenIso = IsoGenSum(iEvent, 0., p->eta(), p->phi(), 0.3, 1.5); //0. since gamma not charged!
	    part->setUserRecord<double>("GenIso", GenIso);
	    
	    //save mother of stable gamma
	    p_mother =p->mother(); 
	    if (p_mother != 0) {
	       mother = p_mother->pdgId();
	       //in case of final state radiation need to access mother of mother of mother...until particle ID really changes
	       while (abs(mother) == 22) {
	          p_mother = p_mother->mother();
	          mother = p_mother->pdgId();
	       }	       
	       part->setUserRecord<int>("mother_id", mother);
            } else {
	       part->setUserRecord<int>("mother_id", -1);
	    }          
	    numGammaMC++;
         }
      }

      GenId++;
   } //end of loop over generated-particles

   if (fDebug > 1)  cout << "MC Found:  " << numMuonMC <<  " mu  " << numEleMC << " e  " << numGammaMC << " gamma" << endl;
   EvtView->setUserRecord<int>("NumMuon", numMuonMC);
   EvtView->setUserRecord<int>("NumEle", numEleMC);
   EvtView->setUserRecord<int>("NumGamma", numGammaMC);
}

// ------------ reading the Generator Jets ------------

void ePaxAnalyzer::analyzeGenJets(const edm::Event& iEvent, pxl::EventView* EvtView) {

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
   for (reco::GenJetCollection::const_iterator genJet = ItCone5_GenJets->begin(); genJet != ItCone5_GenJets->end(); ++genJet) {
      if (JetMC_cuts(genJet)) {
         pxl::Particle* part = EvtView->create<pxl::Particle>();
         part->setName("ItCone5Jet");
         part->setP4(genJet->px(), genJet->py(), genJet->pz(), genJet->energy());
	 //fill additional jet-related infos
	 part->setUserRecord<double>("EmE", genJet->emEnergy());
	 part->setUserRecord<double>("HadE", genJet->hadEnergy());
	 part->setUserRecord<double>("InvE", genJet->invisibleEnergy());
	 part->setUserRecord<double>("AuxE", genJet->auxiliaryEnergy());
         numItCone5JetMC++;

	 //save number of GenJet-constituents fulfilling some cuts
	 numGenJetConstit_withcuts = 0; //reset variable
	 genJetConstit = genJet->getGenConstituents();
	 for (std::vector<const GenParticle*>::iterator constit = genJetConstit.begin(); constit != genJetConstit.end(); ++constit ) {
	    //raise counter if cut passed
	    if( (*constit)->pt() > constit_pT ) numGenJetConstit_withcuts++; 
	 }
	 part->setUserRecord<int>("GenJetConstit", numGenJetConstit_withcuts);
      }
   }
   EvtView->setUserRecord<int>("NumItCone5Jet", numItCone5JetMC);
   if (fDebug > 1) cout << "Found MC Jets:  "  << numItCone5JetMC << " It5  " << endl;
}

// ------------ reading the Generator MET ------------

void ePaxAnalyzer::analyzeGenMET(const edm::Event& iEvent, pxl::EventView* EvtView) {

   edm::Handle<reco::GenMETCollection> GenMet;
   iEvent.getByLabel(fMETMCLabel, GenMet);
   const GenMETCollection *genmetcol = GenMet.product();
   const GenMET genmet = genmetcol->front();  // MET exists only once!
 
   int numMETMC = 0; //means no MET in event

   pxl::Particle* part = EvtView->create<pxl::Particle>();
   part->setName("MET");
   part->setP4(genmet.px(), genmet.py(), genmet.pz(), genmet.energy());
   part->setUserRecord<double>("sumEt", genmet.sumEt());
   part->setUserRecord<double>("mEtSig", genmet.mEtSig());
   //fill additional jet-related infos
   part->setUserRecord<double>("EmE", genmet.emEnergy());
   part->setUserRecord<double>("HadE", genmet.hadEnergy());
   part->setUserRecord<double>("InvE", genmet.invisibleEnergy());
  
   if (METMC_cuts(part)) numMETMC++; 
   EvtView->setUserRecord<int>("NumMET", numMETMC);
   if (numMETMC && fDebug > 1) cout << "Event contains MET" << endl;
}

// ------------ reading the Reconstructed MET ------------
//the stored information should already contain the muon corrections plus several other corrections
//it will have certainly corrections in the future as there are certain functions for uncorrection planned

void ePaxAnalyzer::analyzeRecMET(const edm::Event& iEvent, pxl::EventView* EvtView) {

   edm::Handle<std::vector<pat::MET> > METHandle;
   iEvent.getByLabel(fMETRecoLabel, METHandle);
   std::vector<pat::MET> METs = *METHandle;
   std::vector<pat::MET>::const_iterator met = METs.begin();

   int numMETRec = 0;
   pxl::Particle* part = EvtView->create<pxl::Particle>();
   part->setName("MET");
   part->setP4(met->px(), met->py(), met->pz(), met->energy());
   part->setUserRecord<double>("sumEt", met->sumEt());
   part->setUserRecord<double>("mEtSig", met->mEtSig());
   part->setUserRecord<double>("EmEt", met->emEtFraction());         
   part->setUserRecord<double>("HadEt", met->etFractionHadronic()); 
   part->setUserRecord<double>("MaxEtEm", met->maxEtInEmTowers());   
   part->setUserRecord<double>("MaxEtHad", met->maxEtInHadTowers()); 

   if (MET_cuts(part)) numMETRec++;
   EvtView->setUserRecord<int>("NumMET", numMETRec);
}

// ------------ reading HLT and L1 Trigger Bits ------------

void ePaxAnalyzer::analyzeTrigger(const edm::Event& iEvent, pxl::EventView* EvtView) {
   /*
   //HLT trigger bits
   edm::Handle<edm::TriggerResults> hltresults;
   //HLT producer is called several times within production steps, thus need Input tag with label and process name here
   edm::InputTag hlt = edm::InputTag("TriggerResults","","HLT");
   try {iEvent.getByLabel(hlt, hltresults);} catch (...) { errMsg=errMsg + "  -- No HLTRESULTS";}
   // trigger names
   edm::TriggerNames triggerNames_;
  */

   //HLT: set to false as default
   EvtView->setUserRecord<bool>("HLT1Electron", false);
   EvtView->setUserRecord<bool>("HLT1ElectronRelaxed", false);
   EvtView->setUserRecord<bool>("HLT2Electron", false);
   EvtView->setUserRecord<bool>("HLT2ElectronRelaxed", false);
   EvtView->setUserRecord<bool>("HLT1EMHighEt", false);
   EvtView->setUserRecord<bool>("HLT1EMVeryHighEt", false);
   EvtView->setUserRecord<bool>("HLT1MuonIso", false);
   EvtView->setUserRecord<bool>("HLT1MuonNonIso", false);
   EvtView->setUserRecord<bool>("CandHLT2MuonIso", false);
   EvtView->setUserRecord<bool>("HLT2MuonNonIso", false);
   EvtView->setUserRecord<bool>("HLTNMuonNonIso", false);
   EvtView->setUserRecord<bool>("HLTXElectronMuon", false);
   EvtView->setUserRecord<bool>("HLTXElectronMuonRelaxed", false);
   EvtView->setUserRecord<bool>("HLTXElectron1Jet", false);
   EvtView->setUserRecord<bool>("HLTXElectron2Jet", false);
   EvtView->setUserRecord<bool>("HLTXElectron3Jet", false);
   EvtView->setUserRecord<bool>("HLTXElectron4Jet", false);
   EvtView->setUserRecord<bool>("HLTXMuonJets", false);
   //L1: set to false as default (have all prescale = 1)
   EvtView->setUserRecord<bool>("L1_SingleMu7", false);
   EvtView->setUserRecord<bool>("L1_DoubleMu3", false);
   EvtView->setUserRecord<bool>("L1_SingleIsoEG12", false);
   EvtView->setUserRecord<bool>("L1_SingleEG15", false);
   EvtView->setUserRecord<bool>("L1_DoubleIsoEG8", false);
   EvtView->setUserRecord<bool>("L1_DoubleEG10", false);
}

// ------------ reading Reconstructed Primary Vertices ------------

void ePaxAnalyzer::analyzeRecVertices(const edm::Event& iEvent, pxl::EventView* EvtView) {
  
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(fVertexRecoLabel, vertices);
   
   int numVertices = 0;
   
   for (reco::VertexCollection::const_iterator  vertex = vertices->begin(); vertex != vertices->end(); ++vertex ) {
      //only fill primary vertex if cuts passed
      if (Vertex_cuts(vertex)) { 
         pxl::Vertex* vtx = EvtView->create<pxl::Vertex>();
         vtx->setName("PV");
         vtx->setXYZ(vertex->x(), vertex->y(), vertex->z());
	 // chi2 of vertex-fit
         vtx->setUserRecord<double>("NormChi2", vertex->normalizedChi2() );        
         // number of tracks with origin in that vertex???
	 vtx->setUserRecord<int>("NumTracks", vertex->tracksSize());
         numVertices++;
      }
   }
   EvtView->setUserRecord<int>("NumVertices", numVertices); 
}

// ------------ reading Reconstructed Muons ------------

void ePaxAnalyzer::analyzeRecMuons(const edm::Event& iEvent, pxl::EventView* RecView, pxl::EventView* GenView, std::map<const Particle*, pxl::Particle*> & genmap) {

   // get pat::Muon's from event
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(fMuonRecoLabel, muonHandle);
   std::vector<pat::Muon> muons = *muonHandle;
   //gen particles for PAT-matching
   edm::Handle<reco::GenParticleCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );
   // count muons   
   int numMuonRec = 0;
   // loop over all pat::Muon's but store only store GLOBAL MUONs
   for (std::vector<pat::Muon>::const_iterator muon = muons.begin();  muon != muons.end(); ++muon ) {
      if (Muon_cuts(*muon)) { 
         pxl::Particle* part = RecView->create<pxl::Particle>();
         part->setName("Muon");
         part->setCharge(muon->charge());
         part->setP4(muon->px(), muon->py(), muon->pz(), muon->energy());
         part->setUserRecord<double>("Vtx_X", muon->vx());
         part->setUserRecord<double>("Vtx_Y", muon->vy());
         part->setUserRecord<double>("Vtx_Z", muon->vz()); 
         //store PAT matching info
         const reco::Particle* recogen = muon->genLepton();
         
	 pxl::Particle* pxlgen = genmap[recogen];
	 if (pxlgen != NULL) {
	    part->linkSoft(pxlgen, "pat-match");
	  
	    /*//check stored matching info
	    if (part->getSoftRelations().has(pxlgen)) {
	   	cout << "Soft-Relation ele rec -> gen ok" << endl;
	    }
	    if (pxlgen->getSoftRelations().has(part)) {
	   	cout << "Soft-Relation ele gen -> rec ok" << endl;
	    }
	    cout << "pt of the matched rec-electron: " << part->getPt() << endl;
	    cout << "pt of the matched rec-electron: " << pxlgen->getPt() << endl;
	    //end check*/
	 }
	 
	 //mind that from CMSSW 2.1 on information will be stored in outerTrack, innerTrack, globalTrack
	 //combinedMuon and standAloneMuon and track might then be deprecated sooner or later!
         // since we are only interested in Global Muons we can look at the track given by combinedMuon() or globalTrack() 
	 // both methods return the same! See DataFormats/PatCandidates/interface/Muon.h
	 
	 //save info about quality of track-fit for combined muon (muon system + tracker)
	 reco::TrackRef muontrack = muon->globalTrack();
         part->setUserRecord<double>("NormChi2", muontrack->normalizedChi2());
         part->setUserRecord<int>("ValidHits", muontrack->numberOfValidHits());
         part->setUserRecord<int>("LostHits", muontrack->numberOfLostHits());
	 //error info also used in muon-Met corrections, thus store variable to save info for later re-corrections
	 part->setUserRecord<double>("dPtRelTrack", muontrack->error(0)/(muontrack->qoverp()));
	 part->setUserRecord<double>("dPtRelTrack_off", muontrack->ptError()/muontrack->pt());
	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part->setUserRecord<double>("DSZ", muontrack->dsz());
       
	 //official CaloIso and TrkIso
	 //Def:  aMuon.setCaloIso(aMuon.isolationR03().emEt + aMuon.isolationR03().hadEt + aMuon.isolationR03().hoEt);
	 part->setUserRecord<double>("CaloIso", muon->caloIso());
	 part->setUserRecord<double>("TrkIso", muon->trackIso()); 
	 //save offical isolation information: delta R = 0.3
	 const MuonIsolation& muonIsoR03 = muon->isolationR03();
	 part->setUserRecord<double>("IsoR3SumPt", muonIsoR03.sumPt);
	 part->setUserRecord<double>("IsoR3EmEt", muonIsoR03.emEt);
	 part->setUserRecord<double>("IsoR3HadEt", muonIsoR03.hadEt);
	 part->setUserRecord<double>("IsoR3HoEt", muonIsoR03.hoEt);
	 part->setUserRecord<int>("IsoR3NTracks", muonIsoR03.nTracks);
	 part->setUserRecord<int>("IsoR3NJets", muonIsoR03.nJets);
         //save offical isolation information: delta R = 0.5
	 const MuonIsolation& muonIsoR05 = muon->isolationR05();
	 part->setUserRecord<double>("IsoR5SumPt", muonIsoR05.sumPt);
	 part->setUserRecord<double>("IsoR5EmEt", muonIsoR05.emEt);
	 part->setUserRecord<double>("IsoR5HadEt", muonIsoR05.hadEt);
	 part->setUserRecord<double>("IsoR5HoEt", muonIsoR05.hoEt);
	 part->setUserRecord<int>("IsoR5NTracks", muonIsoR05.nTracks);
	 part->setUserRecord<int>("IsoR5NJets", muonIsoR05.nJets);
	 //save some stuff related to Muon-ID (Calo-info etc.)
	 part->setUserRecord<double>("CaloComp", muon->caloCompatibility());
	 part->setUserRecord<int>("NumChambers", muon->numberOfChambers());
	 part->setUserRecord<int>("NumMatches", muon->numberOfMatches());
	 part->setUserRecord<double>("EMDeposit", muon->calEnergy().em);
	 part->setUserRecord<double>("HCALDeposit", muon->calEnergy().had);
	 // check good muon method
	 part->setUserRecord<bool>("isGood", muon->isGood(reco::Muon::GlobalMuonPromptTight));
         part->setUserRecord<float>("SegComp", muon->segmentCompatibility());
         numMuonRec++;
      }
   }
   RecView->setUserRecord<int>("NumMuon", numMuonRec);
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 
}

// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, pxl::EventView* RecView, pxl::EventView* GenView, EcalClusterLazyTools& lazyTools, std::map<const Particle*, pxl::Particle*> & genmap) {

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

   for (std::vector<pat::Electron>::const_iterator ele = electrons.begin(); ele != electrons.end(); ++ele ) {
      if (Ele_cuts(ele)) {
         if (fDebug > 1) {
	    cout << "Electron Energy scale corrected: " << ele->isEnergyScaleCorrected() 
	         << "  Momentum corrected: " << ele->isMomentumCorrected() << endl;
         }
	 pxl::Particle* part = RecView->create<pxl::Particle>();

//pxl::ParticleFilter recElectronList(recEvtView().getObjects(), kElectronName);
	//for (pxl::ParticleFilterIterator iter(recElectronList); !iter.isDone(); iter.next()) {
	//	pxl::ParticleWkPtr recPa = iter.wkPtr();
	//	int id = recPa.get().findUserRecord<int>(kCandidateId, -1);
		
		//edm::RefToBase<pat::ElectronType> eleRef = ele->originalObjectRef();
  		//const edm::Ptr<reco::Candidate> & eleRef = ele->originalObjectRef();
		//reco::GenParticle* genEle = (*genMatch)[eleRef];
		
		//if (genEle.isNonnull() && genEle.isAvailable() ) {
			//cout << "rec pt " << muoeleRef->pt() << " matched gen particle pt= " << genMuon->pt() << " key " << genMuon.key() << endl;
			//part->set().setUserRecord<int>("MatchedKey", genEle.key());
		//}
		
	//}

         part->setName("Ele");
         part->setCharge(ele->charge());
         part->setP4(ele->px(), ele->py(),ele->pz(), ele->energy());
         part->setUserRecord<double>("Vtx_X", ele->vx());
         part->setUserRecord<double>("Vtx_Y", ele->vy());
         part->setUserRecord<double>("Vtx_Z", ele->vz());
	 part->setUserRecord<float>("EoP", ele->eSuperClusterOverP()); //used for CutBasedElectronID
         part->setUserRecord<float>("HoEm", ele->hadronicOverEm()); //used for CutBasedElectronID
	 part->setUserRecord<float>("DEtaSCVtx", ele->deltaEtaSuperClusterTrackAtVtx()); //used for CutBasedElectronID
	 part->setUserRecord<float>("DPhiSCVtx", ele->deltaPhiSuperClusterTrackAtVtx()); //used for CutBasedElectronID
	 //! the seed cluster eta - track eta at calo from outermost state
	 part->setUserRecord<double>("DEtaSeedTrk", ele->deltaEtaSeedClusterTrackAtCalo()); //ok
	 part->setUserRecord<double>("DPhiSeedTrk", ele->deltaPhiSeedClusterTrackAtCalo()); //ok
	 //part->setUserRecord<double>("SCE", ele->superCluster()->energy());
         // the super cluster energy corrected by EnergyScaleFactor
	 part->setUserRecord<float>("SCE", ele->caloEnergy()); //ok
         // the errors on the supercluster energy and track momentum
         part->setUserRecord<float>("SCEErr", ele->caloEnergyError()); //ok
	 // the errors on the track momentum
         part->setUserRecord<float>("PErr", ele->trackMomentumError());	 
         part->setUserRecord<double>("TrackerP", ele->gsfTrack()->p()); //ok
	 //the seed cluster energy / track momentum at calo from outermost state
	 part->setUserRecord<double>("ESCSeedPout", ele->eSeedClusterOverPout()); //ok
         //part->setUserRecord<double>("NormChi2", ele->gsfTrack()->normalizedChi2()); //why was this left out?
         part->setUserRecord<int>("ValidHits", ele->gsfTrack()->numberOfValidHits()); //ok
	 part->setUserRecord<int>("Class", ele->classification()); //ok

	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part->setUserRecord<double>("DSZ", ele->gsfTrack()->dsz()); //ok
	 //part->setUserRecord<double>("DZ", ele->gsfTrack()->dz()); //not needed since vz()==dz()
         

	 //store PAT matching info

         const reco::Particle* recogen = ele->genLepton();

	 pxl::Particle* pxlgen = genmap[recogen];
	 if (pxlgen != NULL) {
	    part->linkSoft(pxlgen, "pat-match");
	  
	    /*//check stored matching info
	    if (part->getSoftRelations().has(pxlgen)) {
	  	cout << "Soft-Relation ele rec -> gen ok" << endl;
	    }
	    if (pxlgen->getSoftRelations().has(part)) {
	  	cout << "Soft-Relation ele gen -> rec ok" << endl;
	    }
	    cout << "pt of the matched rec-electron: " << part->getPt() << endl;
	    cout << "pt of the matched rec-electron: " << pxlgen->getPt() << endl;
	    //end check*/
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
	 part->setUserRecord<double>("e3x3",  lazyTools.e3x3(*SCSeed) );
	 part->setUserRecord<double>("e5x5",  lazyTools.e5x5(*SCSeed)  );
         std::vector<float> covariances = lazyTools.covariances(*SCSeed, 4.7 );
	 part->setUserRecord<double>("EtaEta", covariances[0] ); //used for CutBasedElectronID
	 part->setUserRecord<double>("EtaPhi", covariances[1] );
	 part->setUserRecord<double>("PhiPhi", covariances[2] );

	 //save eta/phi and DetId info from seed-cluster to prevent duplication of Electron/Photon-Candidates (in final selection)
	 part->setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part->setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 //part->setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
	 //additional data used for cutBasedElectronId in CMSSW 2_0_7
	 part->setUserRecord<double>("eSeed", SCRef->seed()->energy()); //used for CutBasedElectronID
	 part->setUserRecord<double>("pin", ele->trackMomentumAtVtx().R() ); //used for CutBasedElectronID	 
	 part->setUserRecord<double>( "pout", ele->trackMomentumOut().R() ); //used for CutBasedElectronID	
	 //store ID information
	 //Cut based ID is stored as float in 2_1_0 hence it is converted to bool for backwards compability

	 float IDfloat =  ele->leptonID("tight");
	 bool IDbool = false ;
	 if (IDfloat > 0.5) {IDbool = true;}
	 part->setUserRecord<bool>("CutBasedIDTight", IDbool);
	 IDfloat = ele->leptonID("robust");
	 IDbool = false;
	 if (IDfloat > 0.5) {IDbool = true;}	
	 part->setUserRecord<bool>("CutBasedIDRobust", IDbool);
	 
         //save official isolation information
	 double CaloIso = ele->caloIso();
	 double TrkIso = ele->trackIso();
	 part->setUserRecord<double>("CaloIso", CaloIso);
	 part->setUserRecord<double>("TrackIso", TrkIso);

         //save additional isolation information(not working for the PAT version corresponding to 2_0_9 I checked out) 
	 //not working for > 2_0_7 anyway
	 //in CMSSW_2_0_9 the data can be found in lepton.h
	 //FIXME after moving to CMSSW_2_0_9 
	 //part->setUserRecord<float>("egammaTkIso", ele->egammaTkIso());
	 //part->setUserRecord<int>("egammaTkNumIso", ele->egammaTkNumIso());
	 part->setUserRecord<double>("ECALIso", ele->ecalIso());
	 part->setUserRecord<double>("HCALIso", ele->hcalIso());
 
         //FIXME this should somehow be accessible after moving to CMSSW_2_0_9
         // Get the association vector for track number
         //edm::Handle<reco::PMGsfElectronIsoNumCollection> trackNumHandle;
         //iEvent.getByLabel(fElectronTrackNumProducer, trackNumHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 //int numVal = (*trackNumHandle)[ electronIsoRef ];
	 //part->setUserRecord<int>("TrackNum", numVal);
	 	 
         numEleRec++;
      }
      numEleAll++;
   }
   RecView->setUserRecord<int>("NumEle", numEleRec);
}

// ------------ reading Reconstructed Jets ------------
// by default for creating the allLayer0Jets the iterativeCone5CaloJets are used
// discuss if we also need informations from other Jet-Collections

void ePaxAnalyzer::analyzeRecJets(const edm::Event& iEvent, pxl::EventView* EvtView) {

   int numItCone5JetRec = 0;
   //get primary vertex (hopefully correct one) for physics eta
   double VertexZ = 0.;
   if (EvtView->findUserRecord<int>("NumVertices") > 0) {
      pxl::ObjectOwner::TypeIterator<pxl::Vertex> vtx_iter = EvtView->getObjectOwner().begin<pxl::Vertex>();
      VertexZ = (*vtx_iter)->getZ();
   } 
   
   edm::Handle<std::vector<pat::Jet> > jetHandle;
   iEvent.getByLabel(fItCone5JetRecoLabel, jetHandle);
   std::vector<pat::Jet> jets = *jetHandle;
   
   //Get the GenJet collections for PAT matching
   //edm::Handle<reco::GenJetCollection> Kt_GenJets;
   edm::Handle<reco::GenJetCollection> ItCone5_GenJets;
   //iEvent.getByLabel(fKtJetMCLabel, Kt_GenJets);
   iEvent.getByLabel(fItCone5JetMCLabel, ItCone5_GenJets);

   for (std::vector<pat::Jet>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet) {
      if (Jet_cuts(jet)) {
         pxl::Particle* part = EvtView->create<pxl::Particle>();
         part->setName("ItCone5Jet"); //compare with Clemens; this may need further discussion
         part->setP4(jet->px(), jet->py(), jet->pz(), jet->energy());
	 part->setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part->setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part->setUserRecord<int>("N90", jet->n90());
	 part->setUserRecord<int>("N60", jet->n60());
	 std::vector <CaloTowerPtr> caloRefs = jet->getCaloConstituents();
	 part->setUserRecord<int>("NCaloRefs", caloRefs.size());
	 part->setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part->setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
 	 part->setUserRecord<double>("TowersArea", jet->towersArea());
	 part->setUserRecord<double>("PhysicsEta", jet->physicsEta(VertexZ,jet->eta()));

	 /*
         //store PAT matching info
         //FIXME does not work, because implementation of genJet() differs from genLepton()
         int count = 0;
         const reco::Jet* matched = jet->genJet();
       	 for (reco::GenJetCollection::const_iterator pa = ItCone5_GenJets->begin(); pa != ItCone5_GenJets->end(); ++ pa ) {
            if (&(*pa) == matched ) {
               //cout << "MATCH!" << endl; 
               //double ptrec = jet->pt();
               //double ptgen = pa->pt(); 
               //cout << "pt des gen-jets: " << ptgen << endl;         
               //cout << "pt des rec-jets: " << ptrec<< endl;
               part->setUserRecord<double>("MatchedGenId", count); 
               count++;
            } else {
               part->setUserRecord<double>("MatchedGenId", -1);
               count++;           
            }
	 }
	*/
         numItCone5JetRec++;
      }
   }
   EvtView->setUserRecord<int>("NumItCone5Jet", numItCone5JetRec);
   if (fDebug > 1) cout << "Found Rec Jets:  " << numItCone5JetRec << " It5  " << endl;
}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventView* RecView, pxl::EventView* GenView, EcalClusterLazyTools& lazyTools, std::map<const Particle*, pxl::Particle*> & genmap){
   
   // get Photon Collection     
   edm::Handle<std::vector<pat::Photon> > photonHandle;
   iEvent.getByLabel(fGammaRecoLabel, photonHandle);
   std::vector<pat::Photon> photons = *photonHandle;
   
   int numGammaRec = 0;
   int numGammaAll = 0; //need counter for ID ???
   for (std::vector<pat::Photon>::const_iterator photon = photons.begin(); photon != photons.end(); ++photon) {  
      if ( Gamma_cuts(photon) ) { 
         pxl::Particle* part = RecView->create<pxl::Particle>();
         part->setName("Gamma");
         part->setCharge(0);
         part->setP4(photon->px(), photon->py(), photon->pz(), photon->energy());         
	 /// Whether or not the SuperCluster has a matched pixel seed
	 part->setUserRecord<bool>("HasSeed", photon->hasPixelSeed());
	 const SuperClusterRef SCRef = photon->superCluster();
         // Find the entry in the map corresponding to the seed BasicCluster of the SuperCluster
         const BasicClusterRef& SCSeed = SCRef->seed();
	 part->setUserRecord<double>("rawEnergy",  SCRef->rawEnergy() );
 	 part->setUserRecord<double>("preshowerEnergy",  SCRef->preshowerEnergy() );

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
	 part->setUserRecord<double>("e3x3",  lazyTools.e3x3(*SCSeed) );
	 part->setUserRecord<double>("e5x5",  lazyTools.e5x5(*SCSeed)  );
         std::vector<float> covariances = lazyTools.covariances(*SCSeed );
	 part->setUserRecord<double>("EtaEta", covariances[0] ); 
	 part->setUserRecord<double>("EtaPhi", covariances[1] );
	 part->setUserRecord<double>("PhiPhi", covariances[2] );
	 part->setUserRecord<double>("Emax",  lazyTools.eMax(*SCSeed)  );
	 part->setUserRecord<double>("r9", ( lazyTools.e3x3(*SCSeed) )/( SCRef->rawEnergy() + SCRef->preshowerEnergy() ) );
	 part->setUserRecord<double>("r19",  (lazyTools.eMax(*SCSeed) / lazyTools.e3x3(*SCSeed)) );

	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection) adn to reject converted photons
	 part->setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part->setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 //part->setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
	 
	 //set hadronic over electromagnetic energy fraction
	 part->setUserRecord<float>("HoEm", photon->hadronicOverEm());
 	 //save official isolation information
	 part->setUserRecord<double>("HCALIso", photon->hcalIso());
	 part->setUserRecord<double>("ECALIso", photon->ecalIso());
	 part->setUserRecord<double>("TrackIso", photon->trackIso());
	 part->setUserRecord<double>("CALIso", photon->caloIso());

	 //FIXME This part does not work. Should this not be possible in an easy way with PAT?
	 //edm::Ref<reco::PhotonCollection> gammaIsoRef(photons, numGammaAll);
	 //reco::CandidateBaseRef candRef(gammaIsoRef);
	 // Get the association vector for track number
	 //edm::Handle<reco::CandViewIntAssociations> trackNumHandle;
	 //iEvent.getByLabel(fGammaTrackNumProducer, trackNumHandle);
	 //direct access (by object). We have candidate already from HCAL-isolation
	 //int numVal = (*trackNumHandle)[ candRef ];
	 //part->setUserRecord<int>("TrackNum", numVal);

	 //store information about converted state
	 part->setUserRecord<bool>("IsConverted", photon->isConvertedPhoton());  
	 //no way to get this?
	 //part->setUserRecord<int>("ConvertedNtrk", Ntrk_conv);  

	 //FIXME store photon-Id, make sure this is what we want!
	 part->setUserRecord<float>("photonIDTight", photon->isTightPhoton());
	 // store Gamma info corrected for primary vertex (this changes direction but leaves energy of SC unchanged 
	 //get primary vertex (hopefully correct one) for physics eta THIS NEEDS TO BE CHECKED !!!
         math::XYZPoint vtx(0., 0., 0.);
//FIXME: setIndex fr Primary Vertex for fast access!
         if (RecView->findUserRecord<int>("NumVertices") > 0) {
            pxl::ObjectOwner::TypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
            vtx = math::XYZPoint((*vtx_iter)->getX(), (*vtx_iter)->getY(), (*vtx_iter)->getZ());
         } 
	 /////  Set event vertex
	 pat::Photon localPho(*photon);
	 //localPho.setVertex(vtx);  //FIXME this line does not work (missing Cluster Shape info)
	 part->setUserRecord<double>("PhysicsEta", localPho.p4().eta());
	 part->setUserRecord<double>("PhysicsPhi", localPho.p4().phi());
	 part->setUserRecord<double>("PhysicsPt", localPho.p4().pt());


//store PAT matching info
         const reco::Particle* recogen = photon->genPhoton();

	  pxl::Particle* pxlgen = genmap[recogen];
	  if(pxlgen != NULL){
	  part->linkSoft(pxlgen, "pat-match");
	  
	  /*//check stored matching info
	  if (part->getSoftRelations().has(pxlgen)) {
	  	cout << "Soft-Relation ele rec -> gen ok" << endl;
	  }
	  if (pxlgen->getSoftRelations().has(part)) {
	  	cout << "Soft-Relation ele gen -> rec ok" << endl;
	  }
	  cout << "pt of the matched rec-electron: " << part->getPt() << endl;
	  cout << "pt of the matched rec-electron: " << pxlgen->getPt() << endl;
	  //end check*/
	  }
	 
	 //FIXME Pi0 stuff still missing
         numGammaRec++;
      }	 
      numGammaAll++;
   }
   RecView->setUserRecord<int>("NumGamma", numGammaRec);
   if (fDebug > 1) cout << "Rec Gamma: " << numGammaRec << endl;
}

// ------------ method returning the EventClassType ------------

std::string ePaxAnalyzer::getEventClass(pxl::EventView* EvtView) {

   ostringstream EventType;
   //set default values to 0 for Gen-only mode
   EventType << EvtView->findUserRecord<int>("NumEle", 0) <<  "e"
             << EvtView->findUserRecord<int>("NumMuon", 0) << "mu"
             << EvtView->findUserRecord<int>("NumGamma", 0) << "gam"
             << EvtView->findUserRecord<int>("NumKtJet", 0) << "kt"
             << EvtView->findUserRecord<int>("NumMET", 0) << "met";
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

bool ePaxAnalyzer::METMC_cuts(const pxl::Particle* MCmet) const {
   // 
   if (MCmet->getPt() < 30.) return false;
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

bool ePaxAnalyzer::Muon_cuts(const pat::Muon& muon) const {
   // basic preselection cuts
   if (!muon.isGlobalMuon()) return false;
   if (muon.pt() < 15.)  return false;
   if (fabs(muon.eta()) > 3.) return false;
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

bool ePaxAnalyzer::MET_cuts(const pxl::Particle* met) const {
   // 
   if (met->getPt() < 30.) return false;
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
