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
#include <vector>
#include <algorithm>
// include Message Logger for Debug
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h" 
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// necessary objects:
#include "FWCore/Framework/interface/ESHandle.h"

// for electron shapes:
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"

// for ECAL enumerator
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

// for HCAL navigation used in HadOverEm calculation
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"

//for GenParticles
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"

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

//for electron-isolation
#include "DataFormats/Candidate/src/classes.h"

//for Trigger Bits
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
//ok

//For L1 and Hlt objects
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

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
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
 
//includes for trigger info in 2_1_9
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
// L1 Trigger stuff
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"

// special stuff for sim truth of converted photons
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"


using namespace std;
using namespace edm;

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
   // Use SIM info
   fUseSIM = iConfig.getUntrackedParameter<bool>("UseSIM");
   // The labels used in cfg-file 
   //fTruthVertexLabel = iConfig.getUntrackedParameter<string>("TruthVertexLabel");
   fgenParticleCandidatesLabel  = iConfig.getUntrackedParameter<string>("genParticleCandidatesLabel");
   fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fVertexRecoLabel = iConfig.getUntrackedParameter<string>("VertexRecoLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   fbarrelClusterCollection = iConfig.getParameter<edm::InputTag>("barrelClusterCollection");
   fendcapClusterCollection = iConfig.getParameter<edm::InputTag>("endcapClusterCollection");
   freducedBarrelRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
   freducedEndcapRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   // Jet labels - used for Gen and Rec
   fJetMCLabels  = iConfig.getParameter<std::vector<std::string> >("JetMCLabels");
   fJetRecoLabels = iConfig.getParameter<std::vector<std::string> >("JetRecoLabels");
   // MET label
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");

   fL1GlobalTriggerTag =  iConfig.getParameter<edm::InputTag>("L1GlobalTriggerReadoutRecord");
   fL1TriggerObjectMapTag = iConfig.getParameter<edm::InputTag>("L1TriggerObjectMapTag");
   ftriggerResultsTag = iConfig.getParameter<edm::InputTag>("triggerResults");
   fTriggerEvent = iConfig.getParameter<edm::InputTag>("triggerEvent");
   fStoreL3Objects = iConfig.getUntrackedParameter<bool>("StoreL3Objects");
   
   HLTConfigProvider hltConfig_;
   hltConfig_.init("HLT");
   //cout << "Available TriggerNames are: " << endl;
   //hltConfig_.dump("Triggers"); //dump table of available HLT

   // Electrons
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoEle15_L1I")] = "HLT_IsoEle15_L1I";
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoEle18_L1R")] = "HLT_IsoEle18_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoEle15_LW_L1I")] = "HLT_IsoEle15_LW_L1I";
   fHLTMap[hltConfig_.triggerIndex("HLT_Ele15_SW_L1R")] = "HLT_Ele15_SW_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_Ele15_LW_L1R")] = "HLT_Ele15_LW_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_EM80")] = "HLT_EM80";
   fHLTMap[hltConfig_.triggerIndex("HLT_EM200")] = "HLT_EM200";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoEle10_L1I")] = "HLT_DoubleIsoEle10_L1I";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoEle12_L1R")] = "HLT_DoubleIsoEle12_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoEle10_LW_L1I")] = "HLT_DoubleIsoEle10_LW_L1I";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoEle12_LW_L1R")] = "HLT_DoubleIsoEle12_LW_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleEle5_SW_L1R")] = "HLT_DoubleEle5_SW_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleEle10_LW_OnlyPixelM_L1R")] = "HLT_DoubleEle10_LW_OnlyPixelM_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleEle6_Exclusive")] = "HLT_DoubleEle6_Exclusive";
   
   // Photons
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoPhoton15_L1R")] = "HLT_IsoPhoton15_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoPhoton20_L1R")] = "HLT_IsoPhoton20_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_Photon15_L1R")] = "HLT_Photon15_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_Photon25_L1R")] = "HLT_Photon25_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoPhoton20_L1R")] = "HLT_DoubleIsoPhoton20_L1R";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoPhoton20_L1I")] = "HLT_DoubleIsoPhoton20_L1I";      
   fHLTMap[hltConfig_.triggerIndex("HLT_DoublePhoton10_Exclusive")] = "HLT_DoublePhoton10_Exclusive";      
   // Muons
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoMu15")] = "HLT_IsoMu15";
   fHLTMap[hltConfig_.triggerIndex("HLT_Mu15")] = "HLT_Mu15";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleIsoMu3")] = "HLT_DoubleIsoMu3";
   fHLTMap[hltConfig_.triggerIndex("HLT_DoubleMu3")] = "HLT_DoubleMu3";
   // X-Triggers
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoEle8_IsoMu7")] = "HLT_IsoEle8_IsoMu7";
   fHLTMap[hltConfig_.triggerIndex("HLT_IsoEle5_TripleJet30")] = "HLT_IsoEle5_TripleJet30";
   fHLTMap[hltConfig_.triggerIndex("HLT_Mu5_TripleJet30")] = "HLT_Mu5_TripleJet30";
   // L1 Trigger Bits:    
   fL1Map[47] = "L1_SingleMu10";
   fL1Map[51] = "L1_DoubleMu3";
   fL1Map[11] = "L1_SingleIsoEG12";
   fL1Map[17] = "L1_SingleEG10";
   fL1Map[18] = "L1_SingleEG12";	  
   fL1Map[19] = "L1_SingleEG15";
   fL1Map[82] = "L1_DoubleIsoEG8";
   fL1Map[102] = "L1_DoubleEG10";
   
   Matcher = new ParticleMatcher();
   fNumEvt=0;
}

// ------------ MIS Destructor  ------------

ePaxAnalyzer::~ePaxAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete Matcher;
}

// ------------ method called to for each event  ------------

void ePaxAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // set event counter   
   fNumEvt++; 
   // Owner of all Pxl Objects 
   pxl::Event event;

   // event-specific data
   bool IsMC =  !iEvent.isRealData();
   event.setUserRecord<bool>("MC", IsMC);  //distinguish between MC and data
   event.setUserRecord<int>("Run", iEvent.id().run());
   event.setUserRecord<int>("ID", iEvent.id().event());	
   if (fDebug > 0) {
      cout << "Run " << iEvent.id().run() << "   EventID = " << iEvent.id().event() << " is MC = " << !iEvent.isRealData() << endl;  
   }
   // create two ePaxEventViews for Generator/Reconstructed Objects
   pxl::EventView* GenEvtView = event.createIndexed<pxl::EventView>("Gen");
   pxl::EventView* RecEvtView = event.createIndexed<pxl::EventView>("Rec");
   GenEvtView->setName("Gen"); RecEvtView->setName("Rec");
   GenEvtView->setUserRecord<std::string>("Type", "Gen");
   RecEvtView->setUserRecord<std::string>("Type", "Rec");
   
   //maps for matching
   std::map<const Particle*, pxl::Particle*> genmap; 
   std::map<const Particle*, pxl::Particle*> genjetmap;

   //set process name
   GenEvtView->setUserRecord<std::string>("Process", fProcess);
   RecEvtView->setUserRecord<std::string>("Process", fProcess);

   //create object for EcalClusterLazyTools
   EcalClusterLazyTools lazyTools( iEvent, iSetup, freducedBarrelRecHitCollection, freducedEndcapRecHitCollection);
   
   // Generator stuff
   if (IsMC) {
      analyzeGenInfo(iEvent, GenEvtView, genmap);
      analyzeGenRelatedInfo(iEvent, GenEvtView);  // PDFInfo, Process ID, scale, pthat
      analyzeGenJets(iEvent, GenEvtView, genjetmap);
      analyzeGenMET(iEvent, GenEvtView);
      if (fUseSIM){
         analyzeSIM(iEvent, GenEvtView);
      }
   }
   // store Rec Objects only if requested
   if (!fGenOnly) {
      //Trigger bits
      analyzeTrigger(iEvent, RecEvtView);
      // Reconstructed stuff
      analyzeRecVertices(iEvent, RecEvtView);
      analyzeRecMuons(iEvent, RecEvtView, IsMC, genmap);
      analyzeRecElectrons(iEvent, RecEvtView, IsMC, lazyTools, genmap);
      analyzeRecJets(iEvent, RecEvtView, IsMC, genjetmap);
      analyzeRecMET(iEvent, RecEvtView);
      analyzeRecGammas(iEvent, RecEvtView, IsMC, lazyTools, genmap);
   }

   if (IsMC){
      const string met_name = "MET";
      Matcher->matchObjects(GenEvtView, RecEvtView, fJetRecoLabels, met_name);
   }

   // set event class strings
   if( IsMC ){
      GenEvtView->setUserRecord<std::string>("EventClass", getEventClass(GenEvtView));
   }
   RecEvtView->setUserRecord<std::string>("EventClass", getEventClass(RecEvtView));
   
   if (fDebug > 0) {  
      cout << "UserRecord  " <<  "   e   " << "  mu   " << " Gamma ";
      for (std::vector<std::string>::const_iterator jetAlgo = fJetRecoLabels.begin(); jetAlgo != fJetRecoLabels.end(); ++jetAlgo) {
         cout << (*jetAlgo) << " ";
      }
      cout << "  MET  " << endl;
      cout << "Found (MC): " << setw(4) << GenEvtView->findUserRecord<int>("NumEle") 
           << setw(7) << GenEvtView->findUserRecord<int>("NumMuon")
           << setw(7) << GenEvtView->findUserRecord<int>("NumGamma");
      for (std::vector<std::string>::const_iterator jetAlgo = fJetRecoLabels.begin(); jetAlgo != fJetRecoLabels.end(); ++jetAlgo) {
         cout << setw(4) << GenEvtView->findUserRecord<int>("Num"+(*jetAlgo)) << " ";
      } 
      cout << setw(7) << GenEvtView->findUserRecord<int>("NumMET") << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView->findUserRecord<int>("NumEle")    
           << setw(7) << RecEvtView->findUserRecord<int>("NumMuon") 
           << setw(7) << RecEvtView->findUserRecord<int>("NumGamma");
      for (std::vector<std::string>::const_iterator jetAlgo = fJetRecoLabels.begin(); jetAlgo != fJetRecoLabels.end(); ++jetAlgo) {
         cout << setw(4) << RecEvtView->findUserRecord<int>("Num"+(*jetAlgo)) << " ";
      } 
      cout << setw(7) << RecEvtView->findUserRecord<int>("NumMET") << endl;
      cout << "Gen Event Type: " << GenEvtView->findUserRecord<string>("EventClass") << endl;
      cout << "Rec Event Type: " << RecEvtView->findUserRecord<string>("EventClass") << endl;
   }   
   fePaxFile.writeEvent(&event, RecEvtView->findUserRecord<string>("EventClass"));
}

// ------------ reading Generator related Stuff ------------

void ePaxAnalyzer::analyzeGenRelatedInfo(const edm::Event& iEvent, pxl::EventView* EvtView) {
   // 
   // this works at least for RECO. Need to check if this works on AOD or PAT-Ntuplee 
  
   // get and store EventScale aka pt_hat 
   edm::Handle<double> genEventScale;
   iEvent.getByLabel("genEventScale", genEventScale);
   EvtView->setUserRecord<double>("pthat", *genEventScale);  // pt_hat
  
   // get and store EventWeigth
   edm::Handle<double> genEventWeight;
   iEvent.getByLabel("genEventWeight", genEventWeight);
   EvtView->setUserRecord<double>("Weight", *genEventWeight);  

//   edm::Handle<double> ktValue; 
//   iEvent.getByLabel("genEventKTValue", ktValue);
//   cout << "the KT = " << *ktValue << std::endl;

   // read and store PDF Info
   edm::Handle<reco::PdfInfo> pdfstuff;
   iEvent.getByLabel("genEventPdfInfo", pdfstuff);
   EvtView->setUserRecord<float>("x1", pdfstuff->x1);
   EvtView->setUserRecord<float>("x2", pdfstuff->x2);
   EvtView->setUserRecord<float>("Q", pdfstuff->scalePDF);
   EvtView->setUserRecord<int>("f1", pdfstuff->id1);
   EvtView->setUserRecord<int>("f2", pdfstuff->id2);
   EvtView->setUserRecord<float>("pdf1", pdfstuff->pdf1);
   EvtView->setUserRecord<float>("pdf2", pdfstuff->pdf2);
   //cout << "x1: " << pdfstuff->x1 << "   x2: " << pdfstuff->x2 << endl;
   //cout << "pdf1: " << pdfstuff->pdf1 << "   pdf2: " << pdfstuff->pdf2 << endl;
   //cout << "x1*pdf1: " << pdfstuff->x1*pdfstuff->pdf1 << "    x2*pdf2: " << pdfstuff->x2*pdfstuff->pdf2 << endl;
   // store also in the fpdf_vec 
   PDFInf pdf;
   pdf.x1 = pdfstuff->x1;
   pdf.x2 = pdfstuff->x2;
   pdf.f1 = pdfstuff->id1;
   pdf.f2 = pdfstuff->id2;
   pdf.Q = pdfstuff->scalePDF;
   fpdf_vec.push_back(pdf);
   
   if (fDebug > 0) {
      cout << "Event Scale (pthat): " << *genEventScale << ", EventWeight: " << *genEventWeight << endl;
      cout << "PDFInfo: " << endl 
           << "========" << endl;
      cout << "Momentum of first incoming parton: (id/flavour = " << static_cast<int>(pdfstuff->id1) << ") " << pdfstuff->x1 << endl
           << "Momentum of second incoming parton: (id/flavour = " << static_cast<int>(pdfstuff->id2) << ") " << pdfstuff->x2 << endl
	   << "Scale = " << pdfstuff->scalePDF << endl;
   }
}

// ------------ reading the Generator Stuff ------------

void ePaxAnalyzer::analyzeGenInfo(const edm::Event& iEvent, pxl::EventView* EvtView, std::map<const Particle*, pxl::Particle*>& genmap ) {

   //gen particles
   edm::Handle<reco::GenParticleCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );

     
   const GenParticle* p = (const GenParticle*) &(*genParticleHandel->begin()); //this is the incoming proton
   pxl::Vertex* GenVtx = EvtView->create<pxl::Vertex>();
   GenVtx->setName("PV");
   //mind that clearly the following line crashes in case of ParticleGun RelVal like single photon
   // therefore add a protection :-) C.H. 20.04.09
   if (p->daughter(0) != 0) GenVtx->setXYZ(p->daughter(0)->vx(), p->daughter(0)->vy(), p->daughter(0)->vz()); //need daughter since first particle (proton) has position zero
   else GenVtx->setXYZ(p->vx(), p->vy(), p->vz());  // if we do not have pp collisions
   EvtView->setUserRecord<int>("NumVertices", 1);  
  
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
	       while (abs(mother) == 13){
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

   } //end of loop over generated particles

   if (fDebug > 1)  cout << "MC Found:  " << numMuonMC <<  " mu  " << numEleMC << " e  " << numGammaMC << " gamma" << endl;
   EvtView->setUserRecord<int>("NumMuon", numMuonMC);
   EvtView->setUserRecord<int>("NumEle", numEleMC);
   EvtView->setUserRecord<int>("NumGamma", numGammaMC);
}

// ------------ reading the Generator Jets ------------

void ePaxAnalyzer::analyzeGenJets(const edm::Event& iEvent, pxl::EventView* EvtView, std::map<const Particle*, pxl::Particle*> & genjetmap) {

   int jetcoll_i = 0;  // counter for jet collection - needed for Rec --> Gen collection association
   // perform a loop over all desired jet collections:
   for (std::vector<std::string>::const_iterator jet_label = fJetMCLabels.begin(); jet_label != fJetMCLabels.end(); ++jet_label) {
      //Get the GenJet collections
      edm::Handle<reco::GenJetCollection> GenJets;
      iEvent.getByLabel(*jet_label, GenJets);
      // counter 
      int numJetMC = 0;
      vector <const GenParticle*> genJetConstit; //genJet-constituents
      int numGenJetConstit_withcuts = 0;
      double constit_pT = 5.; //here we have a hardcoded cut, but do we really need cfg-parameter for this?...
   
      //Loop over GenJets
      for (reco::GenJetCollection::const_iterator genJet = GenJets->begin(); genJet != GenJets->end(); ++genJet) {
         if (JetMC_cuts(genJet)) {
            pxl::Particle* part = EvtView->create<pxl::Particle>();

	    //cast iterator into GenParticleCandidate
            const GenParticle* p = (const GenParticle*) &(*genJet);
	    genjetmap[p] = part;
            part->setName(fJetRecoLabels[jetcoll_i]);
            part->setP4(genJet->px(), genJet->py(), genJet->pz(), genJet->energy());
	    //fill additional jet-related infos
	    part->setUserRecord<double>("EmE", genJet->emEnergy());
	    part->setUserRecord<double>("HadE", genJet->hadEnergy());
	    part->setUserRecord<double>("InvE", genJet->invisibleEnergy());
	    part->setUserRecord<double>("AuxE", genJet->auxiliaryEnergy());
            numJetMC++;

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
      EvtView->setUserRecord<int>("Num"+fJetRecoLabels[jetcoll_i], numJetMC);
      if (fDebug > 1) cout << "Found MC Jets:  "  << numJetMC << " of Type " << fJetRecoLabels[jetcoll_i] << endl;
      ++jetcoll_i;
   }
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


   if (fDebug > 1) cout << "GenMET before muon corr: Px = " << genmet.px() << "   Py = " << genmet.py() << "   Pt = " << part->getPt() << endl;
   // Perform Muon Corrections!
   // loop over muons and subtract them
   if (EvtView->findUserRecord<int>("NumMuon") > 0) { 
      vector<pxl::Particle*> GenMuons;
      EvtView->getObjectsOfType<pxl::Particle, pxl::PtComparator>(GenMuons, pxl::ParticleNameCriterion("Muon"));   
      for (vector<pxl::Particle*>::const_iterator muon = GenMuons.begin(); muon != GenMuons.end(); ++muon) {
         if (fDebug > 1) cout << "Correcting with " << (*muon)->getName() << " px = " << (*muon)->getPx() 
                              << " Py = " << (*muon)->getPy() << endl;
         *part -= **muon;
      }
      //reset eta-info after muon corrections
      //part->setP4(part->getPx(), part->getPy(), 0., sqrt(part->getPx()*part->getPx()+part->getPy()*part->getPy()));  
      if (fDebug > 1) cout << "GenMET after muon corr: Px = " << part->getPx() << "   Py = " << part->getPy() << "   Pt = " << part->getPt() << endl;     
   }
   if (METMC_cuts(part)) numMETMC++; 
   EvtView->setUserRecord<int>("NumMET", numMETMC);
   if (numMETMC && fDebug > 1) cout << "Event contains MET" << endl;
   //cout << " Our GenMET: " << part->getPt() << endl;
}

//----------------- SIM -------------------
void ePaxAnalyzer::analyzeSIM(const edm::Event& iEvent, pxl::EventView* EvtView) {

   Handle<SimVertexContainer> simVtcs;
   iEvent.getByLabel("g4SimHits" , simVtcs);
   SimVertexContainer::const_iterator simVertex; 

   Handle<SimTrackContainer> simTracks;
   iEvent.getByLabel("g4SimHits",simTracks);
   SimTrackContainer::const_iterator simTrack;
   SimTrackContainer::const_iterator simTrack2;

   vector<unsigned int> ParentVec;

   for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
      //int TrackID         = simTrack->trackId();
      //cout << "TrackID: " << TrackID << endl;
      int TrackType = simTrack->type();
      if ( (TrackType == 11) || (TrackType == -11) ){
         //double TrackPt = sqrt(simTrack->momentum().perp2());
         //cout << "TrackType: " << TrackType << "TrackPt: " << TrackPt << endl;
         int VtxIndex = simTrack->vertIndex();
         unsigned int ParentTrack = (*simVtcs)[VtxIndex].parentIndex();
         vector<unsigned int>::iterator where = find(ParentVec.begin(), ParentVec.end(), ParentTrack);
         if (where == ParentVec.end()){ 
            ParentVec.push_back(ParentTrack);
            //cout << "ParentTrack " << ParentTrack << endl;
            for ( simTrack2 = simTracks->begin(); simTrack2 != simTracks->end(); ++simTrack2){
               if ( simTrack2->trackId() == ParentTrack && (simTrack2->type() == 22) && (sqrt(simTrack2->momentum().perp2()) > 15.0) ){
                  //do not save photons without corresponding gen particle
                  if(!(simTrack2->noGenpart())) {                 
                     //int ParentType = simTrack2->type();
                     //double ParentPt = sqrt(simTrack2->momentum().perp2());
                     //cout << "TrackType: " << TrackType << "TrackPt: " << TrackPt << endl;
                     //cout << "ParentTrack " << ParentTrack << endl;
                     //cout << "found conversion: " << ParentType << " with pt: " << ParentPt << endl;
                     pxl::Particle* part = EvtView->create<pxl::Particle>();
                     part->setName("SIMConvGamma");
                     part->setP4(simTrack2->momentum().px(), simTrack2->momentum().py(),simTrack2->momentum().pz(), simTrack2->momentum().energy() );
                     part->setUserRecord<unsigned int>("TrackId", ParentTrack);
                     //cout << "found conversion with energy: " << simTrack2->momentum().energy() << " pt: " << part->getPt() << " eta: " << part->getEta() << " phi: " << part->getPhi() << endl;
                     //cout << "------------------" << endl;
                  }
               }
            }
         }
      }

   }

//cout << "---------NEW EVENT ---------" << endl;
}


// ------------ reading the Reconstructed MET ------------
//the stored information should already contain the muon corrections plus several other corrections
//it will have certainly corrections in the future as there are certain functions for uncorrection planned

void ePaxAnalyzer::analyzeRecMET(const edm::Event& iEvent, pxl::EventView* EvtView) {

   edm::Handle<std::vector<pat::MET> > METHandle;
   iEvent.getByLabel(fMETRecoLabel, METHandle);
   const std::vector<pat::MET>& METs = *METHandle;
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
   /*cout << "MET has " << met->nCorrections() << " corrections applied" << endl;
   cout << " Fully Corrected MET    : " << met->pt() << " (x: " << met->px() << ", y: " << met->py() << ") " << endl
        << " Uncorrected MET        : " << met->uncorrectedPt(pat::MET::uncorrALL) << "(x: " << met->corEx(pat::MET::uncorrALL) << ", y: " << met->corEy(pat::MET::uncorrALL) << ") " <<  endl
	<< " Non-JES Corrected MET  : " << met->uncorrectedPt(pat::MET::uncorrJES) << "(x: " << met->corEx(pat::MET::uncorrJES) << ", y: " << met->corEy(pat::MET::uncorrJES) << ") " << endl
	<< " Non-Muon Corrected MET : " << met->uncorrectedPt(pat::MET::uncorrMUON) << "(x: " << met->corEx(pat::MET::uncorrMUON) << ", y: " << met->corEy(pat::MET::uncorrMUON) << ") " << endl;
   */
}

// ------------ reading HLT and L1 Trigger Bits ------------

void ePaxAnalyzer::analyzeTrigger(const edm::Event& iEvent, pxl::EventView* EvtView) {

   //reference for some of the following parts: CMSSW/HLTrigger/HLTcore/plugins/HLTEventAnalyzerAOD.cc
   using namespace edm;
   using namespace trigger;
   edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
   iEvent.getByLabel(ftriggerResultsTag, triggerResultsHandle_); 
   edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
   iEvent.getByLabel(fTriggerEvent, triggerEventHandle_);
   
   HLTConfigProvider hltConfig_;
   hltConfig_.init("HLT");
   // Dump all triggers which have fired:
   //const unsigned int n(hltConfig_.size());
   //for (unsigned int i=0; i!=n; ++i) {
   //   const std::string& name = hltConfig_.triggerName(i);
   //   const unsigned int triggerIndex = hltConfig_.triggerIndex(name);
   //   if (triggerResultsHandle_->accept(triggerIndex)) cout << name << " fired ";
      //if (triggerResultsHandle_->wasrun(triggerIndex)) cout << name << " ran ";
   //}
   
   //loop over selected trigger names
   for (std::map<int, std::string>::const_iterator trig_map = fHLTMap.begin(); trig_map != fHLTMap.end(); trig_map++) {
      //save trigger path status
      if (triggerResultsHandle_->wasrun(trig_map->first) && !(triggerResultsHandle_->error(trig_map->first))) {
     	 EvtView->setUserRecord<bool>(trig_map->second, triggerResultsHandle_->accept(trig_map->first));
	 if (fDebug > 0 && triggerResultsHandle_->accept(trig_map->first)) cout << endl << "Trigger: " << trig_map->second << " fired" << endl;
      } else {
     	 if (!triggerResultsHandle_->wasrun(trig_map->first)) cout << "Trigger: " << trig_map->second << " was not executed!" << endl;
     	 if (triggerResultsHandle_->error(trig_map->first)) cout << "An error occured during execution of Trigger: " << trig_map->second << endl;
      }
      //begin cout of saved information for debugging
      if (fDebug > 1) {
     	 cout << "triggerName: " << trig_map->second << "  triggerIndex: " << trig_map->first << endl;
   	 cout << " Trigger path status:"
   	      << " WasRun=" << triggerResultsHandle_->wasrun(trig_map->first)
   	      << " Accept=" << triggerResultsHandle_->accept(trig_map->first)
   	      << " Error =" << triggerResultsHandle_->error(trig_map->first) << endl;
      }
      if (fStoreL3Objects) {
         const vector<string>& moduleLabels(hltConfig_.moduleLabels(trig_map->first));
         //const unsigned int m(hltConfig_.size(trig_map->first));
         const unsigned int moduleIndex(triggerResultsHandle_->index(trig_map->first));
         //cout << " Last active module - label/type: "
         // << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
         // << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
         // << endl;
         //assert (moduleIndex<m);

         // Results from TriggerEvent product - Attention: must look only for
         // modules actually run in this path for this event!
         // Analyze and store the objects which have fired the trigger
         for (unsigned int j = 0; j <= moduleIndex; ++j) {
            const string& moduleLabel = moduleLabels[j];
            const string  moduleType = hltConfig_.moduleType(moduleLabel);
            // check whether the module is packed up in TriggerEvent product
            const unsigned int filterIndex = triggerEventHandle_->filterIndex(InputTag(moduleLabel, "", "HLT"));
	    //cout << "FilterIndex: " << filterIndex << " for module " << moduleLabel << "/" << moduleType << endl;
            if (filterIndex < triggerEventHandle_->sizeFilters()) {
               //cout << "Trigger " << trig_map->second << ": 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
               const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
               const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
               const size_type nI(VIDS.size());
               const size_type nK(KEYS.size());
               assert(nI==nK);
               size_type n(max(nI,nK));
	       if (n > 5) {
	          cout << "Storing only 5 L3 objects for label/type " << moduleLabel << "/" << moduleType << endl;
	          n = 5;
	       }
               //cout << "   " << n  << " accepted 'L3' objects found: " << endl;
               const TriggerObjectCollection& TOC = triggerEventHandle_->getObjects();
	       for (size_type i=0; i!=n; ++i) {
	          const TriggerObject& TO = TOC[KEYS[i]];
                  pxl::Particle* part = EvtView->create<pxl::Particle>();
                  part->setName(moduleLabel);
                  part->setP4(TO.px(), TO.py(), TO.pz(), TO.energy());
                  part->setUserRecord<double>("ID", TO.id());
	          if (fDebug > 0) cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
	                               << " id: " << TO.id() << " pt: " << TO.pt() << " eta: " << TO.eta() << " phi:" << TO.phi() << " m: " << TO.mass()
	                               << endl;
               }
            }
	 }
      } 
   }
   // Store L1 Trigger Bits
   edm::Handle<L1GlobalTriggerReadoutRecord> L1GlobalTrigger;
   iEvent.getByLabel(fL1GlobalTriggerTag, L1GlobalTrigger);
   // L1 Decission
   if (!L1GlobalTrigger.failedToGet()) {
      DecisionWord gtDecisionWord = L1GlobalTrigger->decisionWord();    
      for (std::map<int,std::string>::const_iterator itMap = fL1Map.begin(); itMap != fL1Map.end(); ++itMap) {
         EvtView->setUserRecord<bool>(itMap->second, gtDecisionWord[itMap->first]);
         if (fDebug > 1) cout << "L1 TD: " << itMap->first << " " << itMap->second << " " << gtDecisionWord[itMap->first]<< endl;
      }
   } else cout << "l1GlobalTrigger Product not available!" << endl;
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
	 // errors
	 vtx->setUserRecord<double>("xErr", vertex->xError());
	 vtx->setUserRecord<double>("yErr", vertex->yError());
	 vtx->setUserRecord<double>("zErr", vertex->zError());	  
	 // chi2 of vertex-fit
         vtx->setUserRecord<double>("NormChi2", vertex->normalizedChi2());        
         // number of tracks with origin in that vertex???
	 vtx->setUserRecord<int>("NumTracks", vertex->tracksSize());
	 // is valid?
	 vtx->setUserRecord<bool>("IsValid", vertex->isValid());
	 vtx->setUserRecord<bool>("IsFake", vertex->isFake());
         numVertices++;
      }
   }
   EvtView->setUserRecord<int>("NumVertices", numVertices); 
}

// ------------ reading Reconstructed Muons ------------

void ePaxAnalyzer::analyzeRecMuons(const edm::Event& iEvent, pxl::EventView* RecView, const bool& MC, std::map<const Particle*, pxl::Particle*> & genmap) {
   // get pat::Muon's from event
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(fMuonRecoLabel, muonHandle);
   const std::vector<pat::Muon>& muons = *muonHandle;

   // count muons   
   int numMuonRec = 0;
   // loop over all pat::Muon's but only store GLOBAL MUONs
   for (std::vector<pat::Muon>::const_iterator muon = muons.begin();  muon != muons.end(); ++muon ) {
      if (Muon_cuts(*muon)) { 
         pxl::Particle* part = RecView->create<pxl::Particle>();
         part->setName("Muon");
         part->setCharge(muon->charge());
         part->setP4(muon->px(), muon->py(), muon->pz(), muon->energy());
         part->setUserRecord<double>("Vtx_X", muon->vx());
         part->setUserRecord<double>("Vtx_Y", muon->vy());
         part->setUserRecord<double>("Vtx_Z", muon->vz()); 

         //store PAT matching info if MC
	 if (MC){
         	const reco::Particle* recogen = muon->genLepton();
	 	pxl::Particle* pxlgen = genmap[recogen];
	 	if (pxlgen != NULL) {
	  		part->linkSoft(pxlgen, "pat-match");
	 	}
	 }
	 
	 //mind that from CMSSW 2.1 on information will be stored in outerTrack, innerTrack, globalTrack
	 //combinedMuon and standAloneMuon and track might then be deprecated sooner or later!
         // since we are only interested in Global Muons we can look at the track given by combinedMuon() or globalTrack() 
	 // both methods return the same! See DataFormats/PatCandidates/interface/Muon.h
	 
	 //save info about quality of track-fit for combined muon (muon system + tracker)
	 reco::TrackRef muontrack = muon->globalTrack();
	 reco::TrackRef trackerTrack = muon->innerTrack();
         part->setUserRecord<double>("NormChi2", muontrack->normalizedChi2());
         part->setUserRecord<int>("VHits", muontrack->numberOfValidHits());
         part->setUserRecord<int>("LHits", muontrack->numberOfLostHits());
	 // Tracker Only
	 part->setUserRecord<double>("NormChi2_TM", trackerTrack->normalizedChi2());
         part->setUserRecord<int>("VHits_TM", trackerTrack->numberOfValidHits());
         part->setUserRecord<int>("LHits_TM", trackerTrack->numberOfLostHits());	
	 //error info also used in muon-Met corrections, thus store variable to save info for later re-corrections
	 part->setUserRecord<double>("dPtRelTrack", muontrack->error(0)/(muontrack->qoverp()));
	 part->setUserRecord<double>("dPtRelTrack_off", muontrack->ptError()/muontrack->pt());
         // Save distance to the primary vertex in z and xy plane i.e. impact parameter
         //get primary vertex (hopefully correct one) for physics eta THIS NEEDS TO BE CHECKED !!
         // units given in cm!!! Use nominal beam spot of Summer/Fall08 if no vertex available
         math::XYZPoint vtx(0.0322, 0., 0.);
         if (RecView->findUserRecord<int>("NumVertices") > 0) {
            pxl::ObjectOwner::TypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
            vtx = math::XYZPoint((*vtx_iter)->getX(), (*vtx_iter)->getY(), (*vtx_iter)->getZ());
         } 
	 part->setUserRecord<double>("Dsz", muontrack->dsz(vtx));
         part->setUserRecord<double>("Dxy", muontrack->dxy(vtx));
       
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
         part->setUserRecord<bool>("lastStationTight", muon->isGood(reco::Muon::TMLastStationTight)); 
         part->setUserRecord<float>("SegComp", muon->segmentCompatibility());
         numMuonRec++;
      }
   }
   RecView->setUserRecord<int>("NumMuon", numMuonRec);
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 
}

// ------------ reading Reconstructed Electrons ------------

void ePaxAnalyzer::analyzeRecElectrons(const edm::Event& iEvent, pxl::EventView* RecView, bool& MC, EcalClusterLazyTools& lazyTools, std::map<const Particle*, pxl::Particle*> & genmap) {

   int numEleRec = 0;   
   int numEleAll = 0;   // for matching
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr; //wegen 2_1_0
   BasicClusterRefVector::iterator basicCluster_iterator;

   edm::Handle<std::vector<pat::Electron> > electronHandle;
   iEvent.getByLabel(fElectronRecoLabel, electronHandle);
   const std::vector<pat::Electron> &electrons = *electronHandle;

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

         // Save distance to the primary vertex in z and xy plane i.e. impact parameter
         //get primary vertex (hopefully correct one)
         // units given in cm!!! Use nominal beam spot of Summer/Fall08 if no vertex available
         math::XYZPoint vtx(0.0322, 0., 0.);
         if (RecView->findUserRecord<int>("NumVertices") > 0) {
            pxl::ObjectOwner::TypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
            vtx = math::XYZPoint((*vtx_iter)->getX(), (*vtx_iter)->getY(), (*vtx_iter)->getZ());
         } 
	 part->setUserRecord<double>("Dsz", ele->gsfTrack()->dsz(vtx));
         part->setUserRecord<double>("Dxy", ele->gsfTrack()->dxy(vtx));

	 //store PAT matching info if MC
	 if (MC) {
         	const reco::Particle* recogen = ele->genLepton();
	 	pxl::Particle* pxlgen = genmap[recogen];
	 	if (pxlgen != NULL) {
	   		part->linkSoft(pxlgen, "pat-match");
	 	}
	 }

	 // Get the supercluster (ref) of the Electron
	 // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
	 // a SuperClusterCollection is a std::vector<SuperCluster>
	 // although we get a vector of SuperClusters an electron is only made out of ONE SC
	 // therefore only the first element of the vector should be available!
	 const SuperClusterRef SCRef = ele->superCluster();
	 const BasicClusterRef& SCSeed = SCRef->seed(); 
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
	 part->setUserRecord<unsigned int>("seedId", lazyTools.getMaximum(*SCSeed).first.rawId()); 
	 //additional data used for cutBasedElectronId in CMSSW 2_0_X
	 part->setUserRecord<double>("eSeed", SCRef->seed()->energy()); //used for CutBasedElectronID
	 part->setUserRecord<double>("pin", ele->trackMomentumAtVtx().R() ); //used for CutBasedElectronID	 
	 part->setUserRecord<double>( "pout", ele->trackMomentumOut().R() ); //used for CutBasedElectronID	
	 //store ID information                    
         part->setUserRecord<bool>("RobustHighE", ele->electronID("eidRobustHighEnergy") > 0.5);
         part->setUserRecord<bool>("RobustLoose", ele->electronID("eidRobustLoose") > 0.5);
         part->setUserRecord<bool>("RobustTight", ele->electronID("eidRobustTight") > 0.5);
	 part->setUserRecord<bool>("Loose", ele->electronID("eidLoose") > 0.5);
	 part->setUserRecord<bool>("Tight", ele->electronID("eidTight") > 0.5);
	 //part->setUserRecord<float>("Likeli", ele->electronID("likelihood"));
         //save official isolation information
	 part->setUserRecord<double>("CaloIso", ele->caloIso());
	 part->setUserRecord<double>("TrkIso", ele->trackIso());
	 part->setUserRecord<double>("ECALIso", ele->ecalIso());
	 part->setUserRecord<double>("HCALIso", ele->hcalIso());
 
         //FIXME this should somehow be accessible after moving to CMSSW_2_0_9
	 // no I don't see an easy way how to do that (CH) Do we really need it?
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

void ePaxAnalyzer::analyzeRecJets(const edm::Event& iEvent, pxl::EventView* RecView, bool& MC, std::map<const Particle*, pxl::Particle*>& genjetmap) {

   int numJetRec = 0;
   //get primary vertex (hopefully correct one) for physics eta
   double VertexZ = 0.;
   if (RecView->findUserRecord<int>("NumVertices") > 0) {
      pxl::ObjectOwner::TypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
      VertexZ = (*vtx_iter)->getZ();
   } 

   int jetcoll_i = 0;  // counter for jet collection - needed for Rec --> Gen collection association
   
   // perform a loop over all desired jet collections:
   for (std::vector<std::string>::const_iterator jet_label = fJetRecoLabels.begin(); jet_label != fJetRecoLabels.end(); ++jet_label) {
      // reset counter 
      numJetRec = 0;
      // get RecoJets
      edm::Handle<std::vector<pat::Jet> > jetHandle;
      iEvent.getByLabel("cleanLayer1Jets"+(*jet_label), jetHandle);
      const std::vector<pat::Jet>& RecJets = *jetHandle;
      //Get the GenJet collections for PAT matching
      edm::Handle<reco::GenJetCollection> GenJets;
      // get matching Gen Collection
      iEvent.getByLabel(fJetMCLabels[jetcoll_i], GenJets);    
      // loop over the jets
      for (std::vector<pat::Jet>::const_iterator jet = RecJets.begin(); jet != RecJets.end(); ++jet) {
         if (Jet_cuts(jet)) {
            pxl::Particle* part = RecView->create<pxl::Particle>();
            part->setName((*jet_label)); 
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
            // store b-tag discriminator values:
	    part->setUserRecord<float>("cSVtx", jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
	    part->setUserRecord<float>("cSVtxMVA", jet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	    part->setUserRecord<float>("BProb", jet->bDiscriminator("jetBProbabilityBJetTags"));
	    part->setUserRecord<float>("Prob", jet->bDiscriminator("jetProbabilityBJetTags"));
	    // to be compared with Generator Flavor:
	    part->setUserRecord<int>("Flavour", jet->partonFlavour());
            if (fDebug > 1) {
	       const std::vector<std::pair<std::string, float> > & bTags = jet->getPairDiscri();
	       part->print(0);
	       for (vector<pair<string, float> >::const_iterator btag = bTags.begin(); btag != bTags.end(); ++btag) {
	          cout << "Name: " << btag->first << "   value: " << btag->second << endl;
	       }
	    }
            //store PAT matching info if MC
	    if (MC) {
	       const reco::Particle* recogen = jet->genJet();
               //begin ugly workaround in order to deal with the genJet copies returned by genJet()
   	       for (reco::GenJetCollection::const_iterator genJetit = GenJets->begin(); genJetit != GenJets->end(); ++genJetit) {
                  if (recogen != NULL && abs(genJetit->pt() - recogen->pt()) < 0.1 ){
	             recogen = (const GenParticle*) &(*genJetit);
                     break;
	          }
	       }
	       //end ugly workaround
               pxl::Particle* pxlgen = genjetmap[recogen];
               if (pxlgen != NULL) part->linkSoft(pxlgen, "pat-match");
	    }
            numJetRec++;
         }
      }
      RecView->setUserRecord<int>("Num"+(*jet_label), numJetRec);
      if (fDebug > 1) cout << "Found Rec Jets:  " << numJetRec << " of Type " << *jet_label << endl;
      jetcoll_i++;
   }
}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventView* RecView, bool& MC, EcalClusterLazyTools& lazyTools, std::map<const Particle*, pxl::Particle*> & genmap){
   
   // get Photon Collection     
   edm::Handle<std::vector<pat::Photon> > photonHandle;
   iEvent.getByLabel(fGammaRecoLabel, photonHandle);
   const std::vector<pat::Photon>& photons = *photonHandle;
   
   int numGammaRec = 0;
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
         //use EcalClusterLazyTools to store ClusterShapeVariables
	 part->setUserRecord<double>("e3x3",  lazyTools.e3x3(*SCSeed) );
	 part->setUserRecord<double>("e5x5",  lazyTools.e5x5(*SCSeed)  );
         std::vector<float> covariances = lazyTools.covariances(*SCSeed );
	 part->setUserRecord<double>("EtaEta", covariances[0] ); 
	 part->setUserRecord<double>("EtaPhi", covariances[1] );
	 part->setUserRecord<double>("PhiPhi", covariances[2] );
	 part->setUserRecord<double>("Emax",  lazyTools.eMax(*SCSeed)  );
	 part->setUserRecord<double>("r9", ( lazyTools.e3x3(*SCSeed) )/( SCRef->rawEnergy() + SCRef->preshowerEnergy() ) );
	 // part->setUserRecord<double>("r9", photon->r9()); <== different computation of r9 here :-(
	 part->setUserRecord<double>("r19",  (lazyTools.eMax(*SCSeed) / lazyTools.e3x3(*SCSeed)) );
	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection) adn to reject converted photons
	 part->setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part->setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 part->setUserRecord<unsigned int>("seedId", lazyTools.getMaximum(*SCSeed).first.rawId()); 
	 //set hadronic over electromagnetic energy fraction
	 part->setUserRecord<float>("HoEm", photon->hadronicOverEm());
 	 //save official isolation information
	 // this is the BAD PAT isolation!!!
	 part->setUserRecord<float>("HCALIso", photon->hcalIso());
	 part->setUserRecord<float>("ECALIso", photon->ecalIso());
	 part->setUserRecord<float>("TrkIso", photon->trackIso());
	 part->setUserRecord<float>("CaloIso", photon->caloIso());
	 part->setUserRecord<int>("TrackNum", photon->nTrkSolidCone());
	 // use egamma isolation based on RecHits:
	 part->setUserRecord<float>("ID_HCALIso", photon->isolationHcalRecHit());
	 part->setUserRecord<float>("ID_ECALIso", photon->isolationEcalRecHit());
	 part->setUserRecord<float>("ID_TrkIso", photon->isolationHollowTrkCone());	  
	 //store information about converted state
	 part->setUserRecord<bool>("Converted", photon->isConverted());
	 part->setUserRecord<bool>("AlsoEle", photon->isAlsoElectron()); 
	 //if (photon->isConvertedPhoton() == true) {cout << "is Converted!" << endl;}
	 // store photon-Id
	 part->setUserRecord<bool>("LooseEM", photon->isLooseEM());
	 part->setUserRecord<bool>("Loose", photon->isLoosePhoton());
	 part->setUserRecord<bool>("Tight", photon->isTightPhoton());
	 // is near a gap ! isEEGap == always false not yet implemented ... CMSSW_2_1_9 
	 part->setUserRecord<bool>("Gap", photon->isEBGap() || photon->isEEGap() || photon->isEBEEGap());
	 // store Gamma info corrected for primary vertex (this changes direction but leaves energy of SC unchanged 
	 //get primary vertex (hopefully correct one) for physics eta THIS NEEDS TO BE CHECKED !!!
         math::XYZPoint vtx(0.0322, 0., 0.);
         if (RecView->findUserRecord<int>("NumVertices") > 0) {
            pxl::ObjectOwner::TypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
            vtx = math::XYZPoint((*vtx_iter)->getX(), (*vtx_iter)->getY(), (*vtx_iter)->getZ());
         } 
	 /////  Set event vertex
	 pat::Photon localPho(*photon);
	 localPho.setVertex(vtx);  //FIXME this line does not work (missing Cluster Shape info)
	 part->setUserRecord<double>("PhysEta", localPho.p4().eta());
	 part->setUserRecord<double>("PhysPhi", localPho.p4().phi());
	 part->setUserRecord<double>("PhysPt", localPho.p4().pt());
         //store PAT matching info
         if (MC) {
            const reco::Particle* recogen = photon->genPhoton();
            pxl::Particle* pxlgen = genmap[recogen];
            if (pxlgen != NULL) part->linkSoft(pxlgen, "pat-match");
         }
	 //FIXME Pi0 stuff still missing --> seems not be working in CMSSW_2_1_X
         numGammaRec++;
      }	 
   }
   RecView->setUserRecord<int>("NumGamma", numGammaRec);
   if (fDebug > 1) cout << "Rec Gamma: " << numGammaRec << endl;
}

// ------------ method returning the EventClassType ------------

std::string ePaxAnalyzer::getEventClass(pxl::EventView* EvtView) {

   ostringstream EventType;
   //set default values to 0 for Gen-only mode
   EventType << EvtView->findUserRecord<int>("NumEle") <<  "e"
             << EvtView->findUserRecord<int>("NumMuon") << "mu"
             << EvtView->findUserRecord<int>("NumGamma") << "gam"
             << EvtView->findUserRecord<int>("Num"+fJetRecoLabels[0]) << "jet"
             << EvtView->findUserRecord<int>("NumMET") << "met";
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
   // close output file:
   fePaxFile.close();

   // evaluate PDF Info
   if (fpdf_vec.size() > 0) {
      vector<float> best_fit;
      vector<vector<float> > weights;
      // create for each event an empty vector
      for (unsigned int i = 0; i < fpdf_vec.size(); ++i) weights.push_back(vector<float>());
      const char *lhaPDFPath = getenv("LHAPATH");
      string pdfSet(lhaPDFPath);
      string::size_type loc = pdfSet.find( ":", 0 );
      if (loc != string::npos) pdfSet = pdfSet.substr(0,loc); 
      pdfSet.append("/cteq61.LHgrid");
      cout << "PDF set - " << pdfSet.data() << endl;
      initpdfset_((char *)pdfSet.data(), pdfSet.size());
      // loop over all subpdf's
      for (int subpdf = 0; subpdf < 41; subpdf++) {
         initpdf_(subpdf);
	 //cout << "Initialized sub PDF " << subpdf << endl;
         if (subpdf == 0) {
	    // loop over all PDFInf's
	    for (vector<PDFInf>::const_iterator pdf = fpdf_vec.begin(); pdf != fpdf_vec.end(); ++pdf) {
 	       best_fit.push_back(xfx(pdf->x1, pdf->Q, pdf->f1) * xfx(pdf->x2, pdf->Q, pdf->f2));
               //cout << "xfx1: " <<  xfx(pdf->x1, pdf->Q, pdf->f1) << "   xfx1: " << xfx(pdf->x2, pdf->Q, pdf->f2) << endl;
	    }
	 } else {
	    vector<float>::const_iterator best_fit_iter = best_fit.begin();
	    vector<vector<float> >::iterator weights_iter = weights.begin();
	    // loop over all PDFInf's
	    for (vector<PDFInf>::const_iterator pdf = fpdf_vec.begin(); pdf != fpdf_vec.end(); ++pdf) {
	       weights_iter->push_back(xfx( pdf->x1, pdf->Q, pdf->f1) * xfx( pdf->x2, pdf->Q, pdf->f2) / (*best_fit_iter));
	       ++weights_iter;
	       ++best_fit_iter;
	    }
	 }
      }
      // ReRead the pxlio file and store PDFInfo
      pxl::InputFile Input(fFileName);
      pxl::OutputFile tmpFile("Tmp"+fFileName);
      vector<vector<float> >::const_iterator weights_iter = weights.begin();
      int count = 1;
      // run event loop:
      while (Input.next()) {
         pxl::Event event;
         // read event from disk
         Input.readEvent(&event);
         // get all stored EventViews
         pxl::EventView* GenEvtView = event.getObjectOwner().findObject<pxl::EventView>("Gen");
         pxl::EventView* RecEvtView = event.getObjectOwner().findObject<pxl::EventView>("Rec");
	 string EC_string = RecEvtView->findUserRecord<std::string>("EventClass");
	 unsigned int i = 1;
	 for (vector<float>::const_iterator weight = (*weights_iter).begin(); weight != (*weights_iter).end(); ++weight) {
            //cout << "weight w" << i << "  " << *weight << endl;
            ostringstream aStream;
            aStream << "w" << i;
            string str_i = aStream.str();
            GenEvtView->setUserRecord<float>(str_i, *weight);
            RecEvtView->setUserRecord<float>(str_i, *weight);
            i++;
         }
	 //string EC_string = RecEvtView->findUserRecord<std::string>("EventClass");
         tmpFile.writeEvent(&event, EC_string);
	 ++weights_iter;
	 ++count;
      }
      tmpFile.close();
      Input.close();
      // rename tmporary file
      system(("mv Tmp" + fFileName + " " + fFileName).c_str());     
   }
}


// ------------ method to define MC-MUON-cuts

bool ePaxAnalyzer::MuonMC_cuts(const GenParticle* MCmuon) const {
   //
   if (MCmuon->pt() < 10.) return false;
   if (fabs(MCmuon->eta()) > 3.) return false;
   return true;
}
 


// ------------ method to define MC-Electron-cuts

bool ePaxAnalyzer::EleMC_cuts(const GenParticle* MCele) const {
   //
   if (MCele->pt() < 10.) return false;
   if (fabs(MCele->eta()) > 3.) return false;
   return true;
}

// ------------ method to define MC-Gamma-cuts

bool ePaxAnalyzer::GammaMC_cuts(const GenParticle* MCgamma) const {
   //
   if (MCgamma->pt() < 10.) return false;
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
  double rV = sqrt( (0.0322-vertex->x())*(0.0322-vertex->x()) + vertex->y() * vertex->y() );
  if (fabs(zV)>20. || rV>1. ) return false;
  return true;
}

// ------------ method to define MUON-cuts

bool ePaxAnalyzer::Muon_cuts(const pat::Muon& muon) const {
   // basic preselection cuts
   if (!muon.isGlobalMuon()) return false;
   if (muon.pt() < 10.)  return false;
   if (fabs(muon.eta()) > 3.) return false;
   return true;
}


// ------------ method to define ELECTRON-cuts

bool ePaxAnalyzer::Ele_cuts(std::vector<pat::Electron>::const_iterator ele) const {
   if (ele->pt() < 10.) return false;
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
   if (photon->pt() < 10.) return false;
   if (fabs(photon->eta()) > 3.) return false;
   return true;
}


// ------------ method to define MET-cuts

bool ePaxAnalyzer::MET_cuts(const pxl::Particle* met) const {
   // 
   if (met->getPt() < 30.) return false;
   return true;
}

//------------------------------------------------------------------------------

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
