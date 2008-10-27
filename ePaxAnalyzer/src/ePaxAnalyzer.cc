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
 
//includes for trigger info in 2_1_9
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"



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

   ftriggerResultsTag =iConfig.getParameter<edm::InputTag>("triggerResults");
   ftriggerEventTag = iConfig.getParameter<edm::InputTag>("triggerEvent");

   
   Matcher = new ParticleMatcher();
   fNumEvt=0;
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

   // set event counter   
   fNumEvt++; 
   // Owner of all Pxl Objects 
   pxl::Event* event = new pxl::Event();

   // event-specific data
   bool IsMC =  !iEvent.isRealData();
   event->setUserRecord<bool>("MC", IsMC);  //distinguish between MC and data
   event->setUserRecord<int>("Run", iEvent.id().run());
   event->setUserRecord<int>("ID", iEvent.id().event());	
   if (fDebug > 0) {
      cout << "Run " << iEvent.id().run() << "   EventID = " << iEvent.id().event() << " is MC = " << !iEvent.isRealData() << endl;  
   }
   // create two ePaxEventViews for Generator/Reconstructed Objects
   pxl::EventView* GenEvtView = event->createIndexed<pxl::EventView>("Gen");
   pxl::EventView* RecEvtView = event->createIndexed<pxl::EventView>("Rec");
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
      // FIXME: remove this hardcoded stuff
      const string met_name = "MET";
      const string jet_name = "ItCone5Jet";
      Matcher->matchObjects(GenEvtView, RecEvtView, jet_name, met_name);
   }

   // set event class strings
   GenEvtView->setUserRecord<std::string>("EventClass", getEventClass(GenEvtView));
   RecEvtView->setUserRecord<std::string>("EventClass", getEventClass(RecEvtView));
   
   if (fDebug > 0) {
      cout << "UserRecord  " <<  "   e   " << "  mu   " << " Gamma " << " SISC5/SISC7/IC5/KT4/KT6 " << "  MET  " << endl;
      cout << "Found (MC): " << setw(4) << GenEvtView->findUserRecord<int>("NumEle") 
           << setw(7) << GenEvtView->findUserRecord<int>("NumMuon")
           << setw(7) << GenEvtView->findUserRecord<int>("NumGamma") 
           << setw(4) << GenEvtView->findUserRecord<int>("NumSISC5Jet") << "/" 
	   << GenEvtView->findUserRecord<int>("NumSISC7Jet") << "/" 
           << GenEvtView->findUserRecord<int>("NumIC5Jet") << "/" 
	   << GenEvtView->findUserRecord<int>("NumKT4Jet") << "/" 
           << GenEvtView->findUserRecord<int>("NumKT6Jet")
           << setw(7) << GenEvtView->findUserRecord<int>("NumMET") << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView->findUserRecord<int>("NumEle")    
           << setw(7) << RecEvtView->findUserRecord<int>("NumMuon") 
           << setw(7) << RecEvtView->findUserRecord<int>("NumGamma")
           << setw(4) << RecEvtView->findUserRecord<int>("NumSISC5Jet") << "/" 
	   << RecEvtView->findUserRecord<int>("NumSISC7Jet") << "/" 
           << RecEvtView->findUserRecord<int>("NumIC5Jet") << "/" 
	   << RecEvtView->findUserRecord<int>("NumKT4Jet") << "/" 
           << RecEvtView->findUserRecord<int>("NumKT6Jet")
           << setw(7) << RecEvtView->findUserRecord<int>("NumMET") << endl;
      cout << "Gen Event Type: " << GenEvtView->findUserRecord<string>("EventClass") << endl;
      cout << "Rec Event Type: " << RecEvtView->findUserRecord<string>("EventClass") << endl;
   }   

   fePaxFile.writeEvent(event, RecEvtView->findUserRecord<string>("EventClass"));
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
  
   // read and store PDF Info
   edm::Handle<reco::PdfInfo> pdfstuff;
   iEvent.getByLabel("genEventPdfInfo", pdfstuff);
   EvtView->setUserRecord<float>("x1", pdfstuff->x1);
   EvtView->setUserRecord<float>("x2", pdfstuff->x2);
   EvtView->setUserRecord<float>("Q", pdfstuff->scalePDF);
   EvtView->setUserRecord<int>("f1", pdfstuff->id1);
   EvtView->setUserRecord<int>("f2", pdfstuff->id2);
   
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
   
/*   //TrackingVertexCollection ONLY available in RECO - not even that ... in Tauola Sample not there
   //get primary vertex from TrackingVertex. First vertex should be equal to HepMC-info, but HepMC is depreciated in AOD...
   edm::Handle<TrackingVertexCollection> TruthVertexContainer;
   //iEvent.getByLabel("mergedtruth","MergedTrackTruth", TruthVertexContainer);
   iEvent.getByType(TruthVertexContainer);

   // take only first vertex as this should correspond to generated primary vertex. Seems to be carbon-copy of HepMC-genVertex
   // Question: What about pile-up/min-bias? Will there be additional primary vertices??? How do we find them (BX-id, event-id) and should we save them?
   TrackingVertexCollection::const_iterator EventVertices = TruthVertexContainer->begin(); 
     
   // Store primary Vertex:
   pxl::Vertex* GenTrackingVtx = EvtView->create<pxl::Vertex>();
   GenTrackingVtx->setName("PTV");
   GenTrackingVtx->setXYZ(EventVertices->position().x(), EventVertices->position().y(), EventVertices->position().z());
   // we only have a single PV at Generator Level. Due to EventView Consistency .i.e. GenEvtView should look identical to RecEvtView
   // this variable is explicetly set
   EvtView->setUserRecord<int>("NumPTVertices", 1);		
   // do we need this BX/event identification???
   //GenVtx.setUserRecord<int>("Vtx_BX", EventVertices->eventId().bunchCrossing());
   //GenVtx.setUserRecord<int>("Vtx_event", EventVertices->eventId().event()); */
   
   const GenParticle* p = (const GenParticle*) &(*genParticleHandel->begin()); //this is the incoming proton
   pxl::Vertex* GenVtx = EvtView->create<pxl::Vertex>();
   GenVtx->setName("PV");
   GenVtx->setXYZ(p->daughter(0)->vx(), p->daughter(0)->vy(), p->daughter(0)->vz()); //need daughter since first particle (proton) has position zero
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
            part->setName(fJetRecoLabels[jetcoll_i]+"Jet");
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
      EvtView->setUserRecord<int>("Num"+fJetRecoLabels[jetcoll_i]+"Jet", numJetMC);
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
   */

   //reference for some of the following parts: CMSSW/HLTrigger/HLTcore/plugins/HLTEventAnalyzerAOD.cc

   //define HLT Trigger vector FIXME:there might be a better place to define this
   std::vector<string> HLTVec;
   HLTVec.reserve(20); //preventing reallocation
   HLTVec.push_back("HLT_IsoEle15_L1I");
   HLTVec.push_back("HLT_Ele15_SW_L1R");
   HLTVec.push_back("HLT_Ele15_LW_L1R");
   HLTVec.push_back("HLT_EM80");
   HLTVec.push_back("HLT_EM200");
   HLTVec.push_back("HLT_DoubleEle10_Z");
   HLTVec.push_back("HLT_IsoMu15");
   HLTVec.push_back("HLT_IsoEle15_LW_L1I");
   HLTVec.push_back("HLT_DoubleIsoEle12_L1R");
   HLTVec.push_back("HLT_DoubleEle10_Z");
   HLTVec.push_back("HLT_IsoPhoton30_L1I");
   HLTVec.push_back("HLT_IsoPhoton40_L1R");
   HLTVec.push_back("HLT_Photon25_L1R");
   HLTVec.push_back("HLT_Mu15");
   HLTVec.push_back("HLT_DoubleMu3");
   HLTVec.push_back("HLT_DoubleIsoMu3");

   using namespace edm;
   using namespace trigger;
   edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
   edm::Handle<trigger::TriggerEvent> triggerEventHandle_;

   HLTConfigProvider hltConfig_;
   hltConfig_.init("HLT");

   //cout << "Available TriggerNames are: " << endl;
   //hltConfig_.dump("Triggers"); //dump table of available HLT
 
   iEvent.getByLabel( ftriggerResultsTag, triggerResultsHandle_); 
   iEvent.getByLabel( ftriggerEventTag, triggerEventHandle_); 

   //loop over selected trigger names
   for(unsigned int i=0; i<HLTVec.size() ;++i){
	const std::string triggerName = HLTVec[i];
	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName)); //returns size() if triggerName does not exist
	const unsigned int n(hltConfig_.size()); //for checking availabilty of triggerName

   	if(triggerIndex<n){  //makes sure that triggerName is defined for the event
	
   		//save trigger path status
		EvtView->setUserRecord<bool>(triggerName+"_wasrun",triggerResultsHandle_->wasrun(triggerIndex));
		EvtView->setUserRecord<bool>(triggerName,triggerResultsHandle_->accept(triggerIndex));
		EvtView->setUserRecord<bool>(triggerName+"_error",triggerResultsHandle_->error(triggerIndex));
   		/*//begin cout of saved information for debugging
   		cout << "triggerName: " << triggerName << "triggerIndex: " << triggerIndex << endl;
   		cout << " Trigger path status:"
        		<< " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
        		<< " Accept=" << triggerResultsHandle_->accept(triggerIndex)
        		<< " Error =" << triggerResultsHandle_->error(triggerIndex)
        		<< endl;
   		//end debugging output*/
   	}
   	else { cout << "triggerName not defined for this event!" << endl;}

   
   }


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

void ePaxAnalyzer::analyzeRecMuons(const edm::Event& iEvent, pxl::EventView* RecView, const bool& MC, std::map<const Particle*, pxl::Particle*> & genmap) {


   // get pat::Muon's from event
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(fMuonRecoLabel, muonHandle);
   std::vector<pat::Muon> muons = *muonHandle;

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
	  
	    	//check stored matching info
	    	//if (part->getSoftRelations().has(pxlgen)) {
	   	//	cout << "Soft-Relation muon rec -> gen ok" << endl;
	    	//}
	    	//if (pxlgen->getSoftRelations().has(part)) {
	   	//	cout << "Soft-Relation muon gen -> rec ok" << endl;
	    	//}
	    	//	cout << "pt of the matched rec-muon: " << part->getPt() << endl;
	    	//	cout << "pt of the matched gen-muon: " << pxlgen->getPt() << endl;
	    	//end check*/
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

	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part->setUserRecord<double>("DSZ", ele->gsfTrack()->dsz()); //ok
	 //part->setUserRecord<double>("DZ", ele->gsfTrack()->dz()); //not needed since vz()==dz()
         

	 //store PAT matching info if MC
	 if (MC) {
         	const reco::Particle* recogen = ele->genLepton();
	 	pxl::Particle* pxlgen = genmap[recogen];
	 	if (pxlgen != NULL) {
	   		part->linkSoft(pxlgen, "pat-match");
	  
	    	//check stored matching info
	    	//if (part->getSoftRelations().has(pxlgen)) {
	  	//	cout << "Soft-Relation ele rec -> gen ok" << endl;
	   	//}
	   	//if (pxlgen->getSoftRelations().has(part)) {
	  	//	cout << "Soft-Relation ele gen -> rec ok" << endl;
	    	//}
	   	//cout << "pt of the matched rec-electron: " << part->getPt() << endl;
	    	//cout << "pt of the matched gen-electron: " << pxlgen->getPt() << endl;
	    //end check*/
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
	 //Cut based ID is stored as float in 2_1_0 hence it is converted to bool for backwards compatibility
	 float IDfloat =  ele->leptonID("tight");
	 bool IDbool = false ;
	 if (IDfloat > 0.5) {IDbool = true;}
	 part->setUserRecord<bool>("CutIDTight", IDbool);
	 IDfloat = ele->leptonID("robust");
	 IDbool = false;
	 if (IDfloat > 0.5) {IDbool = true;}	
	 part->setUserRecord<bool>("CutIDRobust", IDbool);

	 
         //save official isolation information
	 double CaloIso = ele->caloIso();
	 double TrkIso = ele->trackIso();
	 part->setUserRecord<double>("CaloIso", CaloIso);
	 part->setUserRecord<double>("TrkIso", TrkIso);
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
      iEvent.getByLabel("selectedLayer1Jets"+(*jet_label), jetHandle);
      std::vector<pat::Jet> RecJets = *jetHandle;
      //Get the GenJet collections for PAT matching
      edm::Handle<reco::GenJetCollection> GenJets;
      // FIXME: get generic matching Gen Collection!
      iEvent.getByLabel(fJetMCLabels[jetcoll_i], GenJets);    
      // loop over the jets
      for (std::vector<pat::Jet>::const_iterator jet = RecJets.begin(); jet != RecJets.end(); ++jet) {
         if (Jet_cuts(jet)) {
            pxl::Particle* part = RecView->create<pxl::Particle>();
            part->setName((*jet_label)+"Jet"); 
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

            //store PAT matching info if MC
	    if (MC) {
	       const reco::Particle* recogen = jet->genJet();
               //begin ugly workaround in order to deal with the genJet copies returned by genJet()
   	       for (reco::GenJetCollection::const_iterator genJetit = GenJets->begin(); genJetit != GenJets->end(); ++genJetit) {
                  if (recogen != NULL && genJetit->pt() == recogen->pt()){
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
      RecView->setUserRecord<int>("Num"+(*jet_label)+"Jet", numJetRec);
      if (fDebug > 1) cout << "Found Rec Jets:  " << numJetRec << " of Type " << *jet_label << endl;
      jetcoll_i++;
   }
}

// ------------ reading Reconstructed Gammas ------------

void ePaxAnalyzer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventView* RecView, bool& MC, EcalClusterLazyTools& lazyTools, std::map<const Particle*, pxl::Particle*> & genmap){
   
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
	 part->setUserRecord<double>("TrkIso", photon->trackIso());
	 part->setUserRecord<double>("CaloIso", photon->caloIso());

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
	  if(MC){
          const reco::Particle* recogen = photon->genPhoton();
	  pxl::Particle* pxlgen = genmap[recogen];
	  	if(pxlgen != NULL){
	  		part->linkSoft(pxlgen, "pat-match");
	  		/*//check stored matching info
	  		//if (part->getSoftRelations().has(pxlgen)) {
	  		//	cout << "Soft-Relation photon rec -> gen ok" << endl;
	  		//}
	  		//if (pxlgen->getSoftRelations().has(part)) {
	  		//	cout << "Soft-Relation photon gen -> rec ok" << endl;
	  		//}
	  		//cout << "pt of the matched rec-photon: " << part->getPt() << endl;
	  		//cout << "pt of the matched gen-photon: " << pxlgen->getPt() << endl;
	  		//end check*/
	  	}
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
   EventType << EvtView->findUserRecord<int>("NumEle") <<  "e"
             << EvtView->findUserRecord<int>("NumMuon") << "mu"
             << EvtView->findUserRecord<int>("NumGamma") << "gam"
             << EvtView->findUserRecord<int>("NumSISC5Jet") << "jet"
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
