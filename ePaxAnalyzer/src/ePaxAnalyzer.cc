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
#include "DataFormats/EgammaCandidates/interface/ConvertedPhoton.h"

//for electron-isolation
#include "DataFormats/Candidate/src/classes.h"
#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoNumCollection.h"
#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoCollection.h"

//for Trigger Bits
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"

//For L1 and Hlt objects
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//math stuff from Physics tools
#include "DataFormats/Math/interface/deltaR.h"

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
   fIsCSASoup = iConfig.getUntrackedParameter<bool>("IsCSASoup");
   fTruthVertexLabel = iConfig.getUntrackedParameter<string>("TruthVertexLabel");
   fgenParticleCandidatesLabel  = iConfig.getUntrackedParameter<string>("genParticleCandidatesLabel");
   fKtJetMCLabel = iConfig.getUntrackedParameter<string>("KtJetMCLabel");
   fItCone5JetMCLabel = iConfig.getUntrackedParameter<string>("ItCone5JetMCLabel");  
   fMETMCLabel = iConfig.getUntrackedParameter<string>("METMCLabel");
   fVertexRecoLabel = iConfig.getUntrackedParameter<string>("VertexRecoLabel");
   fSAMuonRecoLabel = iConfig.getUntrackedParameter<string>("SAMuonRecoLabel");
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   fBarrelClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("barrelClusterShapeAssociation");
   fEndcapClusterShapeAssocProducer = iConfig.getParameter<edm::InputTag>("endcapClusterShapeAssociation");   
   fPixelMatchElectronRecoLabel = iConfig.getUntrackedParameter<string>("PixelMatchElectronRecoLabel");
   fElectronIDAssocProducer = iConfig.getParameter<std::string>("ElectronIDAssocProducer");
   fElectronHcalIsolationProducer = iConfig.getParameter<std::string>("ElectronHcalIsolationProducer");
   fElectronEcalIsolationProducer = iConfig.getParameter<std::string>("ElectronEcalIsolationProducer");
   fElectronTrackIsolationProducer = iConfig.getParameter<std::string>("ElectronTrackIsolationProducer");
   fElectronTrackNumProducer = iConfig.getParameter<std::string>("ElectronTrackNumProducer");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   fGammaHcalIsolationProducer = iConfig.getParameter<std::string>("GammaHcalIsolationProducer");
   fGammaEcalIsolationProducer = iConfig.getParameter<std::string>("GammaEcalIsolationProducer");
   fGammaTrackIsolationProducer = iConfig.getParameter<std::string>("GammaTrackIsolationProducer");
   fGammaTrackNumProducer = iConfig.getParameter<std::string>("GammaTrackNumProducer");
   fKtJetRecoLabel = iConfig.getUntrackedParameter<string>("KtJetRecoLabel");
   fItCone5JetRecoLabel = iConfig.getUntrackedParameter<string>("ItCone5JetRecoLabel");
   fL2L3JESic5JetRecoLabel = iConfig.getUntrackedParameter<string>("L2L3JESic5JetRecoLabel");
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");
   fMETCorrRecoLabel = iConfig.getUntrackedParameter<string>("METCorrRecoLabel");
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
   
   // Store Run and Event ID
   GenEvtView.set().setUserRecord<int>("Run", iEvent.id().run());
   GenEvtView.set().setUserRecord<int>("ID", iEvent.id().event());
   RecEvtView.set().setUserRecord<int>("Run", iEvent.id().run());
   RecEvtView.set().setUserRecord<int>("ID", iEvent.id().event());
   //cout << "Run " << iEvent.id().run() << "   EventID = " << iEvent.id().event() << endl;

   // FIXME: Is this the only thing to be changed to run over the soup?!? 
   // The process is not defined by the cfg file but by the ID given by the csaweightproducer?
   // set physics process 
   if (!fIsCSASoup) {  
      GenEvtView.set().setUserRecord<string>("Process", fProcess);
      RecEvtView.set().setUserRecord<string>("Process", fProcess);
   }

   double processID = 0.;
   double pthat = 0;

   // if running on RECO use HepMC    
   try {
      edm::Handle<edm::HepMCProduct> MC;
      iEvent.getByLabel( "source", MC );
      const HepMC::GenEvent * genEvt = MC->GetEvent();
    
      processID = genEvt->signal_process_id();
      pthat = genEvt->event_scale(); 

      //genEvt->print();

   } catch (...) {   } 
   
   // else use stuff inside AOD
   try {
      edm::Handle< int > genProcessID;
      iEvent.getByLabel( "genEventProcID", genProcessID );
      processID = *genProcessID;
 
      edm::Handle< double > genEventScale;
      iEvent.getByLabel( "genEventScale", genEventScale );
      pthat = *genEventScale;
    } catch (...) { } 
    
   // NOT STORED!
   //HepMC::PdfInfo* pdfstuff = genEvt->pdf_info();
   //if (pdfstuff != 0) cout << "Momentum of first incoming parton:  (id/flavour = " << pdfstuff->id1() << ")  " <<  pdfstuff->x1() << endl
   //     << "Momentum of second incoming parton: (id/flavour = " << pdfstuff->id2() << ")  " <<  pdfstuff->x2() << endl
   //     << "Scale = " << pdfstuff->scalePDF() << endl;
   //else cout << "PDFstuff not set!" << endl;
   
   // store the ID in both event views? 
   GenEvtView.set().setUserRecord<double>("pthat", pthat);
   RecEvtView.set().setUserRecord<double>("pthat", pthat); 

   // get Event Weight in case of a soup else eventweight is == 1
   double weight = 1.;
   if (fIsCSASoup) {
      //cout << "This is a CSA07 Soup" << endl;
      edm::Handle< double> weightHandle;
      iEvent.getByLabel ("csaweightproducer","weight", weightHandle);
      weight = *weightHandle;
      // Store Process ID 
      edm::Handle< int > genProcessID;
      // Alpen Fix
      if (processID == 4) { // 4 in Pythia means external process
         iEvent.getByLabel ("csaweightproducer","AlpgenProcessID", genProcessID);
         processID = *genProcessID;
      }
   }

   //cout << "Weight = " << weight << endl;
   //cout << "Process ID: " << processID << " and Event Scale (pthat): " << pthat << endl;

   GenEvtView.set().setUserRecord<double>("Weight", weight);
   RecEvtView.set().setUserRecord<double>("Weight", weight);
   
   GenEvtView.set().setUserRecord<double>("ProcID", processID);
   RecEvtView.set().setUserRecord<double>("ProcID", processID);

   // Generator stuff
   analyzeGenInfo(iEvent, GenEvtView);
   analyzeGenJets(iEvent, GenEvtView);
   if( fGenOnly == false ){ //only if info is in event
     analyzeGenMET(iEvent, GenEvtView);
     //Trigger bits
     analyzeTrigger(iEvent, RecEvtView);
     // Reconstructed stuff
     analyzeRecVertices(iEvent, RecEvtView);
     analyzeRecMuons(iEvent, RecEvtView);
     analyzeRecElectrons(iEvent, RecEvtView);
     analyzeRecJets(iEvent, RecEvtView);
     analyzeRecMET(iEvent, RecEvtView);
     analyzeRecGammas(iEvent, RecEvtView);
   }

   Matcher->matchObjects(GenEvtView, RecEvtView);

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

// ------------ tracing status 3 particles ------------

void ePaxAnalyzer::catchParticlesWithStatus3Daughters(std::vector<const reco::Candidate*>& cand, const reco::Candidate* p) {

   //bool hasStatus3DaughtersOnly = true;
   for (unsigned int i = 0; i < p->numberOfDaughters(); i++) {
      const Candidate* daughter = p->daughter(i); 
      if (daughter->status() == 3) catchParticlesWithStatus3Daughters(cand, daughter);
      else return;
   }
   cand.push_back(p);
}

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
     // find protons initiating the hard process:
/*     if ((p->pdgId() == 2212) && (fabs(p->pz()) > 6.999e3)) {
        numOfProton++;
        //cout << "Found Proton: P (" << p->px() << ", " << p->py() << ", " << p->pz() << ") status: " << p->status() << " #daugthers: " << p->numberOfDaughters() << endl;
	// Get daughters:
	unsigned int numberOfPartons = 0;
	bool checkStatus = false;
	if ( p->numberOfDaughters() > 1) checkStatus = true;
	for (unsigned int i = 0; i < p->numberOfDaughters(); i++) {
	   const Candidate* daughter = p->daughter(i); 
	   // check if daughter is a quark or a gluon
	   if (abs(daughter->pdgId()) < 7 ||  daughter->pdgId() == 21 ) {
	      // This is needed for Pythia: only in Pythia the proton might have several quark/gluon daughters,
	      // THe code runs with MC@NLO and Pythia now. What is with e.g. Alpgen?
	      // FIXME / CHECKME !!!
	      if (checkStatus && daughter->status() != 3) continue;    
	      numberOfPartons++;
	      if (numOfProton < 3) {
	         // Set X1 (momentum fraction) and Flavour1 
		 EvtView.set().setUserRecord<double>("X"+numOfProton, daughter->pz());
		 if (daughter->pdgId() == 21) EvtView.set().setUserRecord<int>("F"+numOfProton, 0);
		 else EvtView.set().setUserRecord<int>("F"+numOfProton, daughter->pdgId());
	      } else cout <<  " aAHHHHHHHHHHH Found more then 2 Protons with energy large 6.999 TeV)" << endl;
	      //cout << "Found Proton Daughter: " << daughter->pdgId() << " with P (" << daughter->px() << ", " << daughter->py() << ", " << daughter->pz() << ") status: " << daughter->status() << " #daugthers: " << daughter->numberOfDaughters() << endl;
           } //else {
	   //   cout << "daughter is neither a gluon nor a quark: ID: " << daughter->pdgId() << endl;
	   //}  
	}
        if (numberOfPartons != 1) cout << "Found in total " << numberOfPartons << " Partons with Status code 3 which arise from Protons" << endl; 
	assert(numberOfPartons == 1); 
	
        // Now try to find out the hard process:
	// Idea: find all Particles which have status code 3 and all their daughters have status code 3
	// start from both protons and try to find the overlap latter on
	if (numOfProton == 1) catchParticlesWithStatus3Daughters(cand_prot1, p);
	if (numOfProton == 2) catchParticlesWithStatus3Daughters(cand_prot2, p);

     } */
  
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
   } //end of loop over generated-particles

/*   
   // Analyze the Status 3 Particles wiht Status 3 Only daughter
   cout << endl << "   Proton 1: Stat 3 Particles with Stat 3 daughters: " << endl;
   for (vector<const Candidate*>::const_iterator iter = cand_prot1.begin(); iter != cand_prot1.end(); iter++) {
      cout << "   Particle " << (*iter)->pdgId() << "  (" << (*iter)->px() << ", " << (*iter)->py() << ", " << (*iter)->pz() << ") status: " << (*iter)->status() << " #daugthers: " << (*iter)->numberOfDaughters() << endl;
   }
   cout << endl << "   Proton 2: Stat 3 Particles with Stat 3 daughters: " << endl;
   for (vector<const Candidate*>::const_iterator iter = cand_prot2.begin(); iter != cand_prot2.end(); iter++) {
      cout << "   Particle " << (*iter)->pdgId() << "  (" << (*iter)->px() << ", " << (*iter)->py() << ", " << (*iter)->pz() << ") status: " << (*iter)->status() << " #daugthers: " << (*iter)->numberOfDaughters() << endl;
   }

   bool found = false;
   const Candidate* cand1;
   const Candidate* cand2;
   for (vector<const Candidate*>::const_iterator iter1 = cand_prot1.begin(); iter1 != cand_prot1.end(); iter1++) {
      // check if first daughter is identical with one of the daughters of the other vector
      const Candidate* daughter1 = (*iter1)->daughter(0);
      for (vector<const Candidate*>::const_iterator iter2 = cand_prot2.begin(); iter2 != cand_prot2.end(); iter2++) {
         // iterate over all daughters and check if equal to daughter1
         for (unsigned int i = 0; i < (*iter2)->numberOfDaughters(); i ++) {
	    if (daughter1 == (*iter2)->daughter(i)) {
	       found = true;
	       cand2 = (*iter2);
	       break;
	    }
	 }
	 if (found) break;
      }
      if (found) {
       	 cand1 = (*iter1);
	 break;
      }
   }
   if (found) {
      cout << endl << "Cand1: " << cand1->pdgId() << "  (" << cand1->px() << ", " << cand1->py() << ", " << cand1->pz() << ") status: " << cand1->status() << " #daugthers: " << cand1->numberOfDaughters() << endl;
      cout << "Cand2: " << cand2->pdgId() << "  (" << cand2->px() << ", " << cand2->py() << ", " << cand2->pz() << ") status: " << cand2->status() << " #daugthers: " << cand2->numberOfDaughters() << endl;
      if (cand1->numberOfDaughters() == 1) {
         // Q scale is mass of daughter
	 cout <<  "   Q = " << cand1->daughter(0)->mass() << endl;
      } else if (cand1->numberOfDaughters() == 2) {
         const Candidate* daug1 = cand1->daughter(0);
         const Candidate* daug2 = cand2->daughter(1);
         TLorentzVector v3(daug1->px(), daug1->py(), daug1->pz(), daug1->energy());
         TLorentzVector v4(daug2->px(), daug2->py(), daug2->pz(), daug2->energy());
         TVector3 boost = -(v3 + v4).BoostVector();
         v3.Boost(boost);
         double Q_square = v3.Pt()*v3.Pt() + (v3.M2() + v4.M2())/2;
         cout << "    Q^2 = " << Q_square << "   using pthat = " << v3.Pt() << endl;
      } else {
         cout << "AAAAAAAAHHHHHH you should never end up here!" << endl;
      } 

   } else {
      cout << "No matching" << endl;
   }
   cout << endl; */
   
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
   iEvent.getByLabel(fKtJetMCLabel, Kt_GenJets);
   iEvent.getByLabel(fItCone5JetMCLabel, ItCone5_GenJets);

   int numKtJetMC = 0;
   int numItCone5JetMC = 0;
 
   vector <const GenParticleCandidate*> genJetConstit; //genJet-constituents
   int numGenJetConstit_withcuts = 0;
   double constit_pT = 5.; //here we have a hardcoded cut, but do we really need cfg-parameter for this?...
   
   //Loop over KtGenJets
   for( reco::GenJetCollection::const_iterator genJet = Kt_GenJets->begin(); 
         genJet != Kt_GenJets->end(); ++genJet ) {
      if (JetMC_cuts(genJet)) { 
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("KtJet");
         part.set().vector(pxl::set).setPx(genJet->px());
         part.set().vector(pxl::set).setPy(genJet->py());
         part.set().vector(pxl::set).setPz(genJet->pz());
         part.set().vector(pxl::set).setE((genJet->energy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
	 part.set().setUserRecord<double>("AuxE", genJet->auxiliaryEnergy());
	 numKtJetMC++;
	
	 //save number of GenJet-constituents fullfilling sone cuts
	 numGenJetConstit_withcuts = 0; //reset variable
	 genJetConstit = genJet->getConstituents();
	 for( std::vector<const GenParticleCandidate*>::iterator constit = genJetConstit.begin(); 
	      constit != genJetConstit.end(); ++constit ) {
	    
	   //raise counter if cut passed
	   if( (*constit)->pt()>constit_pT ) numGenJetConstit_withcuts++; 
	
	 }//end-of loop over all constituents
	
	 part.set().setUserRecord<int>("GenJetConstit", numGenJetConstit_withcuts);

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
         part.set().vector(pxl::set).setE((genJet->energy()));
	 //fill additional jet-related infos
	 part.set().setUserRecord<double>("EmE", genJet->emEnergy());
	 part.set().setUserRecord<double>("HadE", genJet->hadEnergy());
	 part.set().setUserRecord<double>("InvE", genJet->invisibleEnergy());
	 part.set().setUserRecord<double>("AuxE", genJet->auxiliaryEnergy());
         numItCone5JetMC++;

	 //save number of GenJet-constituents fullfilling sone cuts
	 numGenJetConstit_withcuts = 0; //reset variable
	 genJetConstit = genJet->getConstituents();
	 for( std::vector<const GenParticleCandidate*>::iterator constit = genJetConstit.begin(); 
	      constit != genJetConstit.end(); ++constit ) {
	    
	   //raise counter if cut passed
	   if( (*constit)->pt()>constit_pT ) numGenJetConstit_withcuts++; 
	
	 }//end-of loop over all constituents
	
	 part.set().setUserRecord<int>("GenJetConstit", numGenJetConstit_withcuts);

	 //save also dummy L2L3JESic5Jet for GenJet in order to be consistent with Rec-jet --> we need symmetry in naming between gen and rec
	 pxl::ParticleRef partcorr = EvtView.set().create<pxl::Particle>(part.get());
	 partcorr.set().setName("L2L3JESic5Jet");

      }

   }
   EvtView.set().setUserRecord<int>("NumItCone5Jet", numItCone5JetMC);
   //save also dummy L2L3JESic5Jet for GenJet in order to be consistent with Rec-jet --> we need symmetry in naming between gen and rec
   EvtView.set().setUserRecord<int>("NumL2L3JESic5Jet", numItCone5JetMC);
  
   if (fDebug > 1) cout << "Found MC Jets:  " << numKtJetMC << " Kt  " << numItCone5JetMC << " It5  " << endl;
   
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

   //no Muon-corrections needed if using MET-collection "genMetNoNu" since here muon already included (mu NOT in exclude list)! -> stuff commented out

   /*if (fDebug > 1) cout << "GenMET before muon corr: Px = " << genmet.px() << "   Py = " << genmet.py() << "   Pt = " << part.get().vector().getPt() << endl;
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
   part.set().vector(pxl::set).setMass(0.);   */

   //there is always MET in event, just decide if cuts passed (do this after muon corrections!)
   if (METMC_cuts(part)) { 
     numMETMC++; 
   }
   EvtView.set().setUserRecord<int>("NumMET", numMETMC);
   if (numMETMC && fDebug > 1) cout << "Event contains MET" << endl;

   //save also dummy METCorr for GenMET in order to be consistent with Rec-MET --> we need symmetry in naming between gen and rec
   //pxl::Particle partcorr(part);
   pxl::ParticleRef partcorr = EvtView.set().create<pxl::Particle>(part.get());
   partcorr.set().setName("METCorr");
   EvtView.set().setUserRecord<int>("NumMETCorr", numMETMC);
}



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

  if (&hltresults) {
    int ntrigs=hltresults->size();
    if (ntrigs==0){std::cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << std::endl;}
    
    triggerNames_.init(*hltresults);
   
    // ...Fill the corresponding accepts in branch-variables, save also selected HLT-objects (need match with gen/reco-objects for trigger-efficiencies)
    for (int itrig = 0; itrig != ntrigs; ++itrig){
      
      string trigName=triggerNames_.triggerName(itrig);
      bool accept = hltresults->accept(itrig);
      string hltobject;
      
      if (fDebug > 1) { 
	std::cout << "%HLTInfo --  Number of HLT Triggers: " << ntrigs << std::endl;
	std::cout << "%HLTInfo --  HLTTrigger(" << itrig << "): " << trigName << " = " << accept << std::endl; 
      }

      if((trigName == "HLT1Electron") && (accept == true)){ 
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "hltL1IsoSingleElectronTrackIsolFilter";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLT1ElectronRelaxed") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "hltL1NonIsoSingleElectronTrackIsolFilter";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLT2Electron") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "hltL1IsoDoubleElectronTrackIsolFilter";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLT2ElectronRelaxed") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "hltL1NonIsoDoubleElectronTrackIsolFilter";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLT1MuonIso") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "SingleMuIsoL3IsoFiltered";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLT1MuonNonIso") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "SingleMuNoIsoL3PreFiltered";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "CandHLT2MuonIso") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "DiMuonIsoL3IsoFiltered";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLT2MuonNonIso") && (accept == true)){
	EvtView.set().setUserRecord<bool>(trigName, accept);
	hltobject = "DiMuonNoIsoL3PreFiltered";
	saveHLTobjects(iEvent, EvtView, hltobject);
      }
      if((trigName == "HLTNMuonNonIso") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLT1EMHighEt") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLT1EMVeryHighEt") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXElectronMuon") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXElectronMuonRelaxed") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXElectron1Jet") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXElectron2Jet") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXElectron3Jet") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXElectron4Jet") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      if((trigName == "HLTXMuonJets") && (accept == true)) EvtView.set().setUserRecord<bool>(trigName, accept);
      
    }
  }

  //L1 trigger bits
  // Check if the double tau trigger fired.
  edm::Handle<l1extra::L1ParticleMapCollection> l1mapColl ;
  iEvent.getByLabel( "l1extraParticleMap", l1mapColl ) ;

  for(l1extra::L1ParticleMapCollection::const_iterator  l1bit = l1mapColl->begin(); l1bit != l1mapColl->end(); ++l1bit ) {
    //const L1ParticleMap& doubleTauMap = ( *mapColl )[ L1ParticleMap::L1_SingleMu7 ] ;
    //bool singleTauFired = doubleTauMap.triggerDecision() ;
    const string& trigName = l1bit->triggerName();
    bool accept = l1bit->triggerDecision();

    if (fDebug > 1) std::cout << "%L1Info --  L1Trigger: " << trigName << " = " << accept << std::endl;

    if((trigName == "L1_SingleMu7") && (accept == true)) EvtView.set().setUserRecord<bool>("L1_SingleMu7", accept);
    if((trigName == "L1_DoubleMu3") && (accept == true)) EvtView.set().setUserRecord<bool>("L1_DoubleMu3", accept);
    if((trigName == "L1_SingleIsoEG12") && (accept == true)) EvtView.set().setUserRecord<bool>("L1_SingleIsoEG12", accept);
    if((trigName == "L1_SingleEG15") && (accept == true)) EvtView.set().setUserRecord<bool>("L1_SingleEG15", accept);
    if((trigName == "L1_DoubleIsoEG8") && (accept == true)) EvtView.set().setUserRecord<bool>("L1_DoubleIsoEG8", accept);
    if((trigName == "L1_DoubleEG10") && (accept == true)) EvtView.set().setUserRecord<bool>("L1_DoubleEG10", accept);
  }

}

// ------------ Function to get the HLT Objects through the HLTFilterObjectWithRefs and save as pxl-partciles ------------------
  
void ePaxAnalyzer::saveHLTobjects(const edm::Event& iEvent, pxl::EventViewRef EvtView, string& hltobject){
  edm::Handle<HLTFilterObjectWithRefs> theHltFilterObjectHandle;
  // if( iEvent.getByLabel(hltobject, theHltFilterObjectHandle) ){ //bool works only >17x, hav to use ugly try...catch for the moment
  try {
    iEvent.getByLabel(hltobject, theHltFilterObjectHandle);
    for (unsigned int i=0; i<theHltFilterObjectHandle->size(); i++) {
      edm::RefToBase<reco::Candidate> candref = theHltFilterObjectHandle->getParticleRef(i);
      //save HLT-object info as pxl particle, use label string as name for particle
      pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
      part.set().setName(hltobject);
      part.set().setCharge(candref->charge());
      part.set().vector(pxl::set).setPx(candref->px());
      part.set().vector(pxl::set).setPy(candref->py());
      part.set().vector(pxl::set).setPz(candref->pz());
      part.set().vector(pxl::set).setE(candref->energy());
      /*cout <<"HLT Collection with label "<< hltobject << endl;
      cout << "theHltFilterObjectVector[" << i << "].pt()  = " << candref->pt()  << endl;
      cout << "theHltFilterObjectVector[" << i << "].eta() = " << candref->eta() << endl;
      cout << "theHltFilterObjectVector[" << i << "].phi() = " << candref->phi() << endl;*/
    }
    } catch (...) {
      //cout <<"No HLT Collection with label "<< hltobject << endl;
    }
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
	 //save info on quality of track-fit
         part.set().setUserRecord<double>("NormChi2", muon->combinedMuon()->normalizedChi2());
         part.set().setUserRecord<int>("ValidHits", muon->combinedMuon()->numberOfValidHits());
         part.set().setUserRecord<int>("LostHits", muon->combinedMuon()->numberOfLostHits()); 
	 //error info also used in muon-Met corrections, thus store variable to save info for later re-corrections
	 part.set().setUserRecord<double>("dPtRelTrack", muon->combinedMuon()->error(0)/(muon->combinedMuon()->qoverp()) );
	 part.set().setUserRecord<double>("dPtRelTrack_off", muon->combinedMuon()->ptError()/muon->combinedMuon()->pt() );
	 part.set().setUserRecord<double>("NormChi2_SA", muon->standAloneMuon()->normalizedChi2());
         part.set().setUserRecord<int>("ValidHits_SA", muon->standAloneMuon()->numberOfValidHits());
         part.set().setUserRecord<int>("LostHits_SA", muon->standAloneMuon()->numberOfLostHits()); 
	 part.set().setUserRecord<double>("NormChi2_TM", muon->track()->normalizedChi2());
         part.set().setUserRecord<int>("ValidHits_TM", muon->track()->numberOfValidHits());
         part.set().setUserRecord<int>("LostHits_TM", muon->track()->numberOfLostHits());
	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part.set().setUserRecord<double>("DSZ", muon->combinedMuon()->dsz()); 
	 //part.set().setUserRecord<double>("DZ", muon->combinedMuon()->dz()); //not needed since vz()==dz()
       
	 // TEMPORARY: calculate isolation ourselves
	 double CaloIso = IsoCalSum(iEvent, 0., muon->track()->eta(), muon->track()->phi(), 0.3, 1.5);
	 double TrkIso = IsoTrkSum(iEvent,  muon->track()->pt(), muon->track()->eta(), muon->track()->phi(), 0.3, 1.5);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);
 
	 //save offical isolation information
	 const MuonIsolation& muonIsoR03 = muon->getIsolationR03();
	 part.set().setUserRecord<double>("IsoR03SumPt", muonIsoR03.sumPt);
	 part.set().setUserRecord<double>("IsoR03EmEt", muonIsoR03.emEt);
	 part.set().setUserRecord<double>("IsoR03HadEt", muonIsoR03.hadEt);
	 part.set().setUserRecord<double>("IsoR03HoEt", muonIsoR03.hoEt);
	 part.set().setUserRecord<int>("IsoR03NTracks", muonIsoR03.nTracks);
	 part.set().setUserRecord<int>("IsoR03NJets", muonIsoR03.nJets);
 
	 const MuonIsolation& muonIsoR05 = muon->getIsolationR05();
	 part.set().setUserRecord<double>("IsoR05SumPt", muonIsoR05.sumPt);
	 part.set().setUserRecord<double>("IsoR05EmEt", muonIsoR05.emEt);
	 part.set().setUserRecord<double>("IsoR05HadEt", muonIsoR05.hadEt);
	 part.set().setUserRecord<double>("IsoR05HoEt", muonIsoR05.hoEt);
	 part.set().setUserRecord<int>("IsoR05NTracks", muonIsoR05.nTracks);
	 part.set().setUserRecord<int>("IsoR05NJets", muonIsoR05.nJets);

	 //save some stuff related to Muon-ID (Calo-info etc.)
	 part.set().setUserRecord<double>("CaloCompatibility", muon->getCaloCompatibility());
	 part.set().setUserRecord<int>("NumberOfChambers", muon->numberOfChambers());
	 part.set().setUserRecord<int>("NumberOfMatches", muon->numberOfMatches());
	 part.set().setUserRecord<double>("MuonDepositEM", muon->getCalEnergy().em);
	 part.set().setUserRecord<double>("MuonDepositHCAL", muon->getCalEnergy().had);
	
         numMuonRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumMuon", numMuonRec);
  
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 

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
   int numPixelEleAll = 0; //need counter for EleID 
   reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;

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

	 //save sz-distance to (0,0,0) for checking compability with primary vertex (need dsz or dz???)
	 part.set().setUserRecord<double>("DSZ", ele->gsfTrack()->dsz()); 
	 //part.set().setUserRecord<double>("DZ", ele->gsfTrack()->dz()); //not needed since vz()==dz()
         
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

	 // Read in electron ID association map and save ID decision
	 edm::Handle<reco::ElectronIDAssociationCollection> electronIDAssocHandle;
	 iEvent.getByLabel(fElectronIDAssocProducer, electronIDAssocHandle);
	 //need Ref of eletron for association-map, use counter of ALL electrons for that
	 edm::Ref<reco::PixelMatchGsfElectronCollection> electronRef(pixelelectrons, numPixelEleAll);
	 // Find entry in electron ID map corresponding electron
	 reco::ElectronIDAssociationCollection::const_iterator electronIDAssocItr;
	 electronIDAssocItr = electronIDAssocHandle->find(electronRef);
	 const reco::ElectronIDRef& eleid = electronIDAssocItr->val;
	 bool cutBasedID = eleid->cutBasedDecision();
	 part.set().setUserRecord<bool>("CutBasedID", cutBasedID);

	 // TEMPORARY: calculate isolation ourselves
	 double CaloPt = ( ele->superCluster()->rawEnergy() + ele->superCluster()->rawEnergy()*ele->hadronicOverEm() ) / cosh(ele->trackMomentumAtCalo().eta());
	 double CaloIso = IsoCalSum(iEvent, CaloPt, ele->trackMomentumAtCalo().eta(), ele->trackMomentumAtCalo().phi(), 0.3, 1.5);
	 double TrkIso = IsoTrkSum(iEvent, ele->gsfTrack()->pt(), ele->gsfTrack()->eta(), ele->gsfTrack()->phi(), 0.3, 1.5);
	 //double TrkIso = IsoTrkSum(iEvent, ele->trackMomentumAtVtx().Rho(), ele->trackMomentumAtVtx().Eta(), ele->trackMomentumAtVtx().Phi(), 0.2, 0.1);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);

	 //save official isolation information

	 // Get the association vector for HCAL isolation
	 edm::Handle<reco::CandViewDoubleAssociations> hcalIsolationHandle;
	 iEvent.getByLabel(fElectronHcalIsolationProducer, hcalIsolationHandle);

	 //direct access (by object). Since we do not really use candidates have to manually convert a reco-Ref into CandidateBaseRef
	 edm::Ref<reco::PixelMatchGsfElectronCollection> electronIsoRef(pixelelectrons, numPixelEleAll);
	 reco::CandidateBaseRef candRef(electronIsoRef);
	 double isoVal = (*hcalIsolationHandle)[ candRef ];
	 part.set().setUserRecord<double>("HCALIso", isoVal);

	 // Get the association vector for ECAL isolation
	 edm::Handle<reco::CandViewDoubleAssociations> ecalIsolationHandle;
	 iEvent.getByLabel(fElectronEcalIsolationProducer, ecalIsolationHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 isoVal = (*ecalIsolationHandle)[ candRef ];
	 part.set().setUserRecord<double>("ECALIso", isoVal);

	 // Get the association vector for track isolation
	 edm::Handle<reco::PMGsfElectronIsoCollection> trackIsolationHandle;
	 iEvent.getByLabel(fElectronTrackIsolationProducer, trackIsolationHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 isoVal = (*trackIsolationHandle)[ electronIsoRef ];
	 part.set().setUserRecord<double>("TrackIso", isoVal);

	 // Get the association vector for track number
	 edm::Handle<reco::PMGsfElectronIsoNumCollection> trackNumHandle;
	 iEvent.getByLabel(fElectronTrackNumProducer, trackNumHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 int numVal = (*trackNumHandle)[ electronIsoRef ];
	 part.set().setUserRecord<int>("TrackNum", numVal);

	 	 
         numPixelEleRec++;
      }
      numPixelEleAll++;
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
   edm::Handle<reco::CaloJetCollection> L2L3JESic5jets;

   iEvent.getByLabel(fKtJetRecoLabel, Ktjets);
   iEvent.getByLabel(fItCone5JetRecoLabel, ItCone5jets);
   iEvent.getByLabel(fL2L3JESic5JetRecoLabel, L2L3JESic5jets);
   
   int numKtJetRec = 0;
   int numItCone5JetRec = 0;
   int numL2L3JESic5JetRec = 0;

   //get primary vertex (hopefully correct one) for physics eta
   double VertexZ = 0.;
   if (EvtView().findUserRecord<int>("NumVertices") > 0) {
      pxl::Objects::TypeIterator<pxl::Vertex> iter(EvtView().getObjects()); 
      pxl::VertexWkPtr vtx = iter.object();
      if(vtx.valid()) VertexZ = vtx.get().vector().getZ();
   } 

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
	 part.set().setUserRecord<int>("N90", jet->n90());
	 part.set().setUserRecord<int>("N60", jet->n60());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
	 part.set().setUserRecord<double>("TowersArea", jet->towersArea());
	 part.set().setUserRecord<double>("PhysicsEta", jet->physicsEtaDetailed(VertexZ));
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
	 part.set().setUserRecord<int>("N90", jet->n90());
	 part.set().setUserRecord<int>("N60", jet->n60());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
	 part.set().setUserRecord<double>("TowersArea", jet->towersArea());
	 part.set().setUserRecord<double>("PhysicsEta", jet->physicsEtaDetailed(VertexZ));
         numItCone5JetRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumItCone5Jet", numItCone5JetRec);

   //same as above but this time corrected for JES using L2/L3 sheme
   for(reco::CaloJetCollection::const_iterator jet = L2L3JESic5jets->begin();
           jet != L2L3JESic5jets->end(); ++jet ) {
      if (Jet_cuts(jet)) {
         pxl::ParticleRef part = EvtView.set().create<pxl::Particle>();
         part.set().setName("L2L3JESic5Jet");
         part.set().vector(pxl::set).setPx(jet->px());
         part.set().vector(pxl::set).setPy(jet->py());
         part.set().vector(pxl::set).setPz(jet->pz());
         part.set().vector(pxl::set).setE(jet->energy());
	 part.set().setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
	 part.set().setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
	 part.set().setUserRecord<int>("N90", jet->n90());
	 part.set().setUserRecord<int>("N60", jet->n60());
	 part.set().setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
	 part.set().setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
	 part.set().setUserRecord<double>("TowersArea", jet->towersArea());
	 part.set().setUserRecord<double>("PhysicsEta", jet->physicsEtaDetailed(VertexZ));
         numL2L3JESic5JetRec++;
      }
   }
   EvtView.set().setUserRecord<int>("NumL2L3JESic5Jet", numL2L3JESic5JetRec);

   if (fDebug > 1) cout << "Found Rec Jets:  " << numKtJetRec << " Kt  " << numItCone5JetRec << " It5  " << endl;

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

   //****************************************************************************************************************************
   //*************** store corrected MET ****************************************************************************************
   //****************************************************************************************************************************
   edm::Handle<reco::CaloMETCollection> CaloMetCorr;
   iEvent.getByLabel(fMETCorrRecoLabel, CaloMetCorr);
   const CaloMETCollection *calometcorrcol = CaloMetCorr.product();
   const CaloMET calometcorr = calometcorrcol->front();  // MET exists only once!
 
   int numMETCorrRec = 0;
   pxl::ParticleRef partcorr = EvtView.set().create<pxl::Particle>();
   partcorr.set().setName("METCorr");
   partcorr.set().vector(pxl::set).setPx(calometcorr.px());
   partcorr.set().vector(pxl::set).setPy(calometcorr.py());
   partcorr.set().vector(pxl::set).setPz(0.);
   partcorr.set().vector(pxl::set).setMass(0.);
   partcorr.set().setUserRecord<double>("sumEt", calometcorr.sumEt());
   partcorr.set().setUserRecord<double>("mEtSig", calometcorr.mEtSig());
   partcorr.set().setUserRecord<double>("EmEt", calometcorr.emEtFraction());          //not muon corrected
   partcorr.set().setUserRecord<double>("HadEt", calometcorr.etFractionHadronic());   //not muon corrected
   partcorr.set().setUserRecord<double>("MaxEtEm", calometcorr.maxEtInEmTowers());   //not muon corrected
   partcorr.set().setUserRecord<double>("MaxEtHad", calometcorr.maxEtInHadTowers()); //not muon corrected

   //DO NOT CORRECT FOR MUONS MANUALLY, USE OFFICIAL MUONMET-CORRECTION!!!   

   //there is always MET in event, just decide if cuts passed (do this after muon corrections!)
   if (MET_cuts(partcorr)) { 
     numMETCorrRec++;
   }

   EvtView.set().setUserRecord<int>("NumMETCorr", numMETCorrRec);
   if (numMETCorrRec && fDebug > 1) cout << "Found RecMETCorr" << endl;
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
   /*edm::Handle<HBHERecHitCollection> hbhe;
   //HBHERecHitMetaCollection *HBHE_RecHits = 0;
   iEvent.getByLabel(fHBHELabel, fHBHEInstanceName, hbhe);  
   HBHERecHitMetaCollection HBHE_RecHits(*hbhe); */
   
   int numGammaRec = 0;
   int numGammaAll = 0; //need counter for EleID 
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

	 //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection) adn to reject converted photons
	 part.set().setUserRecord<double>("seedphi", SCRef->seed()->phi());
	 part.set().setUserRecord<double>("seedeta", SCRef->seed()->eta());
	 part.set().setUserRecord<int>("seedId", seedShapeRef->eMaxId().rawId());
         
	 // calculate HadoverEm - now its getting tough ...
	 //CaloConeSelector sel(0.1, theCaloGeom.product(), DetId::Hcal);
         //GlobalPoint pclu(SCRef->x(),SCRef->y(),SCRef->z());
	 double hcalEnergy = 0.;
         /*std::auto_ptr<CaloRecHitMetaCollectionV> chosen=sel.select(pclu,HBHE_RecHits);
         for (CaloRecHitMetaCollectionV::const_iterator i=chosen->begin(); i!=chosen->end(); i++) {
            //std::cout << HcalDetId(i->detid()) << " : " << (*i) << std::endl;
            hcalEnergy += i->energy();
	    }*/
	 edm::Handle<CaloTowerCollection> caloTowers;
	 iEvent.getByLabel("towerMaker",caloTowers);

	 for(CaloTowerCollection::const_iterator tower = caloTowers->begin(); tower!= caloTowers->end(); ++tower){
	   double dR = deltaR(tower->eta(), tower->phi(), photon->eta(), photon->phi());
	   if(dR < 0.1) hcalEnergy +=tower->hadEnergy(); //FIXME: hardcoded cut
	 }

         HoE = hcalEnergy/photon->energy();
         //cout << "H/E : " << HoE << endl;
         part.set().setUserRecord<float>("HoE", HoE);

	 //*************** store Gamma info corrected for primary vertex (this changes direction but leaves energy of SC unchanged *********
	 //get primary vertex (hopefully correct one) for physics eta
         math::XYZPoint vtx(0., 0., 0.);
         if (EvtView().findUserRecord<int>("NumVertices") > 0) {
	    pxl::Objects::TypeIterator<pxl::Vertex> iter(EvtView().getObjects()); 
	    pxl::VertexWkPtr pxlvtx = iter.object();
	    if(pxlvtx.valid()) vtx = math::XYZPoint(pxlvtx.get().vector().getX(), pxlvtx.get().vector().getY(), pxlvtx.get().vector().getZ());
         }
	 /////  Set event vertex
	 reco::Photon localPho(*photon);
	 localPho.setVertex(vtx);
	 part.set().setUserRecord<double>("PhysicsEta", localPho.p4().eta());
	 part.set().setUserRecord<double>("PhysicsPhi", localPho.p4().phi());

	 // TEMPORARY: calculate isolation ourselves
	 double CaloPt = ( hcalEnergy + photon->superCluster()->rawEnergy() ) / cosh(photon->eta());
	 double CaloIso = IsoCalSum(iEvent, CaloPt, photon->eta(), photon->phi(), 0.3, 1.5);
	 double TrkIso = IsoTrkSum(iEvent, 0., photon->eta(), photon->phi(), 0.3, 1.5);
	 part.set().setUserRecord<double>("CaloIso", CaloIso);
	 part.set().setUserRecord<double>("TrkIso", TrkIso);

   

	 //save official isolation information

	 // Get the association vector for HCAL isolation
	 edm::Handle<reco::CandViewDoubleAssociations> hcalIsolationHandle;
	 iEvent.getByLabel(fGammaHcalIsolationProducer, hcalIsolationHandle);

	 //direct access (by object). Since we do not really use candidates have to manually convert a reco-Ref into CandidateBaseRef
	 edm::Ref<reco::PhotonCollection> gammaIsoRef(Photons, numGammaAll);
	 reco::CandidateBaseRef candRef(gammaIsoRef);
	 double isoVal = (*hcalIsolationHandle)[ candRef ];
	 part.set().setUserRecord<double>("HCALIso", isoVal);

	 // Get the association vector for ECAL isolation
	 edm::Handle<reco::CandViewDoubleAssociations> ecalIsolationHandle;
	 iEvent.getByLabel(fGammaEcalIsolationProducer, ecalIsolationHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 isoVal = (*ecalIsolationHandle)[ candRef ];
	 part.set().setUserRecord<double>("ECALIso", isoVal);

	 // Get the association vector for track isolation
	 edm::Handle<reco::CandViewDoubleAssociations> trackIsolationHandle;
	 iEvent.getByLabel(fGammaTrackIsolationProducer, trackIsolationHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 isoVal = (*trackIsolationHandle)[ candRef ];
	 part.set().setUserRecord<double>("TrackIso", isoVal);

	 // Get the association vector for track number
	 edm::Handle<reco::CandViewIntAssociations> trackNumHandle;
	 iEvent.getByLabel(fGammaTrackNumProducer, trackNumHandle);

	 //direct access (by object). We have candidate already from HCAL-isolation
	 int numVal = (*trackNumHandle)[ candRef ];
	 part.set().setUserRecord<int>("TrackNum", numVal);


	 // store neural-network output for gamma/pi0 rejection. 
	 // WARNING: should only be used with isolated unconverted photons --> compare with SC from converted photons !!!
	 /*edm::Handle<reco::PhotonPi0DiscriminatorAssociationMap>  map;
	 iEvent.getByLabel("piZeroDiscriminators","PhotonPi0DiscriminatorAssociationMap",  map);
	 reco::PhotonPi0DiscriminatorAssociationMap::const_iterator mapIter;*/

	 // get the Converted Photon info
	 edm::Handle<ConvertedPhotonCollection> convertedPhotonHandle; 
	 iEvent.getByLabel("convertedPhotons", "", convertedPhotonHandle);
	 const reco::ConvertedPhotonCollection ConvPhotons = *(convertedPhotonHandle.product());
	 bool isPhotConv = false;
	 int Ntrk_conv = 0;
	 for( reco::ConvertedPhotonCollection::const_iterator iCPho = ConvPhotons.begin(); iCPho != ConvPhotons.end(); iCPho++) { 
	   SuperClusterRef it_superConv = (*iCPho).superCluster();// get the SC related to the  Converted Photon candidate
	   if(SCRef == it_superConv) { //check if photon candidate is identical to converted photon
	     isPhotConv = (*iCPho).isConverted(); 
	     Ntrk_conv = (*iCPho).tracks().size();
	     break;
	   }    
	 } // End of Photon Conversion Loop  
	 part.set().setUserRecord<bool>("IsConverted", isPhotConv);  
	 part.set().setUserRecord<int>("ConvertedNtrk", Ntrk_conv);  
	 if (fDebug > 1) {cout<<"isPhotConv: "<<isPhotConv<<endl;
	 cout<<"Ntrk_conv: "<<Ntrk_conv<<endl;}


	 //mapIter = map->find(edm::Ref<reco::PhotonCollection>(Photons, numGammaAll));
	 double nn = 1; //default: photon candidate is not a Pi0
	 //if(mapIter!=map->end()){//check if photon exists in map at all
	 // nn = mapIter->val; //discrimination variable has distribution peaking close to 0 for unconverted neutral pions and close to 1 for unconverted photons
	   if (fDebug > 1) cout<<"Pi0 dicriminant: "<<nn<<endl;
	   if(fabs(photon->eta()) <= 1.442) { //barrel with preshower-det
	     part.set().setUserRecord<double>("Pi0DisBarrel", nn);
	   } else if( (fabs(photon->eta()) >= 1.556 && fabs(photon->eta()) < 1.65) || fabs(photon->eta()) > 2.5) { //endcap without preshower-det  
	     part.set().setUserRecord<double>("Pi0DisEndcapNoPre", nn);
	   } else if(fabs(photon->eta()) >= 1.65 && fabs(photon->eta()) <= 2.5 ) { //endcap with preshower-det  
	     part.set().setUserRecord<double>("Pi0DisEndcap", nn);
	   } 
	   //}

	 	 
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
   if (muon->pt() < 15.)  return false;
   if (fabs(muon->eta()) > 3.) return false;
   return true;
}

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

// ------------ method to define GAMMA-cuts

bool ePaxAnalyzer::Gamma_cuts(reco::PhotonCollection::const_iterator photon) const {
   //
   if (photon->pt() < 15.) return false;
   if (fabs(photon->eta()) > 3.) return false;
   return true;
}

// ------------ method to define MC-MET-cuts

bool ePaxAnalyzer::MET_cuts(const pxl::ParticleRef met) const {
   // 
   if (met.get().vector().getPt() < 30.) return false;
   return true;
}


// TEMPORARY STUFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


//define this as a plug-in
DEFINE_FWK_MODULE(ePaxAnalyzer);
