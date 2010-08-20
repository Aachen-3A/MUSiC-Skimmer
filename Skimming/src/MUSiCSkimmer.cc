// -*- C++ -*-
//
// Package:    MUSiCSkimmer
// Class:      MUSiCSkimmer
// 
/**\class MUSiCSkimmer MUSiCSkimmer.cc PaxDemo/MUSiCSkimmer/src/MUSiCSkimmer.cc

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
#include "MUSiCProject/Skimming/interface/MUSiCSkimmer.h"

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
#include "FWCore/Utilities/interface/Exception.h"

// necessary objects:
#include "FWCore/Framework/interface/ESHandle.h"

// for electron shapes:
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

//HCAL Noise
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

//for GenParticles
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

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
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//test
#include "DataFormats/Common/interface/Ptr.h"

//ECAL
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// special stuff for sim truth of converted photons
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

//Jet Flavour
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

using namespace std;
using namespace edm;

//
// constructors and destructor
//
MUSiCSkimmer::MUSiCSkimmer(const edm::ParameterSet& iConfig) : fFileName(iConfig.getUntrackedParameter<string>("FileName")), fePaxFile(fFileName) {
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
   freducedBarrelRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
   freducedEndcapRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");

   //get the PSet that contains all jet PSets
   ParameterSet jets_pset = iConfig.getParameter< ParameterSet >( "jets" );
   //get the names of all sub-PSets
   vector< string > jet_names;
   jets_pset.getParameterSetNames( jet_names );
   //loop over the names of the jet PSets
   for( vector< string >::const_iterator jet_name = jet_names.begin(); jet_name != jet_names.end(); ++jet_name ){
      jet_def jet;
      jet.name = *jet_name;

      ParameterSet one_jet = jets_pset.getParameter< ParameterSet >( jet.name );
      jet.MCLabel   = one_jet.getParameter< InputTag >( "MCLabel" );
      jet.RecoLabel = one_jet.getParameter< InputTag >( "RecoLabel" );

      jet_infos.push_back( jet );
   }

   //JetIDs
   vector< string > id_names = jets_pset.getParameter< vector< string > >( "IDs" );
   unsigned int num_IDs = JetIDSelectionFunctor::N_QUALITY;
   if( num_IDs < id_names.size() ){
      cout << "Less JetIDs available than requested, using only available." << endl;
   } else if( num_IDs > id_names.size() ){
      cout << "More JetIDs available than requested, using only requested." << endl;
      num_IDs = id_names.size();
   }
   //ATTENTION: The following is REALLY ugly
   //Looping over enums is apparentlt not forseen in C++
   //Seems to be the only way to make the JetIDs configurable
   for( JetIDSelectionFunctor::Quality_t q = JetIDSelectionFunctor::Quality_t(0);
        q < JetIDSelectionFunctor::N_QUALITY;
        q = JetIDSelectionFunctor::Quality_t( q+1 ) ){
      pair< std::string, JetIDSelectionFunctor > ID( id_names[q], JetIDSelectionFunctor( JetIDSelectionFunctor::CRAFT08, q ) );
      jet_ids.push_back( ID );
   }
   

   // MET label
   fMETRecoLabel = iConfig.getUntrackedParameter<string>("METRecoLabel");

   //HCAL noise
   hcal_noise_label = iConfig.getParameter< InputTag >( "HCALNoiseInfo" );


   //get the PSet that contains all trigger PSets
   ParameterSet trigger_pset = iConfig.getParameter< ParameterSet >( "triggers" );

   vector< string > trigger_processes;
   trigger_pset.getParameterSetNames( trigger_processes );

   //loop over the names of the trigger PSets
   for( vector< string >::const_iterator trg_proc = trigger_processes.begin(); trg_proc != trigger_processes.end(); ++trg_proc ){
      trigger_group trigger;
      trigger.name = *trg_proc;
      
      ParameterSet one_trigger = trigger_pset.getParameter< ParameterSet >( trigger.name );
      trigger.process = one_trigger.getParameter< string >( "process" );

      trigger.L1_result = one_trigger.getParameter< InputTag >( "L1_result" );
      trigger.results = InputTag( one_trigger.getParameter< string >( "results" ), "", trigger.process );
      trigger.event   = InputTag( one_trigger.getParameter< string >( "event" ),   "", trigger.process );
      
      trigger.triggers_names = one_trigger.getParameter< vector< string > >("HLTriggers");

      triggers.push_back( trigger );
   }

   fStoreL3Objects = trigger_pset.getUntrackedParameter<bool>("StoreL3Objects");



   //cuts
   ParameterSet cut_pset = iConfig.getParameter< ParameterSet >( "cuts" );
   min_muon_pt = cut_pset.getParameter< double >( "min_muon_pt" );
   min_ele_pt = cut_pset.getParameter< double >( "min_ele_pt" );
   min_gamma_pt = cut_pset.getParameter< double >( "min_gamma_pt" );
   min_jet_pt = cut_pset.getParameter< double >( "min_jet_pt" );
   min_met = cut_pset.getParameter< double >( "min_met" );
   max_eta = cut_pset.getParameter< double >( "max_eta" );
   min_rechit_energy = cut_pset.getParameter< double >( "min_rechit_energy" );
   min_rechit_swiss_cross = cut_pset.getParameter< double >( "min_rechit_swiss_cross" );
   min_rechit_R19 = cut_pset.getParameter< double >( "min_rechit_R19" );
   vertex_minNDOF = cut_pset.getParameter< double >( "vertex_minNDOF" );
   vertex_maxZ = cut_pset.getParameter< double >( "vertex_maxZ" );
   vertex_maxR = cut_pset.getParameter< double >( "vertex_maxR" );
   PV_minNDOF = cut_pset.getParameter< double >( "PV_minNDOF" );
   PV_maxZ = cut_pset.getParameter< double >( "PV_maxZ" );
   PV_maxR = cut_pset.getParameter< double >( "PV_maxR" );



   Matcher = new ParticleMatcher();
   fNumEvt=0;
}

// ------------ MIS Destructor  ------------

MUSiCSkimmer::~MUSiCSkimmer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   delete Matcher;
}

// ------------ method called to for each event  ------------

void MUSiCSkimmer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   // set event counter   
   fNumEvt++; 
   // Owner of all Pxl Objects 
   pxl::Event event;

   // event-specific data
   bool IsMC =  !iEvent.isRealData();
   event.setUserRecord<bool>("MC", IsMC);  //distinguish between MC and data
   event.setUserRecord< unsigned int >( "Run", iEvent.run() );
   event.setUserRecord< unsigned int >( "LumiSection", iEvent.luminosityBlock());
   event.setUserRecord< unsigned int >( "EventNum", iEvent.id().event() );
   event.setUserRecord< unsigned int >( "BX", iEvent.bunchCrossing() );
   event.setUserRecord< unsigned int >( "Orbit", iEvent.orbitNumber() );

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
   std::map<const Candidate*, pxl::Particle*> genmap; 
   std::map<const Candidate*, pxl::Particle*> genjetmap;

   //set process name
   GenEvtView->setUserRecord<std::string>("Process", fProcess);
   RecEvtView->setUserRecord<std::string>("Process", fProcess);

   //create object for EcalClusterLazyTools
   EcalClusterLazyTools lazyTools( iEvent, iSetup, freducedBarrelRecHitCollection, freducedEndcapRecHitCollection);
   
   // Generator stuff
   if (IsMC) {
      analyzeGenInfo(iEvent, GenEvtView, genmap);
      analyzeGenRelatedInfo(iEvent, GenEvtView);  // PDFInfo, Process ID, scale, pthat
      for( vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         analyzeGenJets(iEvent, GenEvtView, genjetmap, *jet_info);
      }
      analyzeGenMET(iEvent, GenEvtView);
      if (fUseSIM){
         analyzeSIM(iEvent, GenEvtView);
      }
   }
   // store Rec Objects only if requested
   if (!fGenOnly) {
      //Trigger bits
      for( vector< trigger_group >::iterator trg = triggers.begin(); trg != triggers.end(); ++trg ){
         analyzeTrigger( iEvent, iSetup, RecEvtView, *trg );
      }
      // Reconstructed stuff
      analyzeRecVertices(iEvent, RecEvtView);
      analyzeRecMuons(iEvent, RecEvtView, IsMC, genmap);
      analyzeRecElectrons(iEvent, RecEvtView, IsMC, lazyTools, genmap);
      for( vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         analyzeRecJets( iEvent, RecEvtView, IsMC, genjetmap, *jet_info );
      }
      analyzeRecMET(iEvent, RecEvtView);
      analyzeRecGammas(iEvent, RecEvtView, IsMC, lazyTools, genmap);
      analyzeECALRecHits( iEvent, iSetup, RecEvtView );
   }

   if (IsMC){
      const string met_name = "MET";
      Matcher->matchObjects(GenEvtView, RecEvtView, jet_infos, met_name);
   }

   if (fDebug > 0) {  
      cout << "UserRecord  " <<  "   e   " << "  mu   " << " Gamma ";
      for( std::vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         cout << jet_info->name << " ";
      }
      cout << "  MET  " << endl;
      cout << "Found (MC): " << setw(4) << GenEvtView->findUserRecord<int>("NumEle") 
           << setw(7) << GenEvtView->findUserRecord<int>("NumMuon")
           << setw(7) << GenEvtView->findUserRecord<int>("NumGamma");
      for( std::vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         cout << setw(4) << GenEvtView->findUserRecord<int>( "Num"+jet_info->name ) << " ";
      } 
      cout << setw(7) << GenEvtView->findUserRecord<int>("NumMET") << endl;
      cout << "     (Rec): " << setw(4) << RecEvtView->findUserRecord<int>("NumEle")    
           << setw(7) << RecEvtView->findUserRecord<int>("NumMuon") 
           << setw(7) << RecEvtView->findUserRecord<int>("NumGamma");
      for( std::vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         cout << setw(4) << RecEvtView->findUserRecord<int>( "Num"+jet_info->name ) << " ";
      } 
      cout << setw(7) << RecEvtView->findUserRecord<int>("NumMET") << endl;
   }   
   fePaxFile.writeEvent(&event);
}




void MUSiCSkimmer::beginRun( const edm::Run &iRun, const edm::EventSetup &iSetup ) {
   //we are in a new run, so we might have a new trigger config
   //hence we're going to check those guys
   for( vector< trigger_group >::iterator trg = triggers.begin(); trg != triggers.end(); ++trg ){
      trigger_group &trigger = *trg;

      //read the new trigger config, test for error and whether something has changed
      bool changed;
      if( ! trigger.config.init( iRun, iSetup, trigger.process, changed ) ){
         cout << "TRIGGER ERROR: Initialization of trigger config failed." << endl;
         throw "TRIGGER ERROR: Initialization of trigger config failed.";
      }

      //the trigger config has actually changed, so read in the new one
      if( changed ){
         cout << "TRIGGER INFO: HLT table changed, building new trigger map." << endl;
         //reset the map
         trigger.trigger_infos.clear();

         for( vector<string>::const_iterator trg_name = trigger.triggers_names.begin(); trg_name != trigger.triggers_names.end(); ++trg_name ){
            //get the number of the trigger path
            unsigned int index = trigger.config.triggerIndex( *trg_name );

            //check if that's a valid number
            if( index < trigger.config.size() ){
               //it is, so store the name and the number
               trigger_def trg;
               trg.name = *trg_name;
               trg.ID = index;
               trg.active = true;
               trigger.trigger_infos.push_back( trg );
            } else {
               //the number is invalid, the trigger path is not in the config
               cout << "TRIGGER WARNING: In run " << iRun.run() << " trigger " << *trg_name << " not found in HLT config, not added to trigger map (so not used)." << endl;
            }
         }
      }
   }
}





// ------------ reading Generator related Stuff ------------

void MUSiCSkimmer::analyzeGenRelatedInfo(const edm::Event& iEvent, pxl::EventView* EvtView) {
   // this works at least for RECO. Need to check if this works on AOD or PAT-Ntuplee 
  
   edm::Handle< GenEventInfoProduct > genEvtInfo;
   iEvent.getByLabel( "generator", genEvtInfo );

   //if the sample is binned, there should be a binning value. so save it, otherwise just save a 0
   EvtView->setUserRecord< double >( "binScale", genEvtInfo->hasBinningValues() ? genEvtInfo->binningValues()[0] : 0 );

   EvtView->setUserRecord< double >( "Weight", genEvtInfo->weight() );

   unsigned int ID = genEvtInfo->signalProcessID();
   EvtView->setUserRecord< unsigned int >( "processID", ID );

   //don't save PDF infos for processes without partons
   if( genEvtInfo->hasPDF() && !( 91 <= ID && ID <= 95 ) ){
      const gen::PdfInfo *pdf = genEvtInfo->pdf();
      fpdf_vec.push_back( *pdf );

      int id1 = pdf->id.first;
      int id2 = pdf->id.second;

      //reset the code for a gluon, at least SHERPA got a problem there
      if( abs( id1 ) == 9 || abs( id1 ) == 21 ) {
         id1 = 0;
         fpdf_vec.back().id.first = 0;
      }
      if( abs( id2 ) == 9 || abs( id2 ) == 21 ) {
         id2 = 0;
         fpdf_vec.back().id.second = 0;
      }

      EvtView->setUserRecord<float>("x1", pdf->x.first);
      EvtView->setUserRecord<float>("x2", pdf->x.second);
      EvtView->setUserRecord<float>("Q", pdf->scalePDF);
      EvtView->setUserRecord<int>("f1", id1);
      EvtView->setUserRecord<int>("f2", id2);
      EvtView->setUserRecord<float>("pdf1", pdf->xPDF.first);
      EvtView->setUserRecord<float>("pdf2", pdf->xPDF.second);


      if( abs( id1 ) > 6 || abs( id2 ) > 6 ){
         throw cms::Exception( "PDF error" ) << "PDF information corrupted in a non-diffractive event." << endl
                                             << "Process ID " << genEvtInfo->signalProcessID() << " is not in list of diffractive processes (91 <= ID <= 95)." << endl
                                             << "Scale: " << pdf->scalePDF << endl
                                             << "x1 = " << pdf->x.first << "; x2 = " << pdf->x.second << endl
                                             << "ID 1: " << id1 << endl
                                             << "ID 2: " << id2 << endl;
      }
   } else {
      gen::PdfInfo pdf;
      pdf.scalePDF = 0;

      fpdf_vec.push_back( pdf );
   }

   

   if (fDebug > 0) {
      cout << "Event Scale (i.e. pthat): " << (genEvtInfo->hasBinningValues() ? genEvtInfo->binningValues()[0] : 0) << ", EventWeight: " << genEvtInfo->weight() << endl;
      if( genEvtInfo->hasPDF() ){
         cout << "PDFInfo: " << endl 
              << "========" << endl;
         cout << "Momentum of first incoming parton: (id/flavour = "  << genEvtInfo->pdf()->id.first  << ") "
              << genEvtInfo->pdf()->x.first  << endl
              << "Momentum of second incoming parton: (id/flavour = " << genEvtInfo->pdf()->id.second << ") "
              << genEvtInfo->pdf()->x.second << endl
              << "Scale = " << genEvtInfo->pdf()->scalePDF << endl;
      } else {
         cout << "No PDFInfo in this event." << endl;
      }
   }
}

// ------------ reading the Generator Stuff ------------

void MUSiCSkimmer::analyzeGenInfo(const edm::Event& iEvent, pxl::EventView* EvtView, std::map<const Candidate*, pxl::Particle*>& genmap ) {
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

void MUSiCSkimmer::analyzeGenJets( const edm::Event &iEvent, pxl::EventView *EvtView, std::map< const Candidate*, pxl::Particle* > &genjetmap, const jet_def &jet_info ) {
   //Get the GenJet collections
   edm::Handle<reco::GenJetCollection> GenJets;
   iEvent.getByLabel( jet_info.MCLabel, GenJets );

   //get the flavours
   edm::Handle< reco::JetFlavourMatchingCollection > algoFlavour, physicsFlavour;
   iEvent.getByLabel( jet_info.name+"GenJetFlavourAlgo", algoFlavour );
   iEvent.getByLabel( jet_info.name+"GenJetFlavourPhysics", physicsFlavour );


   // counter
   size_t jet_index = 0;
   int numJetMC = 0;
   double constit_pT = 5.; //here we have a hardcoded cut, but do we really need cfg-parameter for this?...
   
   //Loop over GenJets
   for( reco::GenJetCollection::const_iterator genJet = GenJets->begin(); genJet != GenJets->end(); ++genJet, jet_index++ ){
      if( JetMC_cuts( genJet ) ){
         //get the reference
         RefToBase< reco::Jet > jetRef( RefToBaseProd< reco::Jet >( GenJets ), jet_index );

         pxl::Particle *part = EvtView->create< pxl::Particle >();
         
         //cast iterator into GenParticleCandidate
         const GenParticle *p = dynamic_cast< const GenParticle* >( &(*genJet) );
         genjetmap[p] = part;
         part->setName( jet_info.name );
         part->setP4(genJet->px(), genJet->py(), genJet->pz(), genJet->energy());
         //fill additional jet-related infos
         part->setUserRecord<double>("EmE", genJet->emEnergy());
         part->setUserRecord<double>("HadE", genJet->hadEnergy());
         part->setUserRecord<double>("InvE", genJet->invisibleEnergy());
         part->setUserRecord<double>("AuxE", genJet->auxiliaryEnergy());
         numJetMC++;
         
         //save number of GenJet-constituents fulfilling some cuts
         int numGenJetConstit_withcuts = 0;
         const vector< const GenParticle* > &genJetConstit = genJet->getGenConstituents();
         for( std::vector< const GenParticle* >::const_iterator constit = genJetConstit.begin(); constit != genJetConstit.end(); ++constit ) {
            //raise counter if cut passed
            if( (*constit)->pt() > constit_pT ) numGenJetConstit_withcuts++; 
         }
         part->setUserRecord< int >( "GenJetConstit", numGenJetConstit_withcuts );

         part->setUserRecord< int >( "algoFlavour",    (*algoFlavour)   [ jetRef ].getFlavour() );
         part->setUserRecord< int >( "physicsFlavour", (*physicsFlavour)[ jetRef ].getFlavour() );
      }
   }
   EvtView->setUserRecord< int >( "Num"+jet_info.name, numJetMC );
   if( fDebug > 1 ) cout << "Found MC Jets:  "  << numJetMC << " of Type " << jet_info.name << endl;
}

// ------------ reading the Generator MET ------------

void MUSiCSkimmer::analyzeGenMET(const edm::Event& iEvent, pxl::EventView* EvtView) {
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
      pxl::ParticleFilter::apply( EvtView->getObjectOwner(), GenMuons, pxl::ParticlePtEtaNameCriterion("Muon") );
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
void MUSiCSkimmer::analyzeSIM(const edm::Event& iEvent, pxl::EventView* EvtView) {
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

void MUSiCSkimmer::analyzeRecMET(const edm::Event& iEvent, pxl::EventView* EvtView) {
   //save HCAL noise infos
   edm::Handle< HcalNoiseSummary > hcal_noise;
   iEvent.getByLabel( hcal_noise_label, hcal_noise );
   EvtView->setUserRecord< bool >( "HCALNoiseLoose", hcal_noise->passLooseNoiseFilter() );
   EvtView->setUserRecord< bool >( "HCALNoiseTight", hcal_noise->passTightNoiseFilter() );
   EvtView->setUserRecord< bool >( "HCALNoiseHighLevel", hcal_noise->passHighLevelNoiseFilter() );

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

void MUSiCSkimmer::analyzeTrigger( const edm::Event &iEvent,
                                   const edm::EventSetup &iSetup,
                                   pxl::EventView* EvtView,
                                   trigger_group &trigger
                                   ){
   edm::Handle<edm::TriggerResults>   triggerResultsHandle;
   iEvent.getByLabel( trigger.results, triggerResultsHandle );
   edm::Handle<trigger::TriggerEvent> triggerEventHandle;
   iEvent.getByLabel( trigger.event, triggerEventHandle );

   //loop over selected trigger names
   for( vector< trigger_def >::iterator trig = trigger.trigger_infos.begin(); trig != trigger.trigger_infos.end(); ++trig ){
      //skip this trigger if it's not active, e.g. because it's prescaled
      if( !trig->active ) continue;

      //get trigger path status
      bool wasrun = triggerResultsHandle->wasrun( trig->ID );
      bool error = triggerResultsHandle->error( trig->ID );
      unsigned int prescale = 0;

      //check that the trigger was run and not in error
      if( wasrun && !error ){
         //get the current prescale value
         prescale = trigger.config.prescaleValue( iEvent, iSetup, trig->name );
         
         //we can only use unprescaled triggers
         if( prescale == 1 ) {
            //unprescaled, so store it
            EvtView->setUserRecord< bool >( trigger.name+"_"+trig->name, triggerResultsHandle->accept( trig->ID ) );

            //debug output
            if( fDebug > 0 && triggerResultsHandle->accept( trig->ID ) )
               cout << endl << "Trigger: " << trig->name << " in menu " << trigger.process << " fired" << endl;
         } else {
            //prescaled!
            //switch it off
            trig->active = false;
            cout << "TRIGGER WARNING: Prescaled " << trig->name << " in menu " << trigger.process << endl;
            cout << "TRIGGER WARNING: Run " << iEvent.run() << " - LS " << iEvent.luminosityBlock() << " - Event " << iEvent.id().event() << endl;
            cout << "TRIGGER WARNING: " << trig->name << " switched to inactive." << endl;
         }
      } else {
         //either error or was not run
         if( !wasrun ) cout << "TRIGGER WARNING: Trigger: " << trig->name << " in menu " << trigger.process << " was not executed!" << endl;
         if( error ) cout << "TRIGGER WARNING: An error occured during execution of Trigger: " << trig->name << " in menu " << trigger.process << endl;
         cout << "TRIGGER WARNING: Run " << iEvent.run() << " - LS " << iEvent.luminosityBlock() << " - Event " << iEvent.id().event() << endl;
      }

      //begin cout of saved information for debugging
      if( fDebug > 1 ){
         cout << "triggerName: " << trig->name << "  triggerIndex: " << trig->ID << endl;
         cout << " Trigger path status:"
              << " WasRun=" << wasrun
              << " Accept=" << triggerResultsHandle->accept( trig->ID )
              << " Error=" << error
              << " Prescale=" << prescale << endl;
      }
      if( fStoreL3Objects ){
         const vector<string> &moduleLabels( trigger.config.moduleLabels( trig->ID ) );
         const unsigned int moduleIndex( triggerResultsHandle->index( trig->ID) );

         // Results from TriggerEvent product - Attention: must look only for
         // modules actually run in this path for this event!
         // Analyze and store the objects which have fired the trigger
         for( unsigned int j = 0; j <= moduleIndex; ++j ){
            const string &moduleLabel = moduleLabels[j];
            // check whether the module is packed up in TriggerEvent product
            const unsigned int filterIndex = triggerEventHandle->filterIndex( InputTag( moduleLabel, "", trigger.process ) );
            if( filterIndex < triggerEventHandle->sizeFilters() ){
               const trigger::Vids &VIDS( triggerEventHandle->filterIds(  filterIndex ) );
               const trigger::Keys &KEYS( triggerEventHandle->filterKeys( filterIndex ) );
               const size_t nI( VIDS.size() );
               const size_t nK( KEYS.size() );
               assert( nI==nK );
               size_t n( max( nI,nK ) );
               if( n > 5 ){
                  cout << "Storing only 5 L3 objects for label/type " << moduleLabel << "/" << trigger.config.moduleType( moduleLabel ) << endl;
                  n = 5;
               }
               const trigger::TriggerObjectCollection &TOC = triggerEventHandle->getObjects();
               for( size_t i = 0; i != n; ++i ){
                  const trigger::TriggerObject &TO = TOC[KEYS[i]];
                  pxl::Particle *part = EvtView->create< pxl::Particle >();
                  part->setName( moduleLabel );
                  part->setP4( TO.px(), TO.py(), TO.pz(), TO.energy() );
                  part->setUserRecord< double >( "ID", TO.id() );
                  if( fDebug > 0 )
                     cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
                          << " id: " << TO.id() << " pt: " << TO.pt() << " eta: "
                          << TO.eta() << " phi:" << TO.phi() << " m: " << TO.mass()
                          << endl;
               }
            }
         }
      } 
   }

   //get the L1 data
   edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
   iEvent.getByLabel( trigger.L1_result, gtReadoutRecord );
   //get the technical trigger word
   const TechnicalTriggerWord &tech_word = gtReadoutRecord->technicalTriggerWord();
   
   //store the important bits
   EvtView->setUserRecord< bool >( trigger.name+"_L1_0", tech_word[ 0 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_36", tech_word[ 36 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_37", tech_word[ 37 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_38", tech_word[ 38 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_39", tech_word[ 39 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_40", tech_word[ 40 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_41", tech_word[ 41 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_42", tech_word[ 42 ] );
   EvtView->setUserRecord< bool >( trigger.name+"_L1_43", tech_word[ 43 ] );
}

// ------------ reading Reconstructed Primary Vertices ------------

void MUSiCSkimmer::analyzeRecVertices(const edm::Event& iEvent, pxl::EventView* EvtView) {
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(fVertexRecoLabel, vertices);

   //get the beamspot
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByLabel("offlineBeamSpot",recoBeamSpotHandle);
   const reco::BeamSpot &beamspot = *recoBeamSpotHandle;

   //store the BeamSpot
   pxl::Vertex* bs = EvtView->create< pxl::Vertex >();
   bs->setName( "BeamSpot" );
   bs->setXYZ( beamspot.x0(), beamspot.y0(), beamspot.z0() );

   //get the PV
   const reco::Vertex &PV = *( vertices->begin() );

   //save the primary vertex postion for later use
   //use the BeamSpot in case the PV is shit
   if( PV_vertex_cuts( PV ) ) {
      the_vertex = PV.position();
   } else {
      the_vertex = beamspot.position();
   }
   
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
         vtx->setUserRecord<double>("chi2", vertex->chi2());
         vtx->setUserRecord<double>("ndof", vertex->ndof());
         // is valid?
         vtx->setUserRecord<bool>("IsValid", vertex->isValid());
         vtx->setUserRecord<bool>("IsFake", vertex->isFake());
         numVertices++;
      }
   }
   EvtView->setUserRecord<int>("NumVertices", numVertices); 
}

// ------------ reading Reconstructed Muons ------------

void MUSiCSkimmer::analyzeRecMuons(const edm::Event& iEvent, pxl::EventView* RecView, const bool& MC, std::map<const Candidate*, pxl::Particle*> & genmap) {
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
            std::map< const Candidate*, pxl::Particle* >::const_iterator it = genmap.find( muon->genLepton() );
            if( it != genmap.end() ){
               part->linkSoft( it->second, "pat-match" );
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
         part->setUserRecord<int>("TrackerVHits", trackerTrack->numberOfValidHits());
         part->setUserRecord<int>("TrackerLHits", trackerTrack->numberOfLostHits());
         //error info also used in muon-Met corrections, thus store variable to save info for later re-corrections
         part->setUserRecord<double>("dPtRelTrack", muontrack->error(0)/(muontrack->qoverp()));
         part->setUserRecord<double>("dPtRelTrack_off", muontrack->ptError()/muontrack->pt());
         // Save distance to the primary vertex in z and xy plane i.e. impact parameter
         //get primary vertex (hopefully correct one) for physics eta THIS NEEDS TO BE CHECKED !!
         // units given in cm!!! Use nominal beam spot of Summer/Fall08 if no vertex available
         math::XYZPoint vtx(0.0322, 0., 0.);
         if (RecView->findUserRecord<int>("NumVertices") > 0) {
            pxl::ObjectOwnerTypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
            vtx = math::XYZPoint((*vtx_iter)->getX(), (*vtx_iter)->getY(), (*vtx_iter)->getZ());
         } 
         part->setUserRecord<double>("Dsz", muontrack->dsz(vtx));
         part->setUserRecord<double>("Dxy", muontrack->dxy(vtx));
       
         //official CaloIso and TrkIso
         //Def:  aMuon.setCaloIso(aMuon.isolationR03().emEt + aMuon.isolationR03().hadEt + aMuon.isolationR03().hoEt);
         part->setUserRecord<double>("CaloIso", muon->caloIso());
         part->setUserRecord<double>("TrkIso", muon->trackIso());
         part->setUserRecord<double>("ECALIso", muon->ecalIso());
         part->setUserRecord<double>("HCALIso", muon->hcalIso());
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
         part->setUserRecord<bool>("isGood", muon::isGoodMuon( *muon, muon::GlobalMuonPromptTight ) );
         part->setUserRecord<bool>("lastStationTight", muon::isGoodMuon( *muon, muon::TMLastStationTight ) ); 
         part->setUserRecord<float>("SegComp", muon::segmentCompatibility( *muon ) );
         numMuonRec++;
      }
   }
   RecView->setUserRecord<int>("NumMuon", numMuonRec);
   if (fDebug > 1) cout << "Rec Muons: " << numMuonRec << endl; 
}

// ------------ reading Reconstructed Electrons ------------

void MUSiCSkimmer::analyzeRecElectrons(const edm::Event& iEvent, pxl::EventView* RecView, bool& MC, EcalClusterLazyTools& lazyTools, std::map<const Candidate*, pxl::Particle*> & genmap) {
   int numEleRec = 0;   
   int numEleAll = 0;   // for matching

   edm::Handle<std::vector<pat::Electron> > electronHandle;
   iEvent.getByLabel(fElectronRecoLabel, electronHandle);
   const std::vector<pat::Electron> &electrons = *electronHandle;

   edm::Handle< EcalRecHitCollection > barrelRecHits;
   iEvent.getByLabel( freducedBarrelRecHitCollection, barrelRecHits );

   for (std::vector<pat::Electron>::const_iterator ele = electrons.begin(); ele != electrons.end(); ++ele ) {
      if (Ele_cuts(ele)) {
         if (fDebug > 1) {
            cout << "Electron Energy scale corrected: " << ele->isEnergyScaleCorrected() 
                 << "  Momentum corrected: " << ele->isMomentumCorrected() << endl;
         }
         pxl::Particle* part = RecView->create<pxl::Particle>();
         part->setName("Ele");
         part->setCharge(ele->charge());

         //make a new LorentzVector with defined energy, and direction, mass 0
         TLorentzVector temp;
         temp.SetE( ele->caloEnergy() );
         temp.SetPz( ele->caloEnergy() );
         temp.SetTheta( ele->theta() );
         temp.SetPhi( ele->phi() );

         part->setP4( temp.Px(), temp.Py(),temp.Pz(), temp.E() );


         //reconstruction algorithm was (at least) driven by ECAL
         part->setUserRecord< bool >( "ecalDriven", ele->ecalDrivenSeed() );

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

         // the super cluster energy corrected by EnergyScaleFactor
         part->setUserRecord<float>("SCE", ele->caloEnergy()); //ok
         part->setUserRecord< double >( "SCeta", ele->caloPosition().eta() );
         // the errors on the supercluster energy
         part->setUserRecord<float>("SCEErr", ele->ecalEnergyError()); //ok

         // the errors on the track momentum
         part->setUserRecord<float>("PErr", ele->trackMomentumError());	 
         part->setUserRecord<double>("TrackerP", ele->gsfTrack()->p()); //ok
         //the seed cluster energy / track momentum at calo from outermost state
         part->setUserRecord<double>("ESCSeedPout", ele->eSeedClusterOverPout()); //ok
         //part->setUserRecord<double>("NormChi2", ele->gsfTrack()->normalizedChi2()); //why was this left out?
         part->setUserRecord<int>("TrackerVHits", ele->gsfTrack()->numberOfValidHits()); //ok
         part->setUserRecord<int>("TrackerLHits", ele->gsfTrack()->numberOfLostHits()); //ok
         part->setUserRecord<int>("Class", ele->classification()); //ok

         // Save distance to the primary vertex in z and xy plane i.e. impact parameter
         //get primary vertex (hopefully correct one)
         // units given in cm!!! Use nominal beam spot of Summer/Fall08 if no vertex available
         math::XYZPoint vtx(0.0322, 0., 0.);
         if (RecView->findUserRecord<int>("NumVertices") > 0) {
            pxl::ObjectOwnerTypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
            vtx = math::XYZPoint((*vtx_iter)->getX(), (*vtx_iter)->getY(), (*vtx_iter)->getZ());
         } 
         part->setUserRecord<double>("Dsz", ele->gsfTrack()->dsz(vtx));
         part->setUserRecord<double>("Dxy", ele->gsfTrack()->dxy(vtx));

         //store PAT matching info if MC
         if( MC ){
            std::map<const Candidate*, pxl::Particle*>::const_iterator it = genmap.find( ele->genLepton() );
            if( it != genmap.end() ){
               part->linkSoft( it->second, "pat-match" );
            }
         }

         // Get the supercluster (ref) of the Electron
         // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
         // a SuperClusterCollection is a std::vector<SuperCluster>
         // although we get a vector of SuperClusters an electron is only made out of ONE SC
         // therefore only the first element of the vector should be available!
         const SuperClusterRef SCRef = ele->superCluster();
         const CaloClusterPtr& SCSeed = SCRef->seed();

         //use EcalClusterLazyTools to store ClusterShapeVariables
         part->setUserRecord< double >( "e1x5",  ele->e1x5() );
         part->setUserRecord< double >( "e2x5",  ele->e2x5Max() );
         part->setUserRecord< double >( "e5x5",  ele->e5x5() );

         part->setUserRecord< double >( "SwissCross", EcalSeverityLevelAlgo::swissCross( SCSeed->seed(), *barrelRecHits, 0, false ) );
         EcalRecHitCollection::const_iterator recHit_it = barrelRecHits->find( SCSeed->seed() );
         if( recHit_it != barrelRecHits->end() ) {
            const EcalRecHit &seedRecHit = *recHit_it;
            unsigned int recoFlag = seedRecHit.recoFlag();
            part->setUserRecord< unsigned int >("recoFlag",  recoFlag  );
         }

         std::vector<float> covariances = lazyTools.covariances(*SCSeed, 4.7 );
         part->setUserRecord<double>("EtaEta", covariances[0] ); //used for CutBasedElectronID
         part->setUserRecord<double>("EtaPhi", covariances[1] );
         part->setUserRecord<double>("PhiPhi", covariances[2] );

         part->setUserRecord< double >( "iEta_iEta", ele->scSigmaIEtaIEta() );

         //save eta/phi and DetId info from seed-cluster to prevent duplication of Electron/Photon-Candidates (in final selection)
         part->setUserRecord<double>("seedphi", SCRef->seed()->phi());
         part->setUserRecord<double>("seedeta", SCRef->seed()->eta());
         part->setUserRecord<unsigned int>("seedId", lazyTools.getMaximum(*SCSeed).first.rawId()); 
         //additional data used for cutBasedElectronId in CMSSW 2_0_X
         part->setUserRecord<double>("eSeed", SCRef->seed()->energy()); //used for CutBasedElectronID
         part->setUserRecord<double>("pin", ele->trackMomentumAtVtx().R() ); //used for CutBasedElectronID	 
         part->setUserRecord<double>( "pout", ele->trackMomentumOut().R() ); //used for CutBasedElectronID	
         //store ID information
         const vector< pair< string, float > > &electronIDs = ele->electronIDs();
         for( vector< pair< string, float > >::const_iterator electronID = electronIDs.begin(); electronID != electronIDs.end(); ++electronID ){
            part->setUserRecord<bool>( electronID->first, electronID->second > 0.5 );
         }
         //save official isolation information
         part->setUserRecord<double>("CaloIso", ele->caloIso());
         part->setUserRecord<double>("TrkIso", ele->trackIso());
         part->setUserRecord<double>("ECALIso", ele->ecalIso());
         part->setUserRecord<double>("HCALIso", ele->hcalIso());

         part->setUserRecord< double >( "TrkIso03", ele->dr03TkSumPt() );
         part->setUserRecord< double >( "ECALIso03", ele->dr03EcalRecHitSumEt() );
         part->setUserRecord< double >( "HCALIso03d1", ele->dr03HcalDepth1TowerSumEt() );
         part->setUserRecord< double >( "HCALIso03d2", ele->dr03HcalDepth2TowerSumEt() );
 
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

void MUSiCSkimmer::analyzeRecJets( const edm::Event &iEvent, pxl::EventView *RecView, bool &MC, std::map< const Candidate*, pxl::Particle* > &genjetmap, const jet_def &jet_info){
   //get primary vertex (hopefully correct one) for physics eta
   double VertexZ = 0.;
   if (RecView->findUserRecord<int>("NumVertices") > 0) {
      pxl::ObjectOwnerTypeIterator<pxl::Vertex> vtx_iter = RecView->getObjectOwner().begin<pxl::Vertex>();
      VertexZ = (*vtx_iter)->getZ();
   } 
   
   int numJetRec = 0;
   // get RecoJets
   edm::Handle< std::vector< pat::Jet > > jetHandle;
   iEvent.getByLabel( jet_info.RecoLabel, jetHandle );
   const std::vector< pat::Jet > &RecJets = *jetHandle;

   //generator flavour matching only available in MC. Surprise!
   edm::Handle< reco::JetFlavourMatchingCollection > physicsFlavour;
   if( MC ){
      iEvent.getByLabel( jet_info.name+"RecoJetFlavourPhysics", physicsFlavour );
   }


   // loop over the jets
   size_t jet_index = 0;
   for( std::vector< pat::Jet >::const_iterator jet = RecJets.begin(); jet != RecJets.end(); ++jet, ++jet_index ){
      if (Jet_cuts(jet)) {
         pxl::Particle* part = RecView->create<pxl::Particle>();
         part->setName( jet_info.name );
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

         //calculate the kinematics with a new vertex
         reco::Candidate::LorentzVector physP4 = reco::Jet::physicsP4( the_vertex, *jet, jet->vertex() );
         part->setUserRecord<double>("PhysEta", physP4.eta());
         part->setUserRecord<double>("PhysPhi", physP4.phi());
         part->setUserRecord<double>("PhysPt",  physP4.pt());

         part->setUserRecord< double >( "fHPD", jet->jetID().fHPD );
         part->setUserRecord< double >( "fRBX", jet->jetID().fRBX );

         // store b-tag discriminator values:
         const vector< pair< string, float > > &btags = jet->getPairDiscri();
         for( vector< pair< string, float > >::const_iterator btag = btags.begin(); btag != btags.end(); ++btag ){
            part->setUserRecord< float >( btag->first, btag->second );
         }
         //jet IDs
         for( jet_id_list::iterator ID = jet_ids.begin(); ID != jet_ids.end(); ++ID ){
            pat::strbitset ret = ID->second.getBitTemplate();
            part->setUserRecord< bool >( ID->first, ID->second( *jet, ret ) );
         }
         if (fDebug > 1) {
            const std::vector<std::pair<std::string, float> > & bTags = jet->getPairDiscri();
            part->print(0);
            for (vector<pair<string, float> >::const_iterator btag = bTags.begin(); btag != bTags.end(); ++btag) {
               cout << "Name: " << btag->first << "   value: " << btag->second << endl;
            }
         }
         //store PAT matching info if MC
         if (MC) {
            // to be compared with Generator Flavor:
            part->setUserRecord< int >( "algoFlavour", jet->partonFlavour() );
            //make a ref, then get the flavour
            RefToBase< reco::Jet > jetRef( RefToBaseProd< reco::Jet >( jetHandle ), jet_index );
            part->setUserRecord< int >( "physicsFlavour", (*physicsFlavour)[ jetRef ].getFlavour() );

            std::map< const Candidate*, pxl::Particle* >::const_iterator it = genjetmap.find( jet->genJet() );
            if( it != genjetmap.end() ){
               part->linkSoft( it->second, "pat-match" );
            }
         }
         numJetRec++;
      }
   }
   RecView->setUserRecord< int >( "Num"+jet_info.name, numJetRec );
   if( fDebug > 1 ) cout << "Found Rec Jets:  " << numJetRec << " of Type " << jet_info.name << endl;
}

// ------------ reading Reconstructed Gammas ------------

void MUSiCSkimmer::analyzeRecGammas(const edm::Event& iEvent, pxl::EventView* RecView, bool& MC, EcalClusterLazyTools& lazyTools, std::map<const Candidate*, pxl::Particle*> & genmap){
   // get Photon Collection     
   edm::Handle<std::vector<pat::Photon> > photonHandle;
   iEvent.getByLabel(fGammaRecoLabel, photonHandle);
   const std::vector<pat::Photon>& photons = *photonHandle;

   edm::Handle< EcalRecHitCollection > barrelRecHits;
   iEvent.getByLabel( freducedBarrelRecHitCollection, barrelRecHits );
   
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
         const CaloClusterPtr& SCSeed = SCRef->seed();
         part->setUserRecord<double>("rawEnergy",  SCRef->rawEnergy() );
         part->setUserRecord<double>("preshowerEnergy",  SCRef->preshowerEnergy() );
         //use EcalClusterLazyTools to store ClusterShapeVariables
         double e3x3 = lazyTools.e3x3( *SCSeed );
         part->setUserRecord<double>("e3x3",  e3x3 );
         part->setUserRecord<double>("e5x5",  lazyTools.e5x5(*SCSeed)  );
         part->setUserRecord<double>( "SwissCross", EcalSeverityLevelAlgo::swissCross( SCSeed->seed(), *barrelRecHits, 0, false ) );
         part->setUserRecord<double>( "SwissCrossNoBorder", EcalSeverityLevelAlgo::swissCross( SCSeed->seed(), *barrelRecHits, 0, true ) );
         EcalRecHitCollection::const_iterator recHit_it = barrelRecHits->find( SCSeed->seed() );
         if( recHit_it != barrelRecHits->end() ) {
            const EcalRecHit &seedRecHit = *recHit_it;
            unsigned int recoFlag = seedRecHit.recoFlag();
            part->setUserRecord< unsigned int >("recoFlag",  recoFlag  );
         }
         std::vector<float> covariances = lazyTools.covariances(*SCSeed );
         part->setUserRecord<double>("EtaEta", covariances[0] ); 
         part->setUserRecord<double>("EtaPhi", covariances[1] );
         part->setUserRecord<double>("PhiPhi", covariances[2] );
         part->setUserRecord<double>("Emax",  lazyTools.eMax(*SCSeed)  );
         part->setUserRecord<double>("r9", e3x3 /( SCRef->rawEnergy() + SCRef->preshowerEnergy() ) );
         // part->setUserRecord<double>("r9", photon->r9()); <== different computation of r9 here :-(
         part->setUserRecord<double>("r19",  (lazyTools.eMax(*SCSeed) / e3x3 ) );
         //save eta/phi and DetId info from seed-cluster to prevent dublication of Electron/Photon-Candidates (in final selection) adn to reject converted photons
         part->setUserRecord<double>("seedphi", SCRef->seed()->phi());
         part->setUserRecord<double>("seedeta", SCRef->seed()->eta());
         part->setUserRecord<unsigned int>("seedId", lazyTools.getMaximum(*SCSeed).first.rawId()); 
         //set hadronic over electromagnetic energy fraction
         part->setUserRecord<float>("HoEm", photon->hadronicOverEm());
         //save official isolation information
         // this is the BAD PAT isolation!!!
         part->setUserRecord<double>("HCALIso", photon->hcalIso());
         part->setUserRecord<double>("ECALIso", photon->ecalIso());
         part->setUserRecord<double>("TrkIso", photon->trackIso());
         part->setUserRecord<double>("CaloIso", photon->caloIso());
         part->setUserRecord<int>("TrackNum", photon->nTrkSolidConeDR04());
         // use egamma isolation based on RecHits:
         part->setUserRecord<float>("ID_HCALIso", photon->hcalTowerSumEtConeDR04());
         part->setUserRecord<float>("ID_ECALIso", photon->ecalRecHitSumEtConeDR04());
         part->setUserRecord<float>("ID_TrkIso", photon->trkSumPtHollowConeDR04());	  
         //store information about converted state
         part->setUserRecord<bool>("Converted", photon->hasConversionTracks());
         //if (photon->isConvertedPhoton() == true) {cout << "is Converted!" << endl;}
         // store photon-Id
         const vector< pair< string, bool > > &photonIDs = photon->photonIDs();
         for( vector< pair< string, bool > >::const_iterator photonID = photonIDs.begin(); photonID != photonIDs.end(); ++photonID ){
            part->setUserRecord<bool>( photonID->first, photonID->second );
         }
         // is near a gap ! isEEGap == always false not yet implemented ... CMSSW_2_1_9 
         part->setUserRecord<bool>("Gap", photon->isEBGap() || photon->isEEGap() || photon->isEBEEGap());

         // store Gamma info corrected for primary vertex (this changes direction but leaves energy of SC unchanged 
         pat::Photon localPho(*photon);
         // Set event vertex
         localPho.setVertex( the_vertex );
         part->setUserRecord<double>("PhysEta", localPho.eta());
         part->setUserRecord<double>("PhysPhi", localPho.phi());
         part->setUserRecord<double>("PhysPt", localPho.pt());

         //store PAT matching info
         if (MC) {
            std::map< const Candidate*, pxl::Particle* >::const_iterator it = genmap.find( photon->genPhoton() );
            if( it != genmap.end() ){
               part->linkSoft( it->second, "pat-match" );
            }
         }
         //FIXME Pi0 stuff still missing --> seems not be working in CMSSW_2_1_X
         numGammaRec++;
      }	 
   }
   RecView->setUserRecord<int>("NumGamma", numGammaRec);
   if (fDebug > 1) cout << "Rec Gamma: " << numGammaRec << endl;
}





void MUSiCSkimmer::analyzeECALRecHits( const edm::Event &iEvent, const edm::EventSetup &iSetup, pxl::EventView *RecView ) {
   //get the ECAL barrel rec hits
   edm::Handle< EcalRecHitCollection > barrelRecHits_h;
   iEvent.getByLabel( freducedBarrelRecHitCollection, barrelRecHits_h );
   const EcalRecHitCollection &barrelRecHits = *barrelRecHits_h;

   //get the calo geometry
   edm::ESHandle< CaloGeometry > geo;
   iSetup.get< CaloGeometryRecord >().get( geo );

   edm::ESHandle<CaloTopology> topo_h;
   iSetup.get< CaloTopologyRecord >().get( topo_h );
   const CaloTopology &topo = *topo_h;

   //we need a totally useless dummy cluster for matrixEnergy()
   reco::BasicCluster dummy;

   //loop over all rec-hits and store the 4-vector of weird ones
   for( EcalRecHitCollection::const_iterator rechit_it = barrelRecHits.begin(); rechit_it != barrelRecHits.end(); ++rechit_it ) {
      const EcalRecHit &rechit = *rechit_it;
      double energy = rechit.energy();
      if( energy > min_rechit_energy ) {
         const EBDetId id( rechit.id() );
         uint32_t recoFlag = rechit.recoFlag();
         double swiss_cross = EcalSeverityLevelAlgo::swissCross( id, barrelRecHits, 0, false );
         double e3x3 = EcalClusterTools::matrixEnergy( dummy, &barrelRecHits, &topo, id, -1, 1, -1, 1 );
         double r19 = energy / e3x3;

         if( recoFlag == EcalRecHit::kOutOfTime
             || swiss_cross > min_rechit_swiss_cross
             || r19 > min_rechit_R19
             ) {
            int ieta = id.ieta();

            pxl::Particle* part = RecView->create<pxl::Particle>();
            part->setName("ECALRecHit");
            part->setCharge(0);

            TLorentzVector temp;
            temp.SetE( energy );
            temp.SetPz( energy );
            temp.SetTheta( geo->getPosition( id ).theta() );
            temp.SetPhi( geo->getPosition( id ).phi() );
            
            part->setP4( temp.Px(), temp.Py(),temp.Pz(), temp.E() );
         
            part->setUserRecord< int >( "ieta", ieta );
            part->setUserRecord< uint32_t >( "recoFlag", recoFlag );
            part->setUserRecord< double >( "SwissCross", swiss_cross );
            part->setUserRecord< double >( "r19", r19 );
         }
      }
   }
}





// ------------ method returning the EventClassType ------------

std::string MUSiCSkimmer::getEventClass(pxl::EventView* EvtView) {
   ostringstream EventType;
   //set default values to 0 for Gen-only mode
   EventType << EvtView->findUserRecord<int>("NumEle") <<  "e"
             << EvtView->findUserRecord<int>("NumMuon") << "mu"
             << EvtView->findUserRecord<int>("NumGamma") << "gam"
             << EvtView->findUserRecord<int>("Num"+jet_infos[0].name) << "jet"
             << EvtView->findUserRecord<int>("NumMET") << "met";
   EventType.flush();
   return EventType.str();
}

// ------------ method called once each job just after ending the event loop  ------------

void MUSiCSkimmer::endJob() {
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

      //load the best fit PDF
      int first_pdf = 0; //stupid c++
      initpdf_( first_pdf );
      for( vector< gen::PdfInfo >::const_iterator pdf = fpdf_vec.begin(); pdf != fpdf_vec.end(); ++pdf ){
         if( pdf->scalePDF != 0 ){
            if( abs( pdf->id.first ) <= 6 && abs( pdf->id.second ) <= 6 ){
               best_fit.push_back( xfx( pdf->x.first, pdf->scalePDF, pdf->id.first ) * xfx( pdf->x.second, pdf->scalePDF, pdf->id.second ) );
            } else {
               best_fit.push_back( 1 );
               throw cms::Exception( "PDF error" ) << "Found an event with weird partons!" << endl
                                                   << "This should not happen! (Error should have been caught before!)" << endl
                                                   << "Details:" << endl
                                                   << "x1: " << pdf->x.first << endl
                                                   << "x2: " << pdf->x.second << endl
                                                   << "Scale: " << pdf->scalePDF << endl
                                                   << "ID 1: " << pdf->id.first << endl
                                                   << "ID 1: " << pdf->id.second << endl;
            }
         } else {
            best_fit.push_back( 1 );
         }
      }

      //loop over all error PDFs
      for( int subpdf = 1; subpdf < 41; subpdf++ ){
         initpdf_(subpdf);
         //cout << "Initialized sub PDF " << subpdf << endl;
         vector<float>::const_iterator best_fit_iter = best_fit.begin();
         vector<vector<float> >::iterator weights_iter = weights.begin();
         // loop over all PDFInf's
         for( vector< gen::PdfInfo >::const_iterator pdf = fpdf_vec.begin(); pdf != fpdf_vec.end(); ++pdf ){
            if( pdf->scalePDF != 0 && abs( pdf->id.first ) <= 6 && abs( pdf->id.second ) <= 6 ){
               weights_iter->push_back( xfx( pdf->x.first, pdf->scalePDF, pdf->id.first) * xfx( pdf->x.second, pdf->scalePDF, pdf->id.second ) / (*best_fit_iter));
            } else {
               weights_iter->push_back( 1 );
            }
            ++weights_iter;
            ++best_fit_iter;
         }
      }
      // ReRead the pxlio file and store PDFInfo
      pxl::InputFile Input(fFileName);
      pxl::OutputFile tmpFile("Tmp"+fFileName);
      vector<vector<float> >::const_iterator weights_iter = weights.begin();
      int count = 1;
      // run event loop:
      while (Input.nextEvent()) {
         pxl::Event event;
         // read event from disk
         Input.readEvent(&event);
         // get all stored EventViews
         pxl::EventView* GenEvtView = event.getObjectOwner().findObject<pxl::EventView>("Gen");
         pxl::EventView* RecEvtView = event.getObjectOwner().findObject<pxl::EventView>("Rec");
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
         tmpFile.writeEvent(&event);
         ++weights_iter;
         ++count;
      }
      tmpFile.close();
      Input.close();
      // rename tmporary file
      system(("mv Tmp" + fFileName + " " + fFileName).c_str());
   }

   //write a single EOF byte at the end of the file
   //that doesn't hurt PXL, but should avoid the "file has zero size" stage-out problem
   system( ("echo -e \\0004 >> "+fFileName).c_str() );
}


// ------------ method to define MC-MUON-cuts

bool MUSiCSkimmer::MuonMC_cuts(const GenParticle* MCmuon) const {
   if( MCmuon->pt() < min_muon_pt ) return false;
   if( fabs( MCmuon->eta() ) > max_eta ) return false;
   return true;
}
 


// ------------ method to define MC-Electron-cuts

bool MUSiCSkimmer::EleMC_cuts(const GenParticle* MCele) const {
   if( MCele->pt() < min_ele_pt ) return false;
   if( fabs( MCele->eta() ) > max_eta ) return false;
   return true;
}

// ------------ method to define MC-Gamma-cuts

bool MUSiCSkimmer::GammaMC_cuts(const GenParticle* MCgamma) const {
   if( MCgamma->pt() < min_gamma_pt) return false;
   if( fabs(MCgamma->eta() ) > max_eta ) return false;
   return true;
}

// ------------ method to define MC-Jet-cuts

bool MUSiCSkimmer::JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const {
   if( MCjet->pt() < min_jet_pt ) return false;
   if( fabs( MCjet->eta() ) > max_eta ) return false;
   return true;
}

// ------------ method to define MC-MET-cuts

bool MUSiCSkimmer::METMC_cuts(const pxl::Particle* MCmet) const {
   if( MCmet->getPt() < min_met ) return false;
   return true; 
}

// ------------ method to define RecVertex-cuts
bool MUSiCSkimmer::Vertex_cuts( reco::VertexCollection::const_iterator vertex ) const {
   return ( vertex->ndof() >= vertex_minNDOF
            && fabs( vertex->z() ) <= vertex_maxZ
            && vertex->position().rho() <= vertex_maxR );
}

bool MUSiCSkimmer::PV_vertex_cuts( const reco::Vertex &vertex ) const {
   return ( vertex.ndof() >= PV_minNDOF
            && fabs( vertex.z() ) <= PV_maxZ
            && vertex.position().rho() <= PV_maxR );
}

// ------------ method to define MUON-cuts

bool MUSiCSkimmer::Muon_cuts(const pat::Muon& muon) const {
   // basic preselection cuts
   if( !muon.isGlobalMuon() ) return false;
   if( muon.pt() < min_muon_pt )  return false;
   if( fabs( muon.eta() ) > max_eta ) return false;
   return true;
}


// ------------ method to define ELECTRON-cuts

bool MUSiCSkimmer::Ele_cuts(std::vector<pat::Electron>::const_iterator ele) const {
   if( ele->pt() < min_ele_pt ) return false;
   if( fabs( ele->eta() ) > max_eta ) return false;
   return true;
}

// ------------ method to define JET-cuts

bool MUSiCSkimmer::Jet_cuts(std::vector<pat::Jet>::const_iterator jet) const {
   if( jet->pt() < min_jet_pt ) return false;
   if( fabs( jet->eta()) > max_eta ) return false;
   return true;
}


// ------------ method to define GAMMA-cuts

bool MUSiCSkimmer::Gamma_cuts(std::vector<pat::Photon>::const_iterator photon) const {
   if( photon->pt() < min_gamma_pt ) return false;
   if( fabs( photon->eta() ) > max_eta ) return false;
   return true;
}


// ------------ method to define MET-cuts

bool MUSiCSkimmer::MET_cuts(const pxl::Particle* met) const {
   if( met->getPt() < min_met ) return false;
   return true;
}

//------------------------------------------------------------------------------

//FIXME compare to PAT-isolation 
double MUSiCSkimmer::IsoGenSum (const edm::Event& iEvent, double ParticleGenPt, double ParticleGenEta, double ParticleGenPhi, double iso_DR, double iso_Seed){
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
DEFINE_FWK_MODULE(MUSiCSkimmer);
