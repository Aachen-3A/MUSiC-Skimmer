// -*- C++ -*-
//
// Package:    MUSiCProject
// Class:      MUSiCSkimmer
//
//\class MUSiCSkimmer MUSiCProject/Skimming/src/MUSiCSkimmer.cc
//
//Description: Data and MC Skimmer for the Model Unspecific Search in CMS
//
//Implementation:
//
// Original Authors: Carsten Hof, Philipp Biallass, Holger Pieta, Paul Papacz
//         Created:  Mo Okt 30 12:03:52 CET 2006
// $Id$
//
//
// Own header file.
#include "MUSiCProject/Skimming/interface/MUSiCSkimmer.h"

// System include files.
#include <iostream>

// Message Logger for Debug etc.
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Exceptions. Do *not* use edm::LogError(), use cms::Exception() instead!
#include "FWCore/Utilities/interface/Exception.h"

// Necessary objects.
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// For GenParticles.
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// Special stuff for sim truth of converted photons.
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

// Pile-Up information.
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// PDF stuff.
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"

// PAT related stuff.
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

// EGamma stuff.
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperClusterFwd.h"

// For Muon stuff.
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

// Jet stuff.
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

// ECAL + HCAL
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"

// For Trigger Bits:
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

// Misc.
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// Math stuff from Physics tools.
#include "DataFormats/Math/interface/deltaR.h"

// Private ParticleMatcher.
#include "MUSiCProject/Skimming/interface/ParticleMatcher.hh"

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
   // Use SIM info
   fUseSIM = iConfig.getUntrackedParameter<bool>("UseSIM");
   // name of the LHgrid for pdf weights
   fLHgridName = iConfig.getUntrackedParameter<string>("LHgridName"),
   // number of pdf error sets in the LHgrid for pdf weights
   fNumLHgridErrorSets = iConfig.getUntrackedParameter<int>("NumLHgridErrorSets"),
   // The labels used in cfg-file 
   //fTruthVertexLabel = iConfig.getUntrackedParameter<string>("TruthVertexLabel");
   fgenParticleCandidatesLabel  = iConfig.getUntrackedParameter<string>("genParticleCandidatesLabel");
   fVertexRecoLabel = iConfig.getUntrackedParameter<string>("VertexRecoLabel");
   fTauRecoLabel = iConfig.getUntrackedParameter< string >( "TauRecoLabel" );
   fMuonRecoLabel = iConfig.getUntrackedParameter<string>("MuonRecoLabel");
   fElectronRecoLabel = iConfig.getUntrackedParameter<string>("ElectronRecoLabel");
   m_gsfElectronsTag = iConfig.getParameter< InputTag >( "gsfElectronsTag" );
   m_inputTagIsoValElectronsPFId = iConfig.getParameter< vector< InputTag > >( "IsoValElectronPF" );
   m_eleEffAreaTargetLabel = iConfig.getUntrackedParameter< string >( "EleEffAreaTargetLabel" );

   freducedBarrelRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
   freducedEndcapRecHitCollection = iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");
   fGammaRecoLabel = iConfig.getUntrackedParameter<string>("GammaRecoLabel");
   m_inputTagIsoValPhotonsPFId = iConfig.getParameter< vector< InputTag > >( "IsoValPhotonPF" );

   m_particleFlowTag = iConfig.getParameter< InputTag >( "particleFlowTag" );

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
      vector< string > id_names = one_jet.getParameter< vector< string > >( "IDs" );
      jet.isPF = one_jet.getParameter< bool >( "isPF" );

      //JetIDs
      unsigned int num_IDs=0;
      if (jet.isPF) num_IDs = PFJetIDSelectionFunctor::N_QUALITY;
      else num_IDs = JetIDSelectionFunctor::N_QUALITY;
      if( num_IDs < id_names.size() ){
         edm::LogWarning( "JetIDs" ) << "Less JetIDs available (" << num_IDs << ") than requested (" << id_names.size() << "), using only available.";
      } else if( num_IDs > id_names.size() ){
         edm::LogWarning( "JetIDs" ) << "More JetIDs available (" << num_IDs << ") than requested (" << id_names.size() << "), using only requested.";
         num_IDs = id_names.size();
      }

      //ATTENTION: The following is REALLY ugly
      //Looping over enums is apparentlt not forseen in C++
      //Seems to be the only way to make the JetIDs configurable
      if (jet.isPF) {
         for( PFJetIDSelectionFunctor::Quality_t q = PFJetIDSelectionFunctor::Quality_t(0);
            q < PFJetIDSelectionFunctor::N_QUALITY;
            q = PFJetIDSelectionFunctor::Quality_t( q+1 ) ){
            pair< std::string, ::Selector<pat::Jet>* > ID( id_names[q], new PFJetIDSelectionFunctor( PFJetIDSelectionFunctor::FIRSTDATA, q ) );
            jet.IDs.push_back( ID );
         }
      } else {
         for( JetIDSelectionFunctor::Quality_t q = JetIDSelectionFunctor::Quality_t(0);
            q < JetIDSelectionFunctor::N_QUALITY;
            q = JetIDSelectionFunctor::Quality_t( q+1 ) ){
            pair< std::string, ::Selector<pat::Jet>* > ID( id_names[q], new JetIDSelectionFunctor( JetIDSelectionFunctor::PURE09, q ) );
            jet.IDs.push_back( ID );
         }
      }

      jet_infos.push_back( jet );
   }



   //MET
   //get the PSet that contains all MET PSets
   ParameterSet METs_pset = iConfig.getParameter< ParameterSet >( "METs" );
   //get the names of all sub-PSets
   vector< string > MET_names;
   METs_pset.getParameterSetNames( MET_names );
   //loop over the names of the MET PSets
   for( vector< string >::const_iterator MET_name = MET_names.begin(); MET_name != MET_names.end(); ++MET_name ){
      collection_def MET;
      MET.name = *MET_name;

      ParameterSet one_MET = METs_pset.getParameter< ParameterSet >( MET.name );
      MET.MCLabel   = one_MET.getParameter< InputTag >( "MCLabel" );
      MET.RecoLabel = one_MET.getParameter< InputTag >( "RecoLabel" );

      MET_infos.push_back( MET );
   }

   //HCAL noise
   hcal_noise_label = iConfig.getParameter< InputTag >( "HCALNoise" );

   // Conversions
   m_conversionsTag = iConfig.getParameter< InputTag >( "conversionsTag" );

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
      if( trigger.process == "auto" ) {
         trigger.results = InputTag( one_trigger.getParameter< string >( "results" ), "" );
         trigger.event   = InputTag( one_trigger.getParameter< string >( "event" ), "" );
      } else {
         trigger.results = InputTag( one_trigger.getParameter< string >( "results" ), "", trigger.process );
         trigger.event   = InputTag( one_trigger.getParameter< string >( "event" ),   "", trigger.process );
      }
      
      trigger.triggers_names = one_trigger.getParameter< vector< string > >("HLTriggers");

      if( not fGenOnly and trigger.triggers_names.size() == 0 ) {
         throw cms::Exception( "Trigger error" ) << "No Trigger names found in configuration! This is only in 'GenOnly' mode allowed." << endl;
      }

      triggers.push_back( trigger );
   }

   fStoreL3Objects = trigger_pset.getUntrackedParameter<bool>("StoreL3Objects");

   // Filters
   // -------
   // This is based on the triggers handling from above because the information
   // from filters that ran are accessed with help of the edm::TriggerResults.
   // Basically it is foreseen (but not used atm.) to use more than one filter combination.
   //
   ParameterSet filter_pset = iConfig.getParameter< ParameterSet >( "filters" );

   vector< string > filter_paths;
   filter_pset.getParameterSetNames( filter_paths );

   for( vector< string >::const_iterator filter_path = filter_paths.begin(); filter_path != filter_paths.end(); ++filter_path ) {
      trigger_group filter;
      filter.name = *filter_path;

      ParameterSet one_filter = filter_pset.getParameter< ParameterSet >( filter.name );
      filter.process = one_filter.getParameter< string >( "process" );

      filter.results = InputTag( one_filter.getParameter< string >( "results" ), "" );

      filter.triggers_names = one_filter.getParameter< vector< string > >( "paths" );

      filters.push_back( filter );
   }


   //cuts
   ParameterSet cut_pset = iConfig.getParameter< ParameterSet >( "cuts" );
   min_tau_pt = cut_pset.getParameter< double >( "min_tau_pt" );
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
   edm::LogInfo( "MUSiCSkimmer|EventInfo" ) << "Run " << iEvent.id().run()
                                            << ", EventID = " << iEvent.id().event()
                                            << ", is MC = " << !iEvent.isRealData();

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

   // create two ePaxEventViews for Generator/Reconstructed Objects
   pxl::EventView* GenEvtView = event.createIndexed<pxl::EventView>("Gen");
   pxl::EventView* RecEvtView = event.createIndexed<pxl::EventView>("Rec");
   GenEvtView->setName("Gen"); RecEvtView->setName("Rec");
   GenEvtView->setUserRecord<std::string>("Type", "Gen");
   RecEvtView->setUserRecord<std::string>("Type", "Rec");
   
   //maps for matching
   std::map< const reco::Candidate*, pxl::Particle* > genmap;
   std::map< const reco::Candidate*, pxl::Particle* > genjetmap;

   //set process name
   GenEvtView->setUserRecord<std::string>("Process", fProcess);
   RecEvtView->setUserRecord<std::string>("Process", fProcess);

   // Generator stuff
   if (IsMC) {
      analyzeGenInfo(iEvent, GenEvtView, genmap);
      analyzeGenRelatedInfo(iEvent, GenEvtView);  // PDFInfo, Process ID, scale, pthat
      for( vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         analyzeGenJets(iEvent, GenEvtView, genjetmap, *jet_info);
      }
      for( vector< collection_def >::const_iterator MET_info = MET_infos.begin(); MET_info != MET_infos.end(); ++MET_info ){
         analyzeGenMET(iEvent, GenEvtView, *MET_info);
      }
      if (fUseSIM){
         analyzeSIM(iEvent, GenEvtView);
      }
   }
   // store Rec Objects only if requested
   if (!fGenOnly) {
      // Median pt per area for each event.
      // See also:
      // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Rho_for_2011_Effective_Areas
      // https://twiki.cern.ch/twiki/bin/view/CMS/Vgamma2011PhotonID#Recommended_cuts
      //
      edm::Handle< double > rho25;
      iEvent.getByLabel( "kt6PFJets", "rho", rho25 );
      m_rhoFastJet = *rho25;

      RecEvtView->setUserRecord< double >( "rho25", m_rhoFastJet );

      //get the calo geometry
      edm::ESHandle< CaloGeometry > geo;
      iSetup.get< CaloGeometryRecord >().get( geo );

      //create object for EcalClusterLazyTools
      EcalClusterLazyTools lazyTools( iEvent, iSetup, freducedBarrelRecHitCollection, freducedEndcapRecHitCollection);

      //Trigger bits
      for( vector< trigger_group >::iterator trg = triggers.begin(); trg != triggers.end(); ++trg ){
         analyzeTrigger( iEvent, iSetup, RecEvtView, *trg );
      }

      //Filters more info see above.
      //
      for( vector< trigger_group >::iterator filt = filters.begin(); filt != filters.end(); ++filt ){
         analyzeFilter( iEvent, iSetup, RecEvtView, *filt );
      }

      // Initalise the helper for 2012 H/E definition and HCAL isolation.
      // Somehow this is quite ugly.
      //
      ElectronHcalHelper::Configuration hcalCfg;

      hcalCfg.hOverEConeSize = 0.15;
      hcalCfg.useTowers      = true;
      hcalCfg.hcalTowers     = edm::InputTag( "towerMaker" );
      hcalCfg.hOverEPtMin    = 0;

      m_hcalHelper = new ElectronHcalHelper( hcalCfg );
      m_hcalHelper->checkSetup( iSetup );
      m_hcalHelper->readEvent( iEvent );

      // Reconstructed stuff
      analyzeRecVertices(iEvent, RecEvtView);
      analyzeRecTaus( iEvent, RecEvtView, IsMC, genmap );
      analyzeRecMuons(iEvent, RecEvtView, IsMC, genmap);
      analyzeRecElectrons( iEvent, RecEvtView, IsMC, lazyTools, genmap, geo );
      for( vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         analyzeRecJets( iEvent, RecEvtView, IsMC, genjetmap, *jet_info );
      }
      for( vector< collection_def >::const_iterator MET_info = MET_infos.begin(); MET_info != MET_infos.end(); ++MET_info ){
         analyzeRecMET( iEvent, RecEvtView, *MET_info );
      }
      analyzeHCALNoise(iEvent, RecEvtView);
      analyzeRecGammas( iEvent, RecEvtView, IsMC, lazyTools, genmap, geo );
   }

   if (IsMC && !fGenOnly){
      const string met_name = "MET";
      Matcher->matchObjects(GenEvtView, RecEvtView, jet_infos, met_name);
   }

   if( not fGenOnly ) {
      stringstream info;
      info << "UserRecord     e     mu    Gamma ";
      for( std::vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         info << jet_info->name << " ";
      }
      info << "  MET  " << endl;
      if( IsMC ) {
         info << "Found (MC): " << setw( 4 ) << GenEvtView->findUserRecord< int >( "NumEle" )
                                << setw( 7 ) << GenEvtView->findUserRecord< int >( "NumMuon" )
                                << setw( 7 ) << GenEvtView->findUserRecord< int >( "NumGamma" );
         for( std::vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
            info << setw( 4 ) << GenEvtView->findUserRecord< int >( "Num" + jet_info->name ) << " ";
         }
         info << setw( 7 ) << GenEvtView->findUserRecord< int >( "NumMET" ) << endl;
      }
      info << "     (Rec): " << setw( 4 ) << RecEvtView->findUserRecord< int >( "NumEle" )
                             << setw( 7 ) << RecEvtView->findUserRecord< int >( "NumMuon" )
                             << setw( 7 ) << RecEvtView->findUserRecord< int >( "NumGamma" );
      for( std::vector< jet_def >::const_iterator jet_info = jet_infos.begin(); jet_info != jet_infos.end(); ++jet_info ){
         info << setw( 4 ) << RecEvtView->findUserRecord< int >( "Num" + jet_info->name ) << " ";
      } 
      info << setw( 7 ) << RecEvtView->findUserRecord< int >( "NumMET" ) << endl;
      edm::LogVerbatim( "MUSiCSkimmer|EventInfo" ) << info.str();
   }   
   fePaxFile.writeEvent(&event);
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

   stringstream info;
   info << "Event Scale (i.e. pthat) = "
        << ( genEvtInfo->hasBinningValues() ? genEvtInfo->binningValues()[0] : 0 )
        << ", EventWeight = " << genEvtInfo->weight() << endl;
   if( genEvtInfo->hasPDF() ){
      info << "PDFInfo: " << endl
           << "========" << endl;
      info << "Momentum of first incoming parton: (id/flavour = "
           << genEvtInfo->pdf()->id.first  << ") "
           << genEvtInfo->pdf()->x.first  << endl
           << "Momentum of second incoming parton: (id/flavour = "
           << genEvtInfo->pdf()->id.second << ") "
           << genEvtInfo->pdf()->x.second << endl
           << "Scale = " << genEvtInfo->pdf()->scalePDF << endl;
   } else {
      info << "No PDFInfo in this event." << endl;
   }
   edm::LogVerbatim( "MUSiCSkimmer|PDFInfo" ) << info.str();
}

// ------------ reading the Generator Stuff ------------

void MUSiCSkimmer::analyzeGenInfo( const edm::Event& iEvent, pxl::EventView* EvtView, std::map< const reco::Candidate*, pxl::Particle* >& genmap ) {
   //gen particles
   edm::Handle<reco::GenParticleCollection> genParticleHandel;
   iEvent.getByLabel(fgenParticleCandidatesLabel , genParticleHandel );

     
   const reco::GenParticle* p = (const reco::GenParticle*) &(*genParticleHandel->begin()); //this is the incoming proton
   pxl::Vertex* GenVtx = EvtView->create<pxl::Vertex>();
   GenVtx->setName("PV");
   //mind that clearly the following line crashes in case of ParticleGun RelVal like single photon
   // therefore add a protection :-) C.H. 20.04.09
   if (p->daughter(0) != 0) GenVtx->setXYZ(p->daughter(0)->vx(), p->daughter(0)->vy(), p->daughter(0)->vz()); //need daughter since first particle (proton) has position zero
   else GenVtx->setXYZ(p->vx(), p->vy(), p->vz());  // if we do not have pp collisions
   EvtView->setUserRecord<int>("NumVertices", 1);  
  
   int numTauMC = 0;
   int numMuonMC = 0;
   int numEleMC = 0;
   int numGammaMC = 0;
   int GenId = 0;

   //save mother of stable particle
   const reco::Candidate* p_mother;
   int mother = 0;
   
   // loop over all particles
   for (reco::GenParticleCollection::const_iterator pa = genParticleHandel->begin(); pa != genParticleHandel->end(); ++ pa ) {
      //cast iterator into GenParticleCandidate
      const reco::GenParticle* p = (const reco::GenParticle*) &(*pa);

      // the following is interesting for GEN studies
      p_mother = p->mother();
      if( p_mother ) {
         if( p->status() == 3 ) {
            int p_id = p->pdgId();
            mother = p_mother->pdgId();
            if( p_id != mother && mother <= 100 ) {
               pxl::Particle* part = EvtView->create< pxl::Particle >();
               part->setName( "S3" );
               part->setP4( p->px(), p->py(), p->pz(), p->energy() );
               part->setUserRecord< int >( "id", p_id );
               part->setUserRecord< int >( "mother_id", mother );
            }
         }
      }

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

      // fill Gen Taus passing some basic cuts -> status? What kind of taus can be found?
      if( abs( p->pdgId() ) == 15 && p->status() == 2 && TauMC_cuts( p ) ) {
         //check whether the tau is final or radiates a tau
         bool isfinal = true;
         for( reco::GenParticle::const_iterator daughter = p->begin(); daughter != p->end(); ++daughter ) {
            if( daughter->pdgId() == 15 ) {
               isfinal = false;
               break;
            }
         }
         if( isfinal ) {
            pxl::Particle *part = EvtView->create< pxl::Particle >();
            genmap[ p ] = part; //fill genmap
            part->setName( "Tau" );
            part->setCharge( p->charge() );
            part->setP4( p->px(), p->py(), p->pz(), p->energy() );
            part->setUserRecord< float >( "Vtx_X", p->vx() );
            part->setUserRecord< float >( "Vtx_Y", p->vy() );
            part->setUserRecord< float >( "Vtx_Z", p->vz() );
            part->setUserRecord< int >( "GenId", GenId );
            numTauMC++;
         }
      }
   } //end of loop over generated particles

   edm::LogInfo( "MUSiCSkimmer|GenInfo" ) << "MC Found: " << numMuonMC <<  " muon(s), " << numEleMC << " electron(s), " << numGammaMC << " gamma(s)";
   EvtView->setUserRecord<int>("NumMuon", numMuonMC);
   EvtView->setUserRecord<int>("NumEle", numEleMC);
   EvtView->setUserRecord<int>("NumGamma", numGammaMC);

   // take care of the pile-up in the event
   //
   Handle< std::vector< PileupSummaryInfo > >  PUInfo;
   iEvent.getByLabel( InputTag( "addPileupInfo" ), PUInfo );

   vector< PileupSummaryInfo >::const_iterator PUiter;

   // loop over all PU info object in an event and get the number of
   // primary vertices for in-time and out-of-time pile-up
   // See also:
   // https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_2011_Reweighting
   //
   for( PUiter = PUInfo->begin(); PUiter != PUInfo->end(); ++PUiter ) {
      int BX      = (*PUiter).getBunchCrossing();
      int num     = (*PUiter).getPU_NumInteractions();

      if( BX == -1 ) {
         EvtView->setUserRecord< int >( "NumVerticesPULastBX", num );
      } else if( BX == 0 ) {
         EvtView->setUserRecord< int >( "NumVerticesPU", num );
         // The true number of interactions (i.e., the mean used in the Poisson
         // distribution) should be the same for in-time and out-of-time
         // pile-up as the actual number is drawn from the same Poisson distribution.
         //
         EvtView->setUserRecord< int >( "NumVerticesPUTrue", (*PUiter).getTrueNumInteractions() );
      } else if( BX == 1 ) {
         EvtView->setUserRecord< int >( "NumVerticesPUNextBX", num );
      }
   }
}

// ------------ reading the Generator Jets ------------

void MUSiCSkimmer::analyzeGenJets( const edm::Event &iEvent, pxl::EventView *EvtView, std::map< const reco::Candidate*, pxl::Particle* > &genjetmap, const jet_def &jet_info ) {
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
         const reco::GenParticle *p = dynamic_cast< const reco::GenParticle* >( &(*genJet) );
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
         const vector< const reco::GenParticle* > &genJetConstit = genJet->getGenConstituents();
         for( std::vector< const reco::GenParticle* >::const_iterator constit = genJetConstit.begin(); constit != genJetConstit.end(); ++constit ) {
            //raise counter if cut passed
            if( (*constit)->pt() > constit_pT ) numGenJetConstit_withcuts++; 
         }
         part->setUserRecord< int >( "GenJetConstit", numGenJetConstit_withcuts );
         part->setUserRecord< int >( "algoFlavour",    (*algoFlavour)   [ jetRef ].getFlavour() );
         part->setUserRecord< int >( "physicsFlavour", (*physicsFlavour)[ jetRef ].getFlavour() );
      }
   }
   EvtView->setUserRecord< int >( "Num"+jet_info.name, numJetMC );
   edm::LogInfo( "MUSiCSkimmer|GenInfo" ) << "MC Found: " << numJetMC << " jet(s) of type: " << jet_info.name;
}

// ------------ reading the Generator MET ------------

void MUSiCSkimmer::analyzeGenMET(const edm::Event& iEvent, pxl::EventView* EvtView, const collection_def &MET_info ) {
   edm::Handle<reco::GenMETCollection> GenMet;
   iEvent.getByLabel(MET_info.MCLabel, GenMet);
   const reco::GenMETCollection *genmetcol = GenMet.product();
   const reco::GenMET genmet = genmetcol->front();  // MET exists only once!

   int numMETMC = 0; //means no MET in event


   pxl::Particle* part = EvtView->create<pxl::Particle>();
   part->setName(MET_info.name);
   part->setP4(genmet.px(), genmet.py(), genmet.pz(), genmet.energy());
   part->setUserRecord<double>("sumEt", genmet.sumEt());
   part->setUserRecord<double>("mEtSig", genmet.mEtSig());
   //fill additional jet-related infos
   part->setUserRecord<double>("EmE", genmet.emEnergy());
   part->setUserRecord<double>("HadE", genmet.hadEnergy());
   part->setUserRecord<double>("InvE", genmet.invisibleEnergy());

   edm::LogInfo( "MUSiCSkimmer|GenInfo" ) << "GenMET before muon corr: Px = " << genmet.px()
                                          << ", Py = " << genmet.py()
                                          << ", Pt = " << part->getPt();

   // Perform Muon Corrections!
   // loop over muons and subtract them
   if (EvtView->findUserRecord<int>("NumMuon") > 0) { 
      vector<pxl::Particle*> GenMuons;
      pxl::ParticleFilter::apply( EvtView->getObjectOwner(), GenMuons, pxl::ParticlePtEtaNameCriterion("Muon") );
      for (vector<pxl::Particle*>::const_iterator muon = GenMuons.begin(); muon != GenMuons.end(); ++muon) {
         edm::LogInfo( "MUSiCSkimmer|GenInfo" ) << "Correcting with " << (*muon)->getName()
                                                << ", Px = " << (*muon)->getPx()
                                                << ", Py = " << (*muon)->getPy();
         *part -= **muon;
      }
      //reset eta-info after muon corrections
      //part->setP4(part->getPx(), part->getPy(), 0., sqrt(part->getPx()*part->getPx()+part->getPy()*part->getPy()));  
      edm::LogInfo( "MUSiCSkimmer|GenInfo" ) << "GenMET after muon corr: Px = " << part->getPx()
                                             << ", Py = " << part->getPy()
                                             << ", Pt = " << part->getPt();
   }
   if (METMC_cuts(part)) numMETMC++;
   EvtView->setUserRecord<int>("NumMET", numMETMC);
   if( numMETMC ) edm::LogInfo( "MUSiCSkimmer|GenInfo" ) << "Event contains MET";
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
void MUSiCSkimmer::analyzeHCALNoise(const edm::Event& iEvent, pxl::EventView* EvtView) {
   //save HCAL noise infos
   edm::Handle< bool > hcal_noise;
   iEvent.getByLabel( hcal_noise_label, hcal_noise );
   EvtView->setUserRecord< bool >( "HCALNoisy", !*hcal_noise );

}


void MUSiCSkimmer::analyzeRecMET(const edm::Event& iEvent, pxl::EventView* EvtView, const collection_def &MET_info ) {
   edm::Handle<std::vector<pat::MET> > METHandle;
   iEvent.getByLabel(MET_info.RecoLabel, METHandle);
   const std::vector<pat::MET>& METs = *METHandle;
   std::vector<pat::MET>::const_iterator met = METs.begin();

   int numMETRec = 0;
   pxl::Particle* part = EvtView->create<pxl::Particle>();
   part->setName(MET_info.name);
   part->setP4(met->px(), met->py(), met->pz(), met->energy());
   part->setUserRecord<double>("sumEt", met->sumEt());
   part->setUserRecord<double>("mEtSig", met->mEtSig());

   if (MET_cuts(part)) numMETRec++;
   EvtView->setUserRecord<int>("NumMET", numMETRec);
   /*cout << "MET has " << met->nCorrections() << " corrections applied" << endl;
     cout << " Fully Corrected MET    : " << met->pt() << " (x: " << met->px() << ", y: " << met->py() << ") " << endl
     << " Uncorrected MET        : " << met->uncorrectedPt(pat::MET::uncorrALL) << "(x: " << met->corEx(pat::MET::uncorrALL) << ", y: " << met->corEy(pat::MET::uncorrALL) << ") " <<  endl
     << " Non-JES Corrected MET  : " << met->uncorrectedPt(pat::MET::uncorrJES) << "(x: " << met->corEx(pat::MET::uncorrJES) << ", y: " << met->corEy(pat::MET::uncorrJES) << ") " << endl
     << " Non-Muon Corrected MET : " << met->uncorrectedPt(pat::MET::uncorrMUON) << "(x: " << met->corEx(pat::MET::uncorrMUON) << ", y: " << met->corEy(pat::MET::uncorrMUON) << ") " << endl;
   */
}







void MUSiCSkimmer::initializeTrigger( const edm::Event &event,
                                      const edm::EventSetup &setup,
                                      trigger_group &trigger,
                                      const std::string &process
                                      ) {
   //read the new trigger config, test for error and whether something has changed
   bool changed;
   if( ! trigger.config.init( event.getRun(), setup, process, changed ) ){
      throw cms::Exception( "TRIGGER ERROR" ) << "Initialization of trigger config failed.";
   }

   //the trigger config has actually changed, so read in the new one
   if( changed ){
      edm::LogInfo( "MUSiCSkimmer|TRIGGERINFO" ) << "TRIGGER INFO: HLT table changed in run " << event.run()
                                                 << ", building new trigger map for process " << process;
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
            edm::LogWarning( "TRIGGERWARNING" ) << "In run " << event.run() << " trigger "<< *trg_name
                                                << " not found in HLT config, not added to trigger map (so not used).";
         }
      }
   }
}


void MUSiCSkimmer::analyzeFilter( const edm::Event &iEvent,
                                  const edm::EventSetup &iSetup,
                                  pxl::EventView *EvtView,
                                  trigger_group &filter
                                  ) {
   initializeTrigger( iEvent, iSetup, filter, filter.process );

   edm::Handle< edm::TriggerResults > filterResultsHandle;
   iEvent.getByLabel( filter.results, filterResultsHandle );

   for( vector< trigger_def >::iterator filt = filter.trigger_infos.begin(); filt != filter.trigger_infos.end(); ++filt ) {
      if( !filt->active ) continue;

      bool wasrun = filterResultsHandle->wasrun( filt->ID );
      bool error  = filterResultsHandle->error( filt->ID );

      if( wasrun && !error ){
        EvtView->setUserRecord< bool >( filter.name + "_" + filt->name, filterResultsHandle->accept( filt->ID ) );

        if( filterResultsHandle->accept( filt->ID ) )
           edm::LogInfo( "MUSiCSkimmer|FilterInfo" ) << "Event in process: '" << filter.process << "' passed filter: '" << filt->name << "'.";
      } else {
         //either error or was not run
         if( !wasrun ) edm::LogWarning( "FilterWarning" ) << "Filter: " << filt->name << " in process " << filter.process << " was not executed!";
         if( error )   edm::LogWarning( "FilterWarning" ) << "An error occured during execution of Filter: " << filt->name << " in process " << filter.process;
      }
   }
}


void MUSiCSkimmer::analyzeTrigger( const edm::Event &iEvent,
                                   const edm::EventSetup &iSetup,
                                   pxl::EventView* EvtView,
                                   trigger_group &trigger
                                   ){
   edm::Handle< trigger::TriggerEvent > triggerEventHandle;
   edm::Handle< edm::TriggerResults >   triggerResultsHandle;
   iEvent.getByLabel( trigger.event, triggerEventHandle );
   //try to find the right trigger if requested
   std::string process;
   if( trigger.process == "auto" ) {
      process = triggerEventHandle.provenance()->processName();
      edm::InputTag trigResultsTag( trigger.results.label(), trigger.results.instance(), process );
      iEvent.getByLabel( trigResultsTag, triggerResultsHandle );
   } else {
      process = trigger.results.process();
      iEvent.getByLabel( trigger.results, triggerResultsHandle );
   }
   //initialize the trigger config
   initializeTrigger( iEvent, iSetup, trigger, process );

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
            if( triggerResultsHandle->accept( trig->ID ) )
               edm::LogInfo( "MUSiCSkimmer|TriggerInfo" ) << "Trigger: " << trig->name << " in menu " << trigger.process << " fired" << endl;
         } else {
            //prescaled!
            //switch it off
            trig->active = false;
            edm::LogWarning( "TRIGGERWARNING" ) << "TRIGGER WARNING: Prescaled " << trig->name << " in menu " << trigger.process
                                                << " in run " << iEvent.run() << " - LS " << iEvent.luminosityBlock()
                                                << " - Event " << iEvent.id().event();
         }
      } else {
         //either error or was not run
         if( !wasrun ) edm::LogWarning( "TRIGGERWARNING" ) << "Trigger: " << trig->name << " in menu " << trigger.process << " was not executed!";
         if( error )   edm::LogWarning( "TRIGGERWARNING" ) << "An error occured during execution of Trigger: " << trig->name << " in menu " << trigger.process;
      }

      edm::LogInfo( "MUSiCSkimmer|TriggerInfo" ) << "triggerName: " << trig->name << "  triggerIndex: " << trig->ID << endl
                                                 << " Trigger path status:"
                                                 << " WasRun   = " << wasrun
                                                 << " Accept   = " << triggerResultsHandle->accept( trig->ID )
                                                 << " Error    = " << error
                                                 << " Prescale = " << prescale;

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
                  edm::LogInfo( "MUSiCSkimmer|TRIGGERINFO" ) << "Storing only 5 L3 objects for label/type "
                                                             << moduleLabel << "/" << trigger.config.moduleType( moduleLabel );
                  n = 5;
               }
               const trigger::TriggerObjectCollection &TOC = triggerEventHandle->getObjects();
               for( size_t i = 0; i != n; ++i ){
                  const trigger::TriggerObject &TO = TOC[KEYS[i]];
                  pxl::Particle *part = EvtView->create< pxl::Particle >();
                  part->setName( moduleLabel );
                  part->setP4( TO.px(), TO.py(), TO.pz(), TO.energy() );
                  part->setUserRecord< double >( "ID", TO.id() );
                  edm::LogVerbatim( "MUSiCSkimmer|TriggerInfo" )  << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
                                                                  << " id: " << TO.id() << " pt: " << TO.pt() << " eta: "
                                                                  << TO.eta() << " phi:" << TO.phi() << " m: " << TO.mass();
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

   //save the BS for further purpose
   the_beamspot = beamspot.position();

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

// ------------ reading Reconstructed Taus------------

void MUSiCSkimmer::analyzeRecTaus( const edm::Event &iEvent, pxl::EventView *RecView, const bool &MC, std::map< const reco::Candidate*, pxl::Particle* > &genmap ) {
   // get pat::Tau's from event
   edm::Handle< std::vector< pat::Tau > > tauHandle;
   iEvent.getByLabel( fTauRecoLabel, tauHandle );
   const std::vector< pat::Tau > &taus = *tauHandle;

   int numTauRec = 0;
   for( std::vector< pat::Tau >::const_iterator tau = taus.begin(); tau != taus.end(); ++tau ) {
      if( Tau_cuts( *tau ) ) {
         pxl::Particle *part = RecView->create< pxl::Particle >();
         part->setName( "Tau" );
         part->setCharge( tau->charge() );
         part->setP4( tau->px(), tau->py(), tau->pz(), tau->energy() );
         part->setUserRecord< float >( "PhiPhi", tau->phiphiMoment() );
         part->setUserRecord< float >( "EtaPhi", tau->etaphiMoment() );
         part->setUserRecord< float >( "EtaEta", tau->etaetaMoment() );
         //Pt of the Leading Charged Hadron of the Jet
         part->setUserRecord< float >( "LeadingHadronPt",    tau->leadPFChargedHadrCand()->pt() );
         part->setUserRecord< float >( "PfAllParticleIso",   tau->userIsolation( "pat::PfAllParticleIso" ) );
         part->setUserRecord< float >( "PfChargedHadronIso", tau->userIsolation( "pat::PfChargedHadronIso" ) );
         part->setUserRecord< float >( "PfNeutralHadronIso", tau->userIsolation( "pat::PfNeutralHadronIso" ) );
         part->setUserRecord< float >( "PfGammaIso",         tau->userIsolation( "pat::PfGammaIso" ) );
         //Saving all discriminators
         for ( std::vector< pat::Tau::IdPair >::const_iterator it = tau->tauIDs().begin(); it != tau->tauIDs().end(); ++it ) {
            part->setUserRecord < float >( it->first, it->second );
         }
         numTauRec++;
      }
   }
}
// ------------ reading Reconstructed Muons ------------

void MUSiCSkimmer::analyzeRecMuons( const edm::Event& iEvent, pxl::EventView* RecView, const bool& MC, std::map< const reco::Candidate*, pxl::Particle*> & genmap ) {
   // get pat::Muon's from event
   edm::Handle<std::vector<pat::Muon> > muonHandle;
   iEvent.getByLabel(fMuonRecoLabel, muonHandle);
   const std::vector<pat::Muon>& muons = *muonHandle;

   // TODO: The following method is depricated and should be replaced with
   // pat::muon::muonTrackFromMap( MuonTrackType ). The new method is only
   // available for samples reprocessed with CMSSW 4_4_X+. So it should be
   // implemented once Summer11 samples are not used anymore.
   //
   //get basic info for the refits
   Handle <reco::TrackToTrackMap> tevMapH1;
   Handle <reco::TrackToTrackMap> tevMapH2;
   Handle <reco::TrackToTrackMap> tevMapH3;
   iEvent.getByLabel("tevMuons", "default", tevMapH1);
   const reco::TrackToTrackMap tevMap1 = *(tevMapH1.product());
   iEvent.getByLabel("tevMuons", "firstHit", tevMapH2);
   const reco::TrackToTrackMap tevMap2 = *(tevMapH2.product());
   iEvent.getByLabel("tevMuons", "picky", tevMapH3);
   const reco::TrackToTrackMap tevMap3 = *(tevMapH3.product());

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

         // Particle-Flow muons out-of-the-box.
         part->setUserRecord< bool >( "isPFMuon", muon->isPFMuon() );
         if( muon->isPFMuon() ) {
            part->setUserRecord< double >( "PFpx", muon->pfP4().Px() );
            part->setUserRecord< double >( "PFpy", muon->pfP4().Py() );
            part->setUserRecord< double >( "PFpz", muon->pfP4().Pz() );
            part->setUserRecord< double >( "PFE",  muon->pfP4().E()  );
         }

         //store PAT matching info if MC
         if (MC){
            std::map< const reco::Candidate*, pxl::Particle* >::const_iterator it = genmap.find( muon->genLepton() );
            if( it != genmap.end() ){
               part->linkSoft( it->second, "pat-match" );
            }
         }

         //check if muon is Global/Tracker/StandAlone -Muon
         part->setUserRecord<bool>("isGlobalMuon", muon->isGlobalMuon());
         part->setUserRecord<bool>("isTrackerMuon", muon->isTrackerMuon());
         part->setUserRecord<bool>("isStandAloneMuon", muon->isStandAloneMuon());

         // Absolute transverse impact parameter.
         part->setUserRecord< double >( "dB", muon->dB() );

         //save info about quality of track-fit for combined muon (muon system + tracker)
         reco::TrackRef muontrack = muon->globalTrack();
         reco::TrackRef trackerTrack = muon->innerTrack();
         reco::TrackRef outerTrack = muon->outerTrack();

         part->setUserRecord<double>("NormChi2", muontrack->normalizedChi2());

         // Store info from HitPattern of the global tracker.

         // Number of lost ( = invalid) hits on track.
         //
         part->setUserRecord< int >( "LHits", muontrack->hitPattern().numberOfLostHits() );

         // Valid hit information
         //
         part->setUserRecord<int>("VHits",muontrack->hitPattern().numberOfValidHits());
         part->setUserRecord<int>("VHitsPixel",muontrack->hitPattern().numberOfValidPixelHits());
         part->setUserRecord<int>("VHitsTracker",muontrack->hitPattern().numberOfValidTrackerHits());
         part->setUserRecord<int>("VHitsMuonSys",muontrack->hitPattern().numberOfValidMuonHits());

         // Store the number of tracker layers with measurement.
         //
         part->setUserRecord< int >( "TrackerLayersWithMeas", muontrack->hitPattern().trackerLayersWithMeasurement() );
         part->setUserRecord< int >( "PixelLayersWithMeas",   muontrack->hitPattern().pixelLayersWithMeasurement() );

         //store the number of muon stations containing segments
         part->setUserRecord< int > ( "NMatchedStations", muon->numberOfMatchedStations() );

         // Store the pt and error from the global track.
         // ( qoverpError() is the same as error(0) for a track. )
         //
         part->setUserRecord< double >( "qoverp",      muontrack->qoverp() );
         part->setUserRecord< double >( "qoverpError", muontrack->qoverpError() );
         part->setUserRecord< double >( "ptError",     muontrack->ptError() );
         part->setUserRecord< double >( "pt",          muontrack->pt() );

         // TODO: These variables are still used in the analysis and should be
         // replaced with those above in the future.
         //
         //error info also used in muon-Met corrections, thus store variable to save info for later re-corrections
         part->setUserRecord<double>("dPtRelTrack", muontrack->error(0)/(muontrack->qoverp()));
         part->setUserRecord<double>("dPtRelTrack_off", muontrack->ptError()/muontrack->pt());

         // Store also the pt error from the tracker track.
         // ( qoverpError() is the same as error(0) for a track. )
         //
         part->setUserRecord< double >( "qoverpTracker",      trackerTrack->qoverp() );
         part->setUserRecord< double >( "qoverpErrorTracker", trackerTrack->qoverpError() );
         part->setUserRecord< double >( "ptErrorTracker",     trackerTrack->ptError() );
         part->setUserRecord< double >( "ptTracker",          trackerTrack->pt() );

         // Save distance to the primary vertex and the beam spot in z and xy plane, respectively
         // (i.e. the impact parameter)
         part->setUserRecord< double >( "Dsz", muontrack->dsz( the_vertex ) );
         part->setUserRecord< double >( "Dxy", muontrack->dxy( the_vertex ) );

         part->setUserRecord< double >( "DszBS", muontrack->dsz( the_beamspot ) );
         part->setUserRecord< double >( "DxyBS", muontrack->dxy( the_beamspot ) );

         part->setUserRecord< double >( "Dz",   trackerTrack->dz( the_vertex ) );
         part->setUserRecord< double >( "DzBS", trackerTrack->dz( the_beamspot ) );

         // TODO: Depricated, see above.
         // Store information for "cocktail" high energy refit.
         reco::TrackRef pmcTrack = muon::tevOptimized( *muon, tevMap1, tevMap2, tevMap3 );
         if( pmcTrack.isAvailable() ) {
            part->setUserRecord< bool >( "validCocktail", true );
            // Same as above but for cocktail
            //
            part->setUserRecord< double >( "pxCocktail", pmcTrack->px() );
            part->setUserRecord< double >( "pyCocktail", pmcTrack->py() );
            part->setUserRecord< double >( "pzCocktail", pmcTrack->pz() );

            part->setUserRecord< double >( "qoverpCocktail",      pmcTrack->qoverp() );
            part->setUserRecord< double >( "qoverpErrorCocktail", pmcTrack->qoverpError() );
            part->setUserRecord< double >( "ptErrorCocktail",     pmcTrack->ptError() );
            part->setUserRecord< double >( "ptCocktail",          pmcTrack->pt() );

            part->setUserRecord< double >( "NormChi2Cocktail", pmcTrack->normalizedChi2() );

            part->setUserRecord< int >( "LHitsCocktail",        pmcTrack->hitPattern().numberOfLostHits() );
            part->setUserRecord< int >( "VHitsCocktail",        pmcTrack->hitPattern().numberOfValidHits() );
            part->setUserRecord< int >( "VHitsPixelCocktail",   pmcTrack->hitPattern().numberOfValidPixelHits() );
            part->setUserRecord< int >( "VHitsTrackerCocktail", pmcTrack->hitPattern().numberOfValidTrackerHits() );
            part->setUserRecord< int >( "VHitsMuonSysCocktail", pmcTrack->hitPattern().numberOfValidMuonHits() );

            part->setUserRecord< double >( "DzCocktail",    pmcTrack->dz( the_vertex ) );
            part->setUserRecord< double >( "DszCocktail",   pmcTrack->dsz( the_vertex ) );
            part->setUserRecord< double >( "DxyCocktail",   pmcTrack->dxy( the_vertex ) );

            part->setUserRecord< double >( "DzBSCocktail",  pmcTrack->dz( the_vertex ) );
            part->setUserRecord< double >( "DszBSCocktail", pmcTrack->dsz( the_beamspot ) );
            part->setUserRecord< double >( "DxyBSCocktail", pmcTrack->dxy( the_beamspot ) );

            part->setUserRecord< int >( "TrackerLayersWithMeasCocktail", muontrack->hitPattern().trackerLayersWithMeasurement() );
            part->setUserRecord< int >( "PixelLayersWithMeasCocktail",   muontrack->hitPattern().pixelLayersWithMeasurement() );
         } else {
            part->setUserRecord< bool >( "validCocktail", false );
         }

         //official CaloIso and TrkIso
         //Def:  aMuon.setCaloIso(aMuon.isolationR03().emEt + aMuon.isolationR03().hadEt + aMuon.isolationR03().hoEt);
         part->setUserRecord<double>("CaloIso", muon->caloIso());
         part->setUserRecord<double>("TrkIso", muon->trackIso());
         part->setUserRecord<double>("ECALIso", muon->ecalIso());
         part->setUserRecord<double>("HCALIso", muon->hcalIso());
         //save offical isolation information: delta R = 0.3
         const reco::MuonIsolation& muonIsoR03 = muon->isolationR03();
         part->setUserRecord<double>("IsoR3SumPt", muonIsoR03.sumPt);
         part->setUserRecord<double>("IsoR3EmEt", muonIsoR03.emEt);
         part->setUserRecord<double>("IsoR3HadEt", muonIsoR03.hadEt);
         part->setUserRecord<double>("IsoR3HoEt", muonIsoR03.hoEt);
         part->setUserRecord<int>("IsoR3NTracks", muonIsoR03.nTracks);
         part->setUserRecord<int>("IsoR3NJets", muonIsoR03.nJets);
         //save offical isolation information: delta R = 0.5
         const reco::MuonIsolation& muonIsoR05 = muon->isolationR05();
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
         part->setUserRecord< bool >( "TMOneStationTight", muon::isGoodMuon( *muon, muon::TMOneStationTight ) );

         // Particle Flow based Isolation (available for samples processed with CMSSW_4_4_X+ otherwise 0).
         const reco::MuonPFIsolation muonPFIso03 = muon->pfIsolationR03();
         const reco::MuonPFIsolation muonPFIso04 = muon->pfIsolationR04();

         // Sum Pt of the charged Hadrons.
         part->setUserRecord< double >( "PFIsoR03ChargedHadrons", muonPFIso03.sumChargedHadronPt );
         part->setUserRecord< double >( "PFIsoR04ChargedHadrons", muonPFIso04.sumChargedHadronPt );
         // Sum Pt of all charged particles (including PF electrons and muons).
         part->setUserRecord< double >( "PFIsoR03ChargeParticles", muonPFIso03.sumChargedParticlePt );
         part->setUserRecord< double >( "PFIsoR04ChargeParticles", muonPFIso04.sumChargedParticlePt );
         // Sum Et of the neutral hadrons.
         part->setUserRecord< double >( "PFIsoR03NeutralHadrons", muonPFIso03.sumNeutralHadronEt );
         part->setUserRecord< double >( "PFIsoR04NeutralHadrons", muonPFIso04.sumNeutralHadronEt );
         // Sum Et of PF photons.
         part->setUserRecord< double >( "PFIsoR03Photons", muonPFIso03.sumPhotonEt );
         part->setUserRecord< double >( "PFIsoR04Photons", muonPFIso04.sumPhotonEt );
         // Sum Pt of the charged particles in the cone of interest but with particles not originating from the primary vertex(for PU corrections).
         part->setUserRecord< double >( "PFIsoR03PU", muonPFIso03.sumPUPt );
         part->setUserRecord< double >( "PFIsoR04PU", muonPFIso04.sumPUPt );
         // TODO: The following two are not available before CMSSW_5_X_Y.
         // Sum of the neutral hadron Et with a higher threshold for the candidates(1 GeV instead of 0.5).
         //part->setUserRecord< double >( "PFIso03NeutralHadronsHighThres", muonPFIso03.sumNeutralHadronEtHighThreshold );
         //part->setUserRecord< double >( "PFIso04NeutralHadronsHighThres", muonPFIso04.sumNeutralHadronEtHighThreshold );
         // Sum of the PF photons Et with higher threshold (1 GeV instead of 0.5).
         //part->setUserRecord< double >( "PFIso03PhotonsHighThres", muonPFIso03.sumPhotonEtHighThreshold );
         //part->setUserRecord< double >( "PFIso04PhotonsHighThres", muonPFIso04.sumPhotonEtHighThreshold );

         numMuonRec++;
      }
   }
   RecView->setUserRecord<int>("NumMuon", numMuonRec);
   edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "Rec Muons: " << numMuonRec;
}

// ------------ reading Reconstructed Electrons ------------

void MUSiCSkimmer::analyzeRecElectrons( const Event &iEvent,
                                        pxl::EventView *RecView,
                                        const bool &MC,
                                        EcalClusterLazyTools &lazyTools,
                                        map< const reco::Candidate*, pxl::Particle*> &genmap,
                                        const ESHandle< CaloGeometry > &geo
                                        ) {
   int numEleRec = 0;
   int numEleAll = 0;   // for matching

   Handle< vector< pat::Electron > > electronHandle;
   iEvent.getByLabel( fElectronRecoLabel, electronHandle );
   const vector< pat::Electron > &patElectrons = *electronHandle;

   Handle< EcalRecHitCollection > barrelRecHits;
   iEvent.getByLabel( freducedBarrelRecHitCollection, barrelRecHits );

   Handle< EcalRecHitCollection > endcapRecHits;
   iEvent.getByLabel( freducedEndcapRecHitCollection, endcapRecHits );

   Handle< reco::ConversionCollection > conversionsHandle;
   iEvent.getByLabel( m_conversionsTag, conversionsHandle );

   for( vector< pat::Electron>::const_iterator patEle = patElectrons.begin(); patEle != patElectrons.end(); ++patEle ) {
      if( Ele_cuts( patEle ) ) {
         edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "Electron Energy scale corrected: " << patEle->isEnergyScaleCorrected() << endl;

         Handle< EcalRecHitCollection > recHits;

         const bool isBarrel = patEle->isEB();
         const bool isEndcap = patEle->isEE();

         if( isBarrel ) recHits = barrelRecHits;
         if( isEndcap ) recHits = endcapRecHits;

         pxl::Particle *pxlEle = RecView->create< pxl::Particle >();
         pxlEle->setName( "Ele" );
         pxlEle->setCharge( patEle->charge() );
         pxlEle->setP4( patEle->px(), patEle->py(), patEle->pz(), patEle->ecalEnergy() );

         pxlEle->setUserRecord< bool >( "isBarrel", isBarrel );
         pxlEle->setUserRecord< bool >( "isEndcap", isEndcap );
         pxlEle->setUserRecord< bool >( "isGap", patEle->isEBGap() || patEle->isEEGap() || patEle->isEBEEGap() );

         //
         // Electron variables orientated to
         // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDInputVariables
         //

         pxlEle->setUserRecord< double >( "PErr",   patEle->trackMomentumError() );
         pxlEle->setUserRecord< double >( "SCeta",  patEle->caloPosition().eta() );
         pxlEle->setUserRecord< double >( "SCEErr", patEle->ecalEnergyError() );

         // Variables for isolation:
         //
         // Track iso deposit with electron footprint removed.
         pxlEle->setUserRecord< double >( "TrkIso03", patEle->dr03TkSumPt() );
         pxlEle->setUserRecord< double >( "TrkIso04", patEle->dr04TkSumPt() ); // (Identical to trackIso()!)
         // ECAL iso deposit with electron footprint removed.
         pxlEle->setUserRecord< double >( "ECALIso03", patEle->dr03EcalRecHitSumEt() );
         pxlEle->setUserRecord< double >( "ECALIso04", patEle->dr04EcalRecHitSumEt() ); // (Identical to ecalIso()!)
         // dr03HcalDepth1TowerSumEt()+dr03HcalDepth2TowerSumEt().
         pxlEle->setUserRecord< double >( "HCALIso03", patEle->dr03HcalTowerSumEt() );
         // dr04HcalDepth1TowerSumEt()+dr04HcalDepth2TowerSumEt().
         pxlEle->setUserRecord< double >( "HCALIso04", patEle->dr04HcalTowerSumEt() ); // (Identical to hcalIso()!)
         // HCAL depth 1 iso deposit with electron footprint removed.
         pxlEle->setUserRecord< double >( "HCALIso03d1", patEle->dr03HcalDepth1TowerSumEt() );
         pxlEle->setUserRecord< double >( "HCALIso04d1", patEle->dr04HcalDepth1TowerSumEt() );
         // HCAL depth 2 iso deposit with electron footprint removed.
         pxlEle->setUserRecord< double >( "HCALIso03d2", patEle->dr03HcalDepth2TowerSumEt() );
         pxlEle->setUserRecord< double >( "HCALIso04d2", patEle->dr04HcalDepth2TowerSumEt() );

         // Track-cluster matching variables.
         //
         // The supercluster eta - track eta position at calo extrapolated from innermost track state.
         pxlEle->setUserRecord< double >( "DEtaSCVtx", patEle->deltaEtaSuperClusterTrackAtVtx() );
         // The supercluster phi - track phi position at calo extrapolated from the innermost track state.
         pxlEle->setUserRecord< double >( "DPhiSCVtx", patEle->deltaPhiSuperClusterTrackAtVtx() );
         // The electron cluster eta - track eta position at calo extrapolated from the outermost state.
         pxlEle->setUserRecord< double >( "DEtaSCCalo", patEle->deltaEtaEleClusterTrackAtCalo() );
         // The seed cluster eta - track eta at calo from outermost state.
         pxlEle->setUserRecord< double >( "DEtaSeedTrk", patEle->deltaEtaSeedClusterTrackAtCalo() );
         // The seed cluster phi - track phi position at calo extrapolated from the outermost track state.
         pxlEle->setUserRecord< double >( "DPhiSeedTrk", patEle->deltaPhiSeedClusterTrackAtCalo() );
         // The seed cluster energy / track momentum at the PCA to the beam spot.
         pxlEle->setUserRecord< double >( "ESCSeedOverP", patEle->eSeedClusterOverP() );
         // The seed cluster energy / track momentum at calo extrapolated from the outermost track state.
         pxlEle->setUserRecord< double >( "ESCSeedPout", patEle->eSeedClusterOverPout() );
         // The supercluster energy / track momentum at the PCA to the beam spot.
         pxlEle->setUserRecord< double >( "EoP", patEle->eSuperClusterOverP() );
         // The electron cluster energy / track momentum at calo extrapolated from the outermost track state.
         pxlEle->setUserRecord< double >( "ESCOverPout", patEle->eEleClusterOverPout() );

         // Calorimeter information.
         //
         pxlEle->setUserRecord< double >( "sigmaIetaIeta", patEle->sigmaIetaIeta() );
         // Energy inside 1x5 in etaxphi around the seed Xtal.
         pxlEle->setUserRecord< double >( "e1x5", patEle->e1x5() );
         // Energy inside 2x5 in etaxphi around the seed Xtal (max bwt the 2 possible sums).
         pxlEle->setUserRecord< double >( "e2x5", patEle->e2x5Max() );
         // Energy inside 5x5 in etaxphi around the seed Xtal.
         pxlEle->setUserRecord< double >( "e5x5", patEle->e5x5() );
         // hcal over ecal seed cluster energy using first hcal depth (hcal is energy of towers within dR = 0.15).
         pxlEle->setUserRecord< double >( "HCALOverECALd1", patEle->hcalDepth1OverEcal() );
         // hadronicOverEm() = hcalDepth1OverEcal() + hcalDepth2OverEcal()
         const double HoEm = patEle->hadronicOverEm();
         pxlEle->setUserRecord< double >( "HoEm", HoEm );
         // Number of basic clusters inside the supercluster - 1.
         pxlEle->setUserRecord< double >( "NumBrems", patEle->numberOfBrems() );

         // Track information
         //
         // The brem fraction from gsf fit:
         // (track momentum in - track momentum out) / track momentum in
         pxlEle->setUserRecord< double >( "fbrem", patEle->fbrem() );
         pxlEle->setUserRecord< double >( "pin",   patEle->trackMomentumAtVtx().R() );
         pxlEle->setUserRecord< double >( "pout",  patEle->trackMomentumOut().R() );

         const reco::GsfTrackRef gsfTrack = patEle->gsfTrack();
         pxlEle->setUserRecord< double >( "TrackerP", gsfTrack->p() );

         pxlEle->setUserRecord< double >( "GSFNormChi2", gsfTrack->normalizedChi2() );
         // Save distance to the primary vertex and the beam spot, respectively (i.e. the impact parameter).
         pxlEle->setUserRecord< double >( "Dz",  gsfTrack->dz( the_vertex ) );
         pxlEle->setUserRecord< double >( "Dsz", gsfTrack->dsz( the_vertex ) );
         pxlEle->setUserRecord< double >( "Dxy", gsfTrack->dxy( the_vertex ) );

         pxlEle->setUserRecord< double >( "DzBS",  gsfTrack->dz( the_beamspot ) );
         pxlEle->setUserRecord< double >( "DszBS", gsfTrack->dsz( the_beamspot ) );
         pxlEle->setUserRecord< double >( "DxyBS", gsfTrack->dxy( the_beamspot ) );

         // Store the number of *expected* crossed layers before the first trajectory's hit.
         // If this number is 0, this is the number of missing hits in that trajectory. (This is for conversion rejection.)
         pxlEle->setUserRecord< int >( "NinnerLayerLostHits", gsfTrack->trackerExpectedHitsInner().numberOfHits() );
         pxlEle->setUserRecord< int >( "TrackerVHits", gsfTrack->numberOfValidHits() );
         pxlEle->setUserRecord< int >( "TrackerLHits", gsfTrack->numberOfLostHits() );

         // True if the electron track had an ecalDriven seed (regardless of the
         // result of the tracker driven seed finding algorithm).
         pxlEle->setUserRecord< bool >( "ecalDriven", patEle->ecalDrivenSeed() );
         // True if ecalDrivenSeed is true AND the electron passes the cut based
         // preselection.
         pxlEle->setUserRecord< bool >( "ecalDrivenEle", patEle->ecalDriven() );

         // Conversion rejection variables.
         //
         // Difference of cot(angle) with the conversion partner track.
         pxlEle->setUserRecord< double >( "convDcot", patEle->convDcot() );
         // Distance to the conversion partner.
         pxlEle->setUserRecord< double >( "convDist", patEle->convDist() );
         // Signed conversion radius.
         pxlEle->setUserRecord< double >( "convRadius", patEle->convRadius() );

         // Vertex coordinates.
         //
         pxlEle->setUserRecord< double >( "Vtx_X", patEle->vx() );
         pxlEle->setUserRecord< double >( "Vtx_Y", patEle->vy() );
         pxlEle->setUserRecord< double >( "Vtx_Z", patEle->vz() );

         // Electron classification: UNKNOWN=-1, GOLDEN=0, BIGBREM=1, OLDNARROW=2, SHOWERING=3, GAP=4.
         pxlEle->setUserRecord< int >( "Class", patEle->classification() );

         // Additional cluster variables for (spike) cleaning:
         //
         // Get the supercluster (ref) of the Electron
         // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
         // a SuperClusterCollection is a std::vector<SuperCluster>
         // although we get a vector of SuperClusters an electron is only made out of ONE SC
         // therefore only the first element of the vector should be available!
         const reco::SuperClusterRef SCRef = patEle->superCluster();

         const double SCenergy = SCRef->energy();
         pxlEle->setUserRecord< double >( "SCE", SCenergy );

         // Get highest energy entry (seed) and SC ID.
         // Use EcalClusterLazyTools to store ClusterShapeVariables.
         //
         const pair< DetId, float > max_hit = lazyTools.getMaximum( *SCRef );
         const DetId seedID = max_hit.first;
         const double eMax  = max_hit.second;
         const double e3x3  = lazyTools.e3x3( *SCRef ); // Energy in 3x3 around most energetic hit.

         pxlEle->setUserRecord< double >( "Emax", eMax );
         pxlEle->setUserRecord< double >( "E2nd", lazyTools.e2nd( *SCRef ) );
         pxlEle->setUserRecord< double >( "e3x3", e3x3 );
         pxlEle->setUserRecord< double >( "r19",  eMax / e3x3 );

         // These are the covariances, if you want the sigmas, you have to sqrt them.
         const vector< float > covariances = lazyTools.covariances( *SCRef, 4.7 );
         pxlEle->setUserRecord< double >( "covEtaEta", covariances[0] );
         pxlEle->setUserRecord< double >( "covEtaPhi", covariances[1] );
         pxlEle->setUserRecord< double >( "covPhiPhi", covariances[2] );

         // SwissCross
         //
         double swissCross         = -1.0;
         double swissCrossNoBorder = -1.0;

         if( isBarrel || isEndcap ) {
            swissCross         = EcalTools::swissCross( seedID, *recHits, 0, false );
            swissCrossNoBorder = EcalTools::swissCross( seedID, *recHits, 0, true );

            EcalRecHitCollection::const_iterator recHit_it = recHits->find( seedID );
            if( recHit_it != recHits->end() ) {
               pxlEle->setUserRecord< unsigned int >( "recoFlag", recHit_it->recoFlag() );
            }
         }

         pxlEle->setUserRecord< double >( "SwissCross",         swissCross );
         pxlEle->setUserRecord< double >( "SwissCrossNoBorder", swissCrossNoBorder );

         // Save eta/phi and DetId info from seed-cluster to prevent duplication of Electron/Photon-Candidates (in final selection).
         pxlEle->setUserRecord< unsigned int >( "seedId", seedID.rawId() );
         pxlEle->setUserRecord< double >( "seedphi", geo->getPosition( seedID ).phi() );
         pxlEle->setUserRecord< double >( "seedeta", geo->getPosition( seedID ).eta() );

         // Store ID information. ID is a string that typically starts with 'eid'.
         const vector< pair< string, float > > &electronIDs = patEle->electronIDs();
         for( vector< pair< string, float > >::const_iterator electronID = electronIDs.begin(); electronID != electronIDs.end(); ++electronID ){
            pxlEle->setUserRecord< bool >( electronID->first, electronID->second > 0.5 );
         }

         // Conversion veto for electron ID.
         // https://twiki.cern.ch/twiki/bin/view/CMS/ConversionTools
         //
         const bool hasMatchedConversion = ConversionTools::hasMatchedConversion( *patEle, conversionsHandle, the_beamspot );
         pxlEle->setUserRecord< bool >( "hasMatchedConversion", hasMatchedConversion );

         // 2012 definition of H/E and related HCAL isolation.
         // See also:
         // https://twiki.cern.ch/twiki/bin/view/CMS/HoverE2012
         //
         vector< CaloTowerDetId > hcalTowersBehindClusters = m_hcalHelper->hcalTowersBehindClusters( *SCRef );

         const double hcalDepth1 = m_hcalHelper->hcalESumDepth1BehindClusters( hcalTowersBehindClusters );
         const double hcalDepth2 = m_hcalHelper->hcalESumDepth2BehindClusters( hcalTowersBehindClusters );
         const double HoverE2012 = ( hcalDepth1 + hcalDepth2 ) / SCenergy;

         const double HCALIsoConeDR03_2012 = patEle->dr03HcalTowerSumEt() + ( HoEm - HoverE2012 ) * SCenergy / cosh( SCRef->eta() );

         pxlEle->setUserRecord< double >( "HoverE2012",           HoverE2012           );
         pxlEle->setUserRecord< double >( "HCALIsoConeDR03_2012", HCALIsoConeDR03_2012 );

         // Need a Ref to access the isolation values in particleFlowBasedIsolation(...).
         //
         pat::ElectronRef eleRef( electronHandle, numEleAll );

         if( eleRef.isNull() ){
            throw cms::Exception( "Reference Error" ) << "Could not create valid edm::Ref() to PAT electron "
                                                      << "(no. " << numEleAll << ")!";
         }

         particleFlowBasedIsolation( iEvent, m_inputTagIsoValElectronsPFId, eleRef, *pxlEle );

         // Store PAT matching info if MC. FIXME: Do we still use this?
         if( MC ) {
            map< const reco::Candidate*, pxl::Particle* >::const_iterator it = genmap.find( patEle->genLepton() );
            if( it != genmap.end() ){
               pxlEle->linkSoft( it->second, "pat-match" );
            }
         }

         numEleRec++;
      }
      numEleAll++;
   }
   RecView->setUserRecord< int >( "NumEle", numEleRec );

   edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "Rec Eles: " << numEleRec;
}

// ------------ reading Reconstructed Jets ------------

void MUSiCSkimmer::analyzeRecJets( const edm::Event &iEvent, pxl::EventView *RecView, bool &MC, std::map< const reco::Candidate*, pxl::Particle* > &genjetmap, const jet_def &jet_info ){
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
         if (jet_info.isPF) {
            part->setUserRecord<double>("chargedHadronEnergyFraction",jet->chargedHadronEnergyFraction());
            part->setUserRecord< double >( "neutralHadronEnergyFraction", jet->neutralHadronEnergyFraction() );
            part->setUserRecord<double>("chargedEmEnergyFraction", jet->chargedEmEnergyFraction());
            part->setUserRecord<double>("neutralEmEnergyFraction", jet->neutralEmEnergyFraction());
            part->setUserRecord<double>("chargedMultiplicity", jet->chargedMultiplicity());
            part->setUserRecord<double>("nconstituents", jet->numberOfDaughters());
         } else {
            part->setUserRecord<double>("EmEFrac", jet->emEnergyFraction());
            part->setUserRecord<double>("HadEFrac", jet->energyFractionHadronic());
            part->setUserRecord<int>("N90", jet->n90());
            part->setUserRecord<int>("N60", jet->n60());
            std::vector <CaloTowerPtr> caloRefs = jet->getCaloConstituents();
            part->setUserRecord<int>("NCaloRefs", caloRefs.size());
            part->setUserRecord<double>("MaxEEm", jet->maxEInEmTowers());
            part->setUserRecord<double>("MaxEHad", jet->maxEInHadTowers());
            part->setUserRecord<double>("TowersArea", jet->towersArea());
         }

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
            edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "BTag name: " << btag->first << ", value: " << btag->second;
         }
         //jet IDs
         for( jet_id_list::const_iterator ID = jet_info.IDs.begin(); ID != jet_info.IDs.end(); ++ID ){
            pat::strbitset ret = ID->second->getBitTemplate();
            part->setUserRecord< bool >( ID->first, (*(ID->second))( *jet, ret ) );
         }

         stringstream info;
         part->print( 0, info );
         edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "PXL Jet Info: " << info.str();

         //store PAT matching info if MC
         if (MC) {
            // to be compared with Generator Flavor:
            part->setUserRecord< int >( "algoFlavour", jet->partonFlavour() );
            //make a ref, then get the flavour
            RefToBase< reco::Jet > jetRef( RefToBaseProd< reco::Jet >( jetHandle ), jet_index );
            part->setUserRecord< int >( "physicsFlavour", (*physicsFlavour)[ jetRef ].getFlavour() );

            std::map< const reco::Candidate*, pxl::Particle* >::const_iterator it = genjetmap.find( jet->genJet() );
            if( it != genjetmap.end() ){
               part->linkSoft( it->second, "pat-match" );
            }
         }
         numJetRec++;
      }
   }
   RecView->setUserRecord< int >( "Num"+jet_info.name, numJetRec );
   edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "Found Rec Jets:  " << numJetRec << " of Type " << jet_info.name;
}

// ------------ reading Reconstructed Gammas ------------

void MUSiCSkimmer::analyzeRecGammas( const Event &iEvent,
                                     pxl::EventView *RecView,
                                     const bool &MC,
                                     EcalClusterLazyTools &lazyTools,
                                     map< const reco::Candidate*, pxl::Particle* > &genmap,
                                     const ESHandle< CaloGeometry > &geo
                                     ) {
   // Get Photon Collection.
   Handle< vector< pat::Photon > > photonHandle;
   iEvent.getByLabel( fGammaRecoLabel, photonHandle );
   const vector< pat::Photon > &patPhotons = *photonHandle;

   Handle< EcalRecHitCollection > barrelRecHits;
   iEvent.getByLabel( freducedBarrelRecHitCollection, barrelRecHits );

   Handle< EcalRecHitCollection > endcapRecHits;
   iEvent.getByLabel( freducedEndcapRecHitCollection, endcapRecHits );

   Handle< reco::ConversionCollection > conversionsHandle;
   iEvent.getByLabel( m_conversionsTag, conversionsHandle );

   Handle< reco::GsfElectronCollection > electronsHandle;
   iEvent.getByLabel( m_gsfElectronsTag, electronsHandle );

   int numGammaRec = 0;
   int numGammaAll = 0;
   for( vector< pat::Photon >::const_iterator patPhoton = patPhotons.begin(); patPhoton != patPhotons.end(); ++patPhoton ) {
      if( Gamma_cuts( patPhoton ) ) {
         Handle< EcalRecHitCollection > recHits;

         const bool isBarrel = patPhoton->isEB();
         const bool isEndcap = patPhoton->isEE();

         if( isBarrel ) recHits = barrelRecHits;
         if( isEndcap ) recHits = endcapRecHits;

         pxl::Particle *pxlPhoton = RecView->create< pxl::Particle >();
         pxlPhoton->setName( "Gamma" );
         pxlPhoton->setCharge( 0 );
         pxlPhoton->setP4( patPhoton->px(), patPhoton->py(), patPhoton->pz(), patPhoton->energy() );

         pxlPhoton->setUserRecord< bool >( "isBarrel", isBarrel );
         pxlPhoton->setUserRecord< bool >( "isEndcap", isEndcap );
         pxlPhoton->setUserRecord< bool >( "Gap", patPhoton->isEBGap() || patPhoton->isEEGap() || patPhoton->isEBEEGap() );

         // Store Photon info corrected for primary vertex (this changes direction but leaves energy of SC unchanged).
         pat::Photon localPho( *patPhoton );
         // Set event vertex
         localPho.setVertex( the_vertex );
         pxlPhoton->setUserRecord< double >( "PhysEta", localPho.eta() );
         pxlPhoton->setUserRecord< double >( "PhysPhi", localPho.phi() );
         pxlPhoton->setUserRecord< double >( "PhysPt",  localPho.pt()  );

         //
         // Photon variables orientated to
         // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDInputVariables
         //

         // Isolation variables.
         //
         // Sum of track pT in a hollow cone of outer radius, inner radius.
         pxlPhoton->setUserRecord< double >( "ID_TrkIso03", patPhoton->trkSumPtHollowConeDR03() );
         pxlPhoton->setUserRecord< double >( "ID_TrkIso",   patPhoton->trkSumPtHollowConeDR04() );
         // Sum of track pT in a cone of dR.
         pxlPhoton->setUserRecord< double >( "ID_TrkIsoSolid03", patPhoton->trkSumPtSolidConeDR03() );
         pxlPhoton->setUserRecord< double >( "ID_TrkIsoSolid",   patPhoton->trkSumPtSolidConeDR04() ); // (Identical to TrkIso()!)
         // EcalRecHit isolation.
         pxlPhoton->setUserRecord< double >( "ID_ECALIso03", patPhoton->ecalRecHitSumEtConeDR03() );
         pxlPhoton->setUserRecord< double >( "ID_ECALIso",   patPhoton->ecalRecHitSumEtConeDR04() ); // (Identical to ecalIso()!)
         // HcalDepth1Tower isolation.
         pxlPhoton->setUserRecord< double >( "ID_HCALIso03", patPhoton->hcalTowerSumEtConeDR03() );
         pxlPhoton->setUserRecord< double >( "ID_HCALIso",   patPhoton->hcalTowerSumEtConeDR04() ); // (Identical to hcalIso()!)
         // Number of tracks in a cone of dR.
         pxlPhoton->setUserRecord< int >( "TrackNum03", patPhoton->nTrkSolidConeDR03() );
         pxlPhoton->setUserRecord< int >( "TrackNum",   patPhoton->nTrkSolidConeDR04() );

         // Calorimeter information.
         //
         pxlPhoton->setUserRecord< double >( "iEta_iEta", patPhoton->sigmaIetaIeta() );
         pxlPhoton->setUserRecord< double >( "r9",        patPhoton->r9() );
         const double HoEm = patPhoton->hadronicOverEm();
         pxlPhoton->setUserRecord< double >( "HoEm", HoEm );

         // Whether or not the SuperCluster has a matched pixel seed (electron veto).
         pxlPhoton->setUserRecord< bool >( "HasSeed", patPhoton->hasPixelSeed() );

         // Store information about converted state.
         pxlPhoton->setUserRecord< bool >( "Converted", patPhoton->hasConversionTracks() );

         // Additional cluster variables for (spike) cleaning:
         //
         // Get the supercluster (ref) of the Electron
         // a SuperClusterRef is a edm::Ref<SuperClusterCollection>
         // a SuperClusterCollection is a std::vector<SuperCluster>
         // although we get a vector of SuperClusters an electron is only made out of ONE SC
         // therefore only the first element of the vector should be available!
         const reco::SuperClusterRef SCRef = patPhoton->superCluster();

         pxlPhoton->setUserRecord< double >( "etaWidth", SCRef->etaWidth() );
         pxlPhoton->setUserRecord< double >( "phiWidth", SCRef->phiWidth() );
         // Set hadronic over electromagnetic energy fraction.

         // Raw uncorrected energy (sum of energies of component BasicClusters).
         pxlPhoton->setUserRecord< double >( "rawEnergy", SCRef->rawEnergy() );
         // Energy deposited in preshower.
         pxlPhoton->setUserRecord< double >( "preshowerEnergy", SCRef->preshowerEnergy() );

         // Get highest energy entry (seed) and SC ID.
         // Use EcalClusterLazyTools to store ClusterShapeVariables.
         //
         const pair< DetId, float > max_hit = lazyTools.getMaximum( *SCRef );
         const DetId seedID = max_hit.first;
         const double eMax = max_hit.second;
         const double e3x3 = patPhoton->e3x3();

         pxlPhoton->setUserRecord< double >( "Emax", eMax );
         pxlPhoton->setUserRecord< double >( "E2nd", lazyTools.e2nd( *SCRef ) );
         pxlPhoton->setUserRecord< double >( "e3x3", e3x3 );
         pxlPhoton->setUserRecord< double >( "e5x5", patPhoton->e5x5() );
         pxlPhoton->setUserRecord< double >( "r19",  eMax / e3x3 );

         // SwissCross
         //
         double swissCross         = -1.0;
         double swissCrossNoBorder = -1.0;

         if( isBarrel || isEndcap ) {
            swissCross         = EcalTools::swissCross( seedID, *recHits, 0, false );
            swissCrossNoBorder = EcalTools::swissCross( seedID, *recHits, 0, true );

            EcalRecHitCollection::const_iterator recHit_it = recHits->find( seedID );
            if( recHit_it != recHits->end() ) {
               pxlPhoton->setUserRecord< unsigned int >( "recoFlag", recHit_it->recoFlag() );
            }
         }

         pxlPhoton->setUserRecord< double >( "SwissCross",         swissCross );
         pxlPhoton->setUserRecord< double >( "SwissCrossNoBorder", swissCrossNoBorder );

         // These are the covariances, if you want the sigmas, you have to sqrt them.
         const vector< float > covariances = lazyTools.covariances( *SCRef, 4.7 );
         pxlPhoton->setUserRecord< double >( "covEtaEta", covariances[0] );
         pxlPhoton->setUserRecord< double >( "covEtaPhi", covariances[1] );
         pxlPhoton->setUserRecord< double >( "covPhiPhi", covariances[2] );

         // Save eta/phi and DetId info from seed-cluster to prevent duplication of Electron/patPhoton-Candidates (in final selection)
         // and to reject converted photons.
         //
         pxlPhoton->setUserRecord< unsigned int >( "seedId", seedID.rawId() );
         pxlPhoton->setUserRecord< double >( "seedphi", geo->getPosition( seedID ).phi() );
         pxlPhoton->setUserRecord< double >( "seedeta", geo->getPosition( seedID ).eta() );

         // Store ID information.
         const vector< pair< string, bool > > &photonIDs = patPhoton->photonIDs();
         for( vector< pair< string, bool > >::const_iterator photonID = photonIDs.begin(); photonID != photonIDs.end(); ++photonID ){
            pxlPhoton->setUserRecord< bool >( photonID->first, photonID->second );
         }

         // Conversion safe electron veto for photon ID.
         // ("Conversion-safe" since it explicitly checks for the presence of a
         // reconstructed conversion.)
         // See also:
         // https://twiki.cern.ch/twiki/bin/view/CMS/ConversionTools
         //
         const bool hasMatchedPromptElectron = ConversionTools::hasMatchedPromptElectron( SCRef,
                                                                                          electronsHandle,
                                                                                          conversionsHandle,
                                                                                          the_beamspot );
         pxlPhoton->setUserRecord< bool >( "hasMatchedPromptElectron", hasMatchedPromptElectron );

         // 2012 definition of H/E and related HCAL isolation.
         // See also:
         // https://twiki.cern.ch/twiki/bin/view/CMS/HoverE2012
         //
         const vector< CaloTowerDetId > hcalTowersBehindClusters = m_hcalHelper->hcalTowersBehindClusters( *SCRef );

         const double hcalDepth1 = m_hcalHelper->hcalESumDepth1BehindClusters( hcalTowersBehindClusters );
         const double hcalDepth2 = m_hcalHelper->hcalESumDepth2BehindClusters( hcalTowersBehindClusters );
         const double HoverE2012 = ( hcalDepth1 + hcalDepth2 ) / SCRef->energy();

         const double HCALIsoConeDR03_2012 = patPhoton->hcalTowerSumEtConeDR03() + ( HoEm - HoverE2012 ) * SCRef->energy() / cosh( SCRef->eta() );

         pxlPhoton->setUserRecord< double >( "HoverE2012"           , HoverE2012           );
         pxlPhoton->setUserRecord< double >( "HCALIsoConeDR03_2012" , HCALIsoConeDR03_2012 );

         // Need a Ref to access the isolation values in particleFlowBasedIsolation(...).
         //
         pat::PhotonRef phoRef( photonHandle, numGammaAll );

         if( phoRef.isNull() ){
            throw cms::Exception( "Reference Error" ) << "Could not create valid edm::Ref() to PAT photon "
                                                      << "(no. " << numGammaAll << ")!";
         }

         particleFlowBasedIsolation( iEvent, m_inputTagIsoValPhotonsPFId, phoRef, *pxlPhoton );

         // Store PAT matching info if MC. FIXME: Do we still use this?
         if( MC ) {
            std::map< const reco::Candidate*, pxl::Particle* >::const_iterator it = genmap.find( patPhoton->genPhoton() );
            if( it != genmap.end() ){
               pxlPhoton->linkSoft( it->second, "pat-match" );
            }
         }
         numGammaRec++;
      }
      numGammaAll++;
   }
   RecView->setUserRecord< int >( "NumGamma", numGammaRec );

   edm::LogInfo( "MUSiCSkimmer|RecInfo" ) << "Rec Gammas: " << numGammaRec;
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
      pdfSet.append("/");
      pdfSet.append(fLHgridName);
      edm::LogInfo( "MUSiCSkimmer|PDFINFO" ) << "PDF set - " << pdfSet.data();
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
      for( int subpdf = 1; subpdf <= fNumLHgridErrorSets; subpdf++ ){
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
// ------------ method to define MC-TAU-cuts

bool MUSiCSkimmer::TauMC_cuts( const reco::GenParticle *MCtau ) const {
   if( MCtau->pt() < min_tau_pt ) return false;
   if( fabs( MCtau->eta() ) > max_eta ) return false;
   return true;
}

// ------------ method to define MC-MUON-cuts

bool MUSiCSkimmer::MuonMC_cuts( const reco::GenParticle *MCmuon ) const {
   if( MCmuon->pt() < min_muon_pt ) return false;
   if( fabs( MCmuon->eta() ) > max_eta ) return false;
   return true;
}
 


// ------------ method to define MC-Electron-cuts

bool MUSiCSkimmer::EleMC_cuts( const reco::GenParticle *MCele ) const {
   if( MCele->pt() < min_ele_pt ) return false;
   if( fabs( MCele->eta() ) > max_eta ) return false;
   return true;
}

// ------------ method to define MC-Gamma-cuts

bool MUSiCSkimmer::GammaMC_cuts( const reco::GenParticle *MCgamma ) const {
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


// ------------ method to define TAU-cuts

bool MUSiCSkimmer::Tau_cuts( const pat::Tau &tau ) const {
   // basic preselection cuts
   if( tau.pt() < min_tau_pt )  return false;
   if( fabs( tau.eta() ) > max_eta ) return false;
   return true;
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
      const reco::GenParticle* p = (const reco::GenParticle*) &(*pa);

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


// Accessing ParticleFlow based isolation (both methods):
// (See also: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation)
//
template< typename T >
void MUSiCSkimmer::particleFlowBasedIsolation( const Event &iEvent,
                                               const vector< InputTag > &inputTagIsoValPFId,
                                               const Ref< T > &ref,
                                               pxl::Particle &part,
                                               const bool &useIsolator ) const {

   const size_t numIsoVals = inputTagIsoValPFId.size();

   // typedef in MUSiCSkimmer.h
   IsoDepositVals isoValPFId( numIsoVals );

   for( size_t i = 0; i < numIsoVals; ++i ) {
      iEvent.getByLabel( inputTagIsoValPFId.at( i ), isoValPFId.at( i ) );
   }

   const double pfIsoCharged = ( *isoValPFId.at( 0 ) )[ ref ];
   const double pfIsoPhoton  = ( *isoValPFId.at( 1 ) )[ ref ];
   const double pfIsoNeutral = ( *isoValPFId.at( 2 ) )[ ref ];

   part.setUserRecord< double >( "PFIso03ChargedHadron", pfIsoCharged );
   part.setUserRecord< double >( "PFIso03NeutralHadron", pfIsoNeutral );
   part.setUserRecord< double >( "PFIso03Photon",        pfIsoPhoton  );

   // PU corrected isolation for electrons, according to:
   // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation#Example_for_photons
   // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/test/ElectronIsoAnalyzer.cc
   //
   if( ref->isElectron() ) {
      ElectronEffectiveArea::ElectronEffectiveAreaTarget eleEffAreaTarget;
      if     ( m_eleEffAreaTargetLabel == "NoCorr"     ) eleEffAreaTarget = ElectronEffectiveArea::kEleEANoCorr;
      else if( m_eleEffAreaTargetLabel == "Data2011"   ) eleEffAreaTarget = ElectronEffectiveArea::kEleEAData2011;
      else if( m_eleEffAreaTargetLabel == "Data2012"   ) eleEffAreaTarget = ElectronEffectiveArea::kEleEAData2012;
      else if( m_eleEffAreaTargetLabel == "Summer11MC" ) eleEffAreaTarget = ElectronEffectiveArea::kEleEASummer11MC;
      else if( m_eleEffAreaTargetLabel == "Fall11MC"   ) eleEffAreaTarget = ElectronEffectiveArea::kEleEAFall11MC;
      else throw cms::Exception( "Configuration" ) << "Unknown effective area " << m_eleEffAreaTargetLabel << endl;

      const ElectronEffectiveArea::ElectronEffectiveAreaType eleEffAreaGammaPlusNeutralHad = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03;

      const double absEta           = fabs( ref->superCluster()->eta() );
      const double effArea          = ElectronEffectiveArea::GetElectronEffectiveArea( eleEffAreaGammaPlusNeutralHad, absEta, eleEffAreaTarget );
      const double PFIsoPUCorrected = pfIsoCharged + max( 0.0, ( pfIsoPhoton + pfIsoNeutral ) - effArea * m_rhoFastJet );

      part.setUserRecord< double >( "EffectiveArea",      effArea          );
      part.setUserRecord< double >( "PFIso03PUCorrected", PFIsoPUCorrected );
   }

   if( useIsolator ) {
      // We need the PFCandidates ...
      //
      Handle< reco::PFCandidateCollection > pfCandidates;
      iEvent.getByLabel( m_particleFlowTag, pfCandidates );
      const PFCandidateCollection thePFCollection = *pfCandidates;

      // ... and the vertices.
      //
      Handle< reco::VertexCollection > vertices;
      iEvent.getByLabel( fVertexRecoLabel, vertices );

      // Primary Vertex
      const reco::VertexRef vtxRef( vertices, 0 );

      // Initialise isolator.
      //
      PFIsolationEstimator isolator;
      isolator.initializeElectronIsolation( kTRUE );
      isolator.setConeSize( 0.3 );

      // Analysing electrons:
      //
      if( ref->isElectron() ) isolator.fGetIsolation( &*ref, &thePFCollection, vtxRef, vertices );

      // Analysing photons:
      //
      if( ref->isPhoton() ) isolator.fGetIsolation( &*ref, &thePFCollection, vtxRef, vertices );

      part.setUserRecord< double >( "PFIso03ChargedHadronFromIsolator", isolator.getIsolationCharged() );
      part.setUserRecord< double >( "PFIso03NeutralHadronFromIsolator", isolator.getIsolationNeutral() );
      part.setUserRecord< double >( "PFIso03PhotonFromIsolator",        isolator.getIsolationPhoton()  );
   }
}

#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(MUSiCSkimmer);
