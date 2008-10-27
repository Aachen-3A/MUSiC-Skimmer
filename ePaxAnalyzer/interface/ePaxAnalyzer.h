#ifndef ePaxAnalyzer_H
#define ePaxAnalyzer_H
//
// class decleration
//

// CMSSW includes
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/SiStripElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/SiStripElectron.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

 
// ROOT stuff

#include "TFile.h"
#include "TTree.h"
#include "TMatrixT.h"

//PAT related stuff
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

//for ClusterShape variables
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

// ePax stuff
// Has to be included as the last header otherwise there will be a warning concerning the 
// zlib. According to Steffen there are two different zlib and ROOT can only deal with one of them
// but ePax can deal with both of them
#include "ePaxPxl/ePax/interface/ePax.h"
#include "ePaxDemo/ePaxAnalyzer/interface/ParticleMatcher.hh"


class ePaxAnalyzer : public edm::EDAnalyzer {
public:

   // why explicit?
   explicit ePaxAnalyzer(const edm::ParameterSet&);
   ~ePaxAnalyzer();
  

private:

   virtual void beginJob(const edm::EventSetup&);
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();
   virtual void analyzeGenInfo(const edm::Event&, pxl::EventView*, std::map<const Particle*, pxl::Particle*>&);
   virtual void analyzeGenRelatedInfo(const edm::Event&, pxl::EventView*);
   virtual void analyzeGenJets(const edm::Event&, pxl::EventView*, std::map<const Particle*, pxl::Particle*>&);
   virtual void analyzeGenMET(const edm::Event&, pxl::EventView*);
   
   virtual void analyzeTrigger(const edm::Event&, pxl::EventView*); //not complete!
   virtual void analyzeRecVertices(const edm::Event&, pxl::EventView*);
   virtual void analyzeRecMuons(const edm::Event&, pxl::EventView*, const bool&, std::map<const Particle*, pxl::Particle*>&);
   virtual void analyzeRecElectrons(const edm::Event&, pxl::EventView*, bool&, EcalClusterLazyTools&, std::map<const Particle*, pxl::Particle*>&);
   virtual void analyzeRecJets(const edm::Event&, pxl::EventView*, bool&, std::map<const Particle*, pxl::Particle*>&);
   virtual void analyzeRecMET(const edm::Event&, pxl::EventView*);
   virtual void analyzeRecGammas(const edm::Event&, pxl::EventView*, bool&, EcalClusterLazyTools&, std::map<const Particle*, pxl::Particle*>&);

   bool MuonMC_cuts(const GenParticle* MCmuon) const;
   bool EleMC_cuts(const GenParticle* MCele) const;
   bool GammaMC_cuts(const GenParticle* MCgamma) const;
   bool JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const;
   bool METMC_cuts(const pxl::Particle* MCmet) const;
   bool Vertex_cuts(reco::VertexCollection::const_iterator vertex) const; 
   bool Muon_cuts(const pat::Muon& muon) const;
   bool Ele_cuts(std::vector<pat::Electron>::const_iterator ele) const;
   bool Gamma_cuts(std::vector<pat::Photon>::const_iterator photon) const;
   bool Jet_cuts(std::vector<pat::Jet>::const_iterator jet) const;
   bool MET_cuts(const pxl::Particle* met) const;
   std::string getEventClass(pxl::EventView* EvtView);

   double IsoGenSum (const edm::Event& iEvent, double ParticleGenPt, double ParticleGenEta, double ParticleGenPhi, double iso_DR, double iso_Seed);

   // ----------member data ---------------------------

   int fNumEvt;// used to count the number of events
   int fDebug; 
   std::string fFileName; 
   std::string fProcess;
   bool fGenOnly;
   // The labels used in cfg-file 
   // Generator 
   std::string fgenParticleCandidatesLabel;
   std::string fMETMCLabel;
   std::string fVertexRecoLabel;
   // Muon
   std::string fMuonRecoLabel;
   // Electron
   std::string fElectronRecoLabel;
   // Photon
   std::string fGammaRecoLabel;
   // Jets labels
   std::vector<std::string> fJetRecoLabels;
   std::vector<std::string> fJetMCLabels;
   // MET labels
   std::string fMETRecoLabel;
   // Cluster
   edm::InputTag fBarrelClusterShapeAssocProducer;
   edm::InputTag fEndcapClusterShapeAssocProducer;
   std::string fHBHELabel;
   std::string fHBHEInstanceName;
   edm::InputTag fbarrelClusterCollection; 
   edm::InputTag fendcapClusterCollection;
   edm::InputTag freducedBarrelRecHitCollection;
   edm::InputTag freducedEndcapRecHitCollection;   
   // Trigger
   edm::InputTag ftriggerResultsTag; 
   edm::InputTag ftriggerEventTag; 

   ParticleMatcher* Matcher;
   // to be used for ePax output 
   pxl::OutputFile fePaxFile;

};
#endif
