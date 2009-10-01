#ifndef MUSiCSkimmer_H
#define MUSiCSkimmer_H

// LHAPDF stuff
extern "C" {
   void initpdfset_ (char *, int len);
   void initpdfsetm_(int &, char *);
   void initpdf_(int &);
   void evolvepdf_(double &, double &, double *);
   void numberpdf_(int &);
}

struct PDFInf {
   int   f1; 
   int   f2; 
   double  x1;
   double  x2;
   double  Q; 
};

// CMSSW includes
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
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
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
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
#include "MUSiCProject/ePaxPxl/interface/ePaxPxl.h"
#include "MUSiCProject/Skimming/interface/ParticleMatcher.hh"


class MUSiCSkimmer : public edm::EDAnalyzer {
public:

   // why explicit?
   explicit MUSiCSkimmer(const edm::ParameterSet&);
   ~MUSiCSkimmer();
  

 private:

   virtual void beginJob(const edm::EventSetup&);
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();
   virtual void analyzeGenInfo(const edm::Event&, pxl::EventView*, std::map<const Candidate*, pxl::Particle*>&);
   virtual void analyzeGenRelatedInfo(const edm::Event&, pxl::EventView*);
   virtual void analyzeGenJets(const edm::Event&, pxl::EventView*, std::map<const Candidate*, pxl::Particle*>&);
   virtual void analyzeGenMET(const edm::Event&, pxl::EventView*);

   virtual void analyzeSIM(const edm::Event&, pxl::EventView*);
   
   virtual void analyzeTrigger( const edm::Event &iEvent,
                                pxl::EventView *EvtView,
                                const std::string &processName,
                                const edm::InputTag &HLTResults,
                                const edm::InputTag &L3Objets,
                                const std::map<int, std::string> &HLT_bit_to_name
                                );
   virtual void analyzeRecVertices(const edm::Event&, pxl::EventView*);
   virtual void analyzeRecMuons(const edm::Event&, pxl::EventView*, const bool&, std::map<const Candidate*, pxl::Particle*>&);
   virtual void analyzeRecElectrons(const edm::Event&, pxl::EventView*, bool&, EcalClusterLazyTools&, std::map<const Candidate*, pxl::Particle*>&);
   virtual void analyzeRecJets(const edm::Event&, pxl::EventView*, bool&, std::map<const Candidate*, pxl::Particle*>&);
   virtual void analyzeRecMET(const edm::Event&, pxl::EventView*);
   virtual void analyzeRecGammas(const edm::Event&, pxl::EventView*, bool&, EcalClusterLazyTools&, std::map<const Candidate*, pxl::Particle*>&);

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
   bool fUseSIM;
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
   std::string   triggerProcess;
   edm::InputTag ftriggerResultsTag; 
   edm::InputTag fTriggerEvent;
   std::map<int, std::string> fHLTMap;

   // Second trigger
   std::string   triggerProcess2;
   edm::InputTag ftriggerResultsTag2;
   edm::InputTag fTriggerEvent2;
   std::map<int, std::string> fHLTMap2;

   bool fStoreL3Objects;
 
   ParticleMatcher* Matcher;
   // to be used for ePax output 
   pxl::OutputFile fePaxFile;
   std::vector<PDFInf> fpdf_vec;
   double xfx(const double &x, const double &Q, int fl) {  
      double f[13], mx = x, mQ = Q;
      evolvepdf_(mx, mQ, f);
      return f[fl+6];
   };

};
#endif
