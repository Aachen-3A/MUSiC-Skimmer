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

// EDM related stuff
// (Don't need headers here, forward declarations are enough!)
namespace edm {
   class Event;
   class EventSetup;
   class ParameterSet;
   template< class T > class ESHandle;
}

// GEN stuff.
namespace gen {
   class PdfInfo;
}

// RECO forward declarations.
namespace reco {
   class Candidate;
   class GenParticle;
   class Vertex;
}

// PAT related stuff.
namespace pat {
   class Electron;
   class Muon;
   class Photon;
   class Jet;
   class Tau;
}

// No namespace for these classes.
class CaloGeometry;
class ElectronHcalHelper;
class EcalClusterLazyTools;
class ParticleMatcher;

// STL
#include <map>
#include <string>
#include <vector>

// CMSSW includes
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

// GenJetCollection and related typedefs and forward declarations.
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// Beam Spot.
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

// HLT stuff.
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Jet selectors (need to be included here, otherwise there are (namespace) problems).
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

// Vertex typedefs / forward declarations.
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Private collection defintions.
#include "MUSiCProject/Skimming/interface/collection_def.h"

// PXL stuff
// Has to be included as the last header otherwise there will be a warning concerning the
// zlib. According to Steffen there are two different zlib and ROOT can only deal with one of them
// but PXL can deal with both of them
#include "MUSiCProject/Pxl/interface/Pxl.h"

class MUSiCSkimmer : public edm::EDAnalyzer {
public:

   // Why explicit?
   // The explicit keyword prevents the constructor from being invoked
   // implicitly as a conversion (what can be done for constructors with one
   // argument).
   explicit MUSiCSkimmer(const edm::ParameterSet&);
   ~MUSiCSkimmer();
  

 private:
   //information about one single trigger
   struct trigger_def {
      std::string name;
      unsigned int ID;
      bool active;
   };
   //information about one trigger group
   struct trigger_group {
      std::string   name;
      std::string   process;
      edm::InputTag L1_result;
      edm::InputTag results;
      edm::InputTag event;
      HLTConfigProvider config;
      std::vector< std::string > triggers_names;
      std::vector< trigger_def > trigger_infos;
   };


   virtual void analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup );
   virtual void endJob();
   virtual void analyzeGenInfo( const edm::Event &iEvent, pxl::EventView *EvtView, std::map< const reco::Candidate*, pxl::Particle* > &genmap );
   virtual void analyzeGenRelatedInfo(const edm::Event&, pxl::EventView*);
   virtual void analyzeGenJets( const edm::Event &iEvent, pxl::EventView *GenEvtView, std::map< const reco::Candidate*, pxl::Particle* > &genjetmap, const jet_def &jet_info );
   virtual void analyzeGenMET(const edm::Event&, pxl::EventView*, const collection_def &MET_info );

   virtual void analyzeSIM(const edm::Event&, pxl::EventView*);
   
   virtual void initializeTrigger( const edm::Event &event,
                                   const edm::EventSetup &setup,
                                   trigger_group &trigger,
                                   const std::string &process
                                   );

   virtual void analyzeFilter( const edm::Event &iEvent,
                               const edm::EventSetup &iSetup,
                               pxl::EventView *EvtView,
                               trigger_group &filter
                               );

   virtual void analyzeTrigger( const edm::Event &iEvent,
                                const edm::EventSetup &iSetup,
                                pxl::EventView *EvtView,
                                trigger_group &trigger
                                );
   virtual void analyzeRecVertices(const edm::Event&, pxl::EventView*);
   virtual void analyzeRecTaus( const edm::Event &iEvent, pxl::EventView *RecView, const bool &MC, std::map< const reco::Candidate*, pxl::Particle*> &genmap );
   virtual void analyzeRecMuons( const edm::Event &iEvent, pxl::EventView *RecView, const bool &MC, std::map< const reco::Candidate*, pxl::Particle* > &genmap );
   virtual void analyzeRecElectrons( const edm::Event &iEvent,
                                     pxl::EventView *RecView,
                                     const bool &MC,
                                     EcalClusterLazyTools &lazyTools,
                                     std::map< const reco::Candidate*, pxl::Particle* > &genmap,
                                     const edm::ESHandle< CaloGeometry > &geo
                                     );
   virtual void analyzeRecJets( const edm::Event &iEvent, pxl::EventView *RecView, bool &MC, std::map< const reco::Candidate*, pxl::Particle* > &genjetmap, const jet_def &jet_info );
   virtual void analyzeRecMET(const edm::Event&, pxl::EventView*, const collection_def &MET_info);
   virtual void analyzeRecGammas( const edm::Event &iEvent,
                                  pxl::EventView *RecView,
                                  const bool &MC,
                                  EcalClusterLazyTools &lazyTools,
                                  std::map< const reco::Candidate*, pxl::Particle* > &genmap,
                                  const edm::ESHandle< CaloGeometry > &geo
                                  );
   virtual void analyzeHCALNoise(const edm::Event&, pxl::EventView*);

   bool TauMC_cuts( const reco::GenParticle *MCtau ) const;
   bool MuonMC_cuts( const reco::GenParticle* MCmuon ) const;
   bool EleMC_cuts( const reco::GenParticle* MCele ) const;
   bool GammaMC_cuts( const reco::GenParticle* MCgamma ) const;
   bool JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const;
   bool METMC_cuts(const pxl::Particle* MCmet) const;
   bool Vertex_cuts(reco::VertexCollection::const_iterator vertex) const;
   bool PV_vertex_cuts( const reco::Vertex &vertex) const;
   bool Tau_cuts (const pat::Tau &tau) const;
   bool Muon_cuts(const pat::Muon& muon) const;
   bool Ele_cuts(std::vector<pat::Electron>::const_iterator ele) const;
   bool Gamma_cuts(std::vector<pat::Photon>::const_iterator photon) const;
   bool Jet_cuts(std::vector<pat::Jet>::const_iterator jet) const;
   bool MET_cuts(const pxl::Particle* met) const;
   std::string getEventClass(pxl::EventView* EvtView);

   double IsoGenSum (const edm::Event& iEvent, double ParticleGenPt, double ParticleGenEta, double ParticleGenPhi, double iso_DR, double iso_Seed);
   // Generic function to write ParticleFlow based isolation into (PXL) photons and
   // electrons. Could be extended to other particles as well.
   //
   template< typename T > void particleFlowBasedIsolation( const edm::Event &iEvent,
                                                           const std::vector< edm::InputTag > &inputTagIsoValPFId,
                                                           const edm::Ref< T > &ref,
                                                           pxl::Particle &part,
                                                           const bool &useIsolator = true ) const;

   // ----------member data ---------------------------

   int fNumEvt;// used to count the number of events
   int fDebug; 
   std::string fFileName; 
   std::string fProcess;
   bool fGenOnly;
   bool fUseSIM;
   std::string fLHgridName;
   int fNumLHgridErrorSets;
   // The labels used in cfg-file 
   // Generator 
   std::string fgenParticleCandidatesLabel;
   std::string fMETMCLabel;
   double m_rhoFastJet;
   std::string fVertexRecoLabel;
   //Tau
   std::string fTauRecoLabel;
   std::string fPFTauDiscriminator;
   // Muon
   std::string fMuonRecoLabel;
   // Electron
   std::string fElectronRecoLabel;
   // GSF Electrons for vetoing.
   edm::InputTag m_gsfElectronsTag;
   // for PF isolation
   typedef std::vector< edm::Handle< edm::ValueMap< double > > > IsoDepositVals;
   std::vector< edm::InputTag > m_inputTagIsoValElectronsPFId;
   std::string m_eleEffAreaTargetLabel;
   // Photon
   std::string fGammaRecoLabel;
   // for PF isolation
   std::vector< edm::InputTag > m_inputTagIsoValPhotonsPFId;
   edm::InputTag m_particleFlowTag;
   // Jets
   std::vector< jet_def > jet_infos;
   //JetIDs
   typedef std::vector< std::pair< std::string, Selector<pat::Jet>* > > jet_id_list;
   // MET labels
   std::vector< collection_def > MET_infos;
   // Cluster
   edm::InputTag freducedBarrelRecHitCollection;
   edm::InputTag freducedEndcapRecHitCollection;
   //HCAL noise
   edm::InputTag hcal_noise_label;

   // Conversions for vetoing.
   edm::InputTag m_conversionsTag;

   //all triggers
   std::vector< trigger_group > triggers;
   std::vector< trigger_group > filters;

   bool fStoreL3Objects;
 
   ParticleMatcher* Matcher;
   // to be used for ePax output 
   pxl::OutputFile fePaxFile;
   std::vector<gen::PdfInfo> fpdf_vec;
   double xfx(const double &x, const double &Q, int fl) {  
      double f[13], mx = x, mQ = Q;
      evolvepdf_(mx, mQ, f);
      return f[fl+6];
   };

   //cuts
   double min_muon_pt,
      min_ele_pt,
      min_gamma_pt,
      min_jet_pt,
      min_met,
      min_tau_pt,
      max_eta,
      min_rechit_energy,
      min_rechit_swiss_cross,
      min_rechit_R19,
      vertex_minNDOF,
      vertex_maxZ,
      vertex_maxR,
      PV_minNDOF,
      PV_maxZ,
      PV_maxR;

   //vertex for physics eta, phi, pt
   reco::BeamSpot::Point the_vertex;
   reco::BeamSpot::Point the_beamspot;
};
#endif
