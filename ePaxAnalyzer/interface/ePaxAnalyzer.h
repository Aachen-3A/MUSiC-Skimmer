#ifndef ePaxAnalyzer_H
#define ePaxAnalyzer_H
//
// class decleration
//

// CMSSW includes
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
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
// ROOT stuff

#include "TFile.h"
#include "TTree.h"

// ePax stuff
// Has to be included as the last header otherwise there will be a warning concerning the 
// zlib. According to Steffen there are two different zlib and ROOT can only deal with one of them
// but ePax can deal with both of them
#include "ePaxPxl/ePax/interface/ePax.h"

class ePaxAnalyzer : public edm::EDAnalyzer {
public:

   // why explicit?
   explicit ePaxAnalyzer(const edm::ParameterSet&);
   ~ePaxAnalyzer();
  

private:

   virtual void beginJob(const edm::EventSetup&);
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();
   virtual void analyzeGenInfo(const edm::Event&, ePaxEventViewRef);
   virtual void analyzeGenJets(const edm::Event&, ePaxEventViewRef);
   virtual void analyzeGenMET(const edm::Event&, ePaxEventViewRef);

   virtual void analyzeRecMuons(const edm::Event&, ePaxEventViewRef);
   virtual void analyzeRecElectrons(const edm::Event&, ePaxEventViewRef);
   virtual void analyzeRecJets(const edm::Event&, ePaxEventViewRef);
   virtual void analyzeRecMET(const edm::Event&, ePaxEventViewRef);
   virtual void analyzeRecGammas(const edm::Event&, ePaxEventViewRef);

   bool MuonMC_cuts(HepMC::GenEvent::particle_const_iterator MCmuon) const;
   bool EleMC_cuts(HepMC::GenEvent::particle_const_iterator MCele) const;
   bool GammaMC_cuts(HepMC::GenEvent::particle_const_iterator MCgamma) const;
   bool JetMC_cuts(reco::GenJetCollection::const_iterator MCjet) const;
   bool METMC_cuts(const reco::GenMET MCmet) const;
   bool Muon_cuts(reco::MuonCollection::const_iterator muon) const;
   bool Ele_cuts(reco::ElectronCollection::const_iterator ele) const;
   bool Gamma_cuts(reco::PhotonCollection::const_iterator photon) const;
   bool Jet_cuts(reco::CaloJetCollection::const_iterator jet) const;
   bool MET_cuts(const reco::MET met) const;
  
   // ----------member data ---------------------------

   int fNumEvt;// used to count the number of events
   int fDebug; 
   std::string fFileName; 
   // The labels used in cfg-file 
   std::string fHepMCLabel;
   std::string fKtJetMCLabel;
   std::string fItCone5JetMCLabel;
   std::string fMidCone5JetMCLabel;
   std::string fMidCone7JetMCLabel;   
   std::string fMETMCLabel;
   std::string fMuonRecoLabel;
   std::string fSAMuonRecoLabel;
   std::string fElectronRecoLabel;
   std::string fPixelMatchElectronRecoLabel;
   std::string fGammaRecoLabel;
   std::string fKtJetRecoLabel;
   std::string fItCone5JetRecoLabel;
   std::string fMidCone5JetRecoLabel;
   std::string fMidCone7JetRecoLabel;
   std::string fMETRecoLabel;
    
   // to be used for ePax output 
   pxl::oDiskFile fePaxFile;

};
#endif


/*

ToDO:
- Gleich  Objecte *sortiert* pt? in die Arrays einfuellen
- Check charge convention: is pid() == 13 Mu- oder Mu+ ??
- check angles (phi, theta for correct range and referenz)
- Put product lables into cfg ..
- other Jets
- primary Vertex
- Angkes ...
- TRIGGER?
- what kind of products are available? 
- Are Muons, Electrons, ... always != 0-Pointer?!?   CXheck needed???
*/
