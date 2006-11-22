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
   virtual void analyzeGenInfo(const edm::Event&, ePax::ePaxEventViewRef);
   virtual void analyzeGenJets(const edm::Event&, ePax::ePaxEventViewRef);
   virtual void analyzeGenMET(const edm::Event&, ePax::ePaxEventViewRef);

   virtual void analyzeRecMuons(const edm::Event&, ePax::ePaxEventViewRef);
   virtual void analyzeRecElectrons(const edm::Event&, ePax::ePaxEventViewRef);
   virtual void analyzeRecJets(const edm::Event&, ePax::ePaxEventViewRef);
   virtual void analyzeRecMET(const edm::Event&, ePax::ePaxEventViewRef);
   virtual void analyzeRecGammas(const edm::Event&, ePax::ePaxEventViewRef);

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
   // The labels used in cfg-file 
   std::string fHepMCLabel;
   std::string fJetMCLabel;
   std::string fMETMCLabel;
   std::string fMuonRecoLabel;
   std::string fElectronRecoLabel;
   std::string fGammaRecoLabel;
   std::string fJetRecoLabel;
   std::string fMETRecoLabel;
    
  
   // to be used for ePax output 
   iotl::oDiskFile fePaxFile;

   //maximal number of physical objects, check that not too small!!!
/*   static const int kMaxGen = 100;
   static const int kMaxRec = 100;
   static const int kMaxJets = 100;
   static const int kMaxMET = 1;

   int fNumMuonMC;
   double fMuonMC_Px[kMaxGen];
   double fMuonMC_Py[kMaxGen];
   double fMuonMC_Pz[kMaxGen];
   double fMuonMC_Pt[kMaxGen];
   double fMuonMC_E[kMaxGen];
   double fMuonMC_Phi[kMaxGen];
   double fMuonMC_Theta[kMaxGen];
   double fMuonMC_Eta[kMaxGen];
   int    fMuonMC_Charge[kMaxGen];
   double fMuonMC_Vtx_X[kMaxGen];
   double fMuonMC_Vtx_Y[kMaxGen];
   double fMuonMC_Vtx_Z[kMaxGen];

   int fNumEleMC;
   double fEleMC_Px[kMaxGen];
   double fEleMC_Py[kMaxGen];
   double fEleMC_Pz[kMaxGen];
   double fEleMC_Pt[kMaxGen];
   double fEleMC_E[kMaxGen];
   double fEleMC_Phi[kMaxGen];
   double fEleMC_Theta[kMaxGen];
   double fEleMC_Eta[kMaxGen];
   int    fEleMC_Charge[kMaxGen];
   double fEleMC_Vtx_X[kMaxGen];
   double fEleMC_Vtx_Y[kMaxGen];
   double fEleMC_Vtx_Z[kMaxGen];
   
   int fNumGammaMC;
   double fGammaMC_Px[kMaxGen];
   double fGammaMC_Py[kMaxGen];
   double fGammaMC_Pz[kMaxGen];
   double fGammaMC_Pt[kMaxGen];
   double fGammaMC_E[kMaxGen];
   double fGammaMC_Phi[kMaxGen];
   double fGammaMC_Theta[kMaxGen];
   double fGammaMC_Eta[kMaxGen];
   
   int fNumMETMC;
   double fMETMC_SumET[kMaxMET];
   double fMETMC_METSig[kMaxMET];
   double fMETMC_Ez[kMaxMET];
   double fMETMC_MET[kMaxMET];
   double fMETMC_MEx[kMaxMET];
   double fMETMC_MEy[kMaxMET];
   double fMETMC_Phi[kMaxMET];
   double fMETMC_EmEnergy[kMaxMET];
   double fMETMC_HadEnergy[kMaxMET];

   int fNumKtJetMC;
   double fKtJetMC_Px[kMaxJets];
   double fKtJetMC_Py[kMaxJets];
   double fKtJetMC_Pz[kMaxJets];
   double fKtJetMC_Pt[kMaxJets];
   double fKtJetMC_E[kMaxJets];
   double fKtJetMC_Phi[kMaxJets];
   double fKtJetMC_Theta[kMaxJets];
   double fKtJetMC_Eta[kMaxJets];

   int fNumMuonRec;
   double fMuonRec_Px[kMaxRec];
   double fMuonRec_Py[kMaxRec];
   double fMuonRec_Pz[kMaxRec];
   double fMuonRec_Pt[kMaxRec];
   double fMuonRec_E[kMaxRec];
   double fMuonRec_Phi[kMaxRec];
   double fMuonRec_Theta[kMaxRec];
   double fMuonRec_Eta[kMaxRec];
   int    fMuonRec_Charge[kMaxRec];
   double fMuonRec_Vtx_X[kMaxRec];
   double fMuonRec_Vtx_Y[kMaxRec];
   double fMuonRec_Vtx_Z[kMaxRec];

   int fNumEleRec;
   double fEleRec_Px[kMaxRec];
   double fEleRec_Py[kMaxRec];
   double fEleRec_Pz[kMaxRec];
   double fEleRec_Pt[kMaxRec];
   double fEleRec_E[kMaxRec];
   double fEleRec_Phi[kMaxRec]; 
   double fEleRec_Theta[kMaxRec];
   double fEleRec_Eta[kMaxRec];
   int    fEleRec_Charge[kMaxRec];
   double fEleRec_Vtx_X[kMaxRec];
   double fEleRec_Vtx_Y[kMaxRec];
   double fEleRec_Vtx_Z[kMaxRec];

   int fNumMETRec;
   double fMETRec_SumET[kMaxMET];
   double fMETRec_METSig[kMaxMET];
   double fMETRec_Ez[kMaxMET];
   double fMETRec_MET[kMaxMET];
   double fMETRec_MEx[kMaxMET];
   double fMETRec_MEy[kMaxMET];
   double fMETRec_Phi[kMaxMET];
  
   int fNumGammaRec;
   double fGammaRec_Px[kMaxRec];
   double fGammaRec_Py[kMaxRec];
   double fGammaRec_Pz[kMaxRec];
   double fGammaRec_Pt[kMaxRec];
   double fGammaRec_E[kMaxRec];
   double fGammaRec_Phi[kMaxRec];
   double fGammaRec_Theta[kMaxRec];
   double fGammaRec_Eta[kMaxRec];

   int fNumKtJetRec;
   double fKtJetRec_Px[kMaxJets];
   double fKtJetRec_Py[kMaxJets];
   double fKtJetRec_Pz[kMaxJets];
   double fKtJetRec_Pt[kMaxJets];
   double fKtJetRec_E[kMaxJets]; 
   double fKtJetRec_Phi[kMaxJets];
   double fKtJetRec_Theta[kMaxJets];
   double fKtJetRec_Eta[kMaxJets];*/

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
