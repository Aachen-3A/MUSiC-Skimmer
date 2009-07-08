runOnData = False

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/antonius/test/photon_jets_48878A79-CCBC-DD11-85F3-0022199A2E95.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/antonius/test/photon_jets_B89E5A6F-CFBC-DD11-888C-00E0814002A9.root'
'/store/user/pieta/test_310/test_1.root'
#'/store/mc/Fall08/PhotonJets200toInf-madgraph/GEN-SIM-RECO/IDEAL_V9_reco-v1/0026/00462FCA-C5FC-DD11-9B76-001A9227D3D1.root'

#'file:/home/home1/institut_3a/dietzlaursonn/MUSiC/CMSSW_2_2_11/src/ePaxDemo/ePaxAnalyzer/python/0033A31F-C9FC-DD11-8989-001E8CCCE114.root',
#'file:/home/home1/institut_3a/dietzlaursonn/MUSiC/CMSSW_2_2_11/src/ePaxDemo/ePaxAnalyzer/python/00462FCA-C5FC-DD11-9B76-001A9227D3D1.root'

#'/store/mc/Summer08/Zee/GEN-SIM-RECO/IDEAL_V9_v1/0004/F291CFA2-C988-DD11-A691-001EC9DB312B.root',
#'/store/mc/Summer08/Zee/GEN-SIM-RECO/IDEAL_V9_v1/0004/F21F86A0-8B89-DD11-8ABF-00E08140E9C7.root'
	)
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

process.load('PhysicsTools.HepMCCandAlgos.genEventScale_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genEventWeight_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genEventPdfInfo_cfi')

# rerun photon ID
process.load("RecoEgamma.PhotonIdentification.photonId_cff")
# for ClusterShape inspired by egamma hypernews
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP31X_v1::All')
process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")
#fix missing PAT cleaning
process.patDefaultSequence = cms.Sequence(
    process.beforeLayer1Objects *
    process.allLayer1Objects *
    process.selectedLayer1Objects *
    process.cleanLayer1Objects
    )


# Remove unneccessary stuff:
if runOnData:
    process.patHighLevelReco_withoutPFTau.remove( process.patJetFlavourId )
    process.patDefaultSequence_withoutTrigMatch.remove( process.patMCTruth )
    process.allLayer1Muons.addGenMatch        = False
    process.allLayer1Jets.addGenPartonMatch   = False
    process.allLayer1Jets.addGenJetMatch      = False
    process.allLayer1Jets.getJetMCFlavour     = False
    process.allLayer1METs.addGenMET           = False



# this might be commented in in order to safe the edm root file containing the PAT Products
# Output module configuration
# process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
#    # save only events passing the full path
#    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#    outputCommands = cms.untracked.vstring('drop *')
# )
# process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
# from PhysicsTools.PatAlgos.patEventContent_cff import *
# process.out.outputCommands += patEventContent
# process.out.outputCommands.extend(["keep *_selectedLayer1Jets*_*_*"])

# insert following lines in order to get info about the HLT content (don't forget to include into path)
# import HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi
# process.hltAnalyzer = HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi.hltEventAnalyzerAOD.clone()

# add some other jet Collection (default one is: itCone5)
from PhysicsTools.PatAlgos.tools.jetTools import *
# which cleaner to use?!?
addJetCollection(process,cms.InputTag("iterativeCone5CaloJets"),'IC5',
                        doJTA=True,doBTagging=True,jetCorrLabel=('IC5','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,cms.InputTag("sisCone5CaloJets"),'SISC5',
                        doJTA=True,doBTagging=True,jetCorrLabel=('SC5','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,cms.InputTag("sisCone7CaloJets"),'SISC7',
                        doJTA=True,doBTagging=True,jetCorrLabel=('SC7','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,cms.InputTag("kt4CaloJets"),'KT4',
                        doJTA=True,doBTagging=True,jetCorrLabel=('KT4','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,cms.InputTag("kt6CaloJets"),'KT6',
                        doJTA=True,doBTagging=True,jetCorrLabel=('KT6','Calo'),doType1MET=True,doL1Counters=False)
#process.load("PhysicsTools.HepMCCandAlgos.genEventKTValue_cfi")

import os
cmsbase = os.environ.get('CMSSW_BASE')
execfile(cmsbase + "/src/ePaxDemo/ePaxAnalyzer/python/configurePAT_cff")

#process.pTHat = cms.EDFilter("PtHatFilter",
#         pt_hat_lower_bound = cms.double(200.),
#         pt_hat_upper_bound = cms.double(500.)
#)

# cut on the factorization scale e.g. suitable for cuts on inv. mass of resonances in Pythia!
#process.FacScale = cms.EDFilter("FacScaleFilter",
#         fac_scale_lower_bound = cms.double(200.),
#         fac_scale_upper_bound = cms.double(500.)
#)

process.ePaxAnalysis = cms.EDAnalyzer("ePaxAnalyzer",
                                      # label of file:
                                      FileName =  cms.untracked.string("PhotonJets200toInf-madgraph_NEWID.pxlio"),
                                      # Debugging: 0 = off, 1 = human readable, 2 = insane
                                      debug = cms.untracked.int32(0),
                                      Process = cms.untracked.string("PhotonJets200toInf"),
                                      # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
                                      GenOnly = cms.untracked.bool(False),
                                      # UseSIM true means to use SIM info for finding converted photons
                                      UseSIM = cms.untracked.bool(True),
                                      #labels of source
                                      genParticleCandidatesLabel = cms.untracked.string("genParticles"),
                                      METMCLabel = cms.untracked.string("genMetCalo"),  # muon-correction needed ---> yes!
                                      #FIXME make sure that this is the correct Collection! (BS = with beam spot constraints?)
                                      VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
                                      MuonRecoLabel = cms.untracked.string("cleanLayer1Muons"),
                                      ElectronRecoLabel = cms.untracked.string("cleanLayer1Electrons"),
                                      GammaRecoLabel = cms.untracked.string("cleanLayer1Photons"),
                                      # Jet labels: used for Gen AND REC Jets , order of used algorithms must be identical , first entry is used for matching
                                      JetMCLabels = cms.vstring("sisCone5GenJets", "kt4GenJets","kt6GenJets", "sisCone7GenJets", "iterativeCone5GenJets"),
                                      JetRecoLabels = cms.vstring( "SISC5" ,"KT4", "KT6", "SISC7", "IC5"),
                                      L1GlobalTriggerReadoutRecord = cms.InputTag("hltGtDigis"),
                                      #L1GlobalTriggerReadoutRecord = cms.InputTag("gtDigis"),
                                      L1TriggerObjectMapTag = cms.InputTag("hltL1GtObjectMap"),
                                      # MET
                                      METRecoLabel = cms.untracked.string("layer1METs"),
                                      reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                      reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                      barrelClusterCollection = cms.InputTag("correctedHybridSuperClusters","electronPixelSeeds"),
                                      endcapClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","electronPixelSeeds"),
                                      triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                      triggerEvent = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                      L1Triggers = cms.vstring( 'L1_SingleMuOpen', 'L1_SingleMu0', 'L1_SingleMu3', 'L1_SingleMu7', 'L1_SingleMu20', 'L1_DoubleMuOpen', 'L1_DoubleMu3',
                                                                'L1_SingleEG5', 'L1_SingleEG8', 'L1_DoubleEG5' ),
                                      CacheL1TriggerBits = cms.bool( True ),
                                      HLTriggers = cms.vstring( 'HLT_L1Mu20', 'HLT_L2Mu9', 'HLT_L2Mu11', 'HLT_Mu3', 'HLT_Mu5', 'HLT_Mu9', 'HLT_DoubleMu0', 'HLT_DoubleMu3', 'HLT_L1DoubleMuOpen', 'HLT_IsoMu3',
                                                                'HLT_Ele10_LW_L1R', 'HLT_Ele10_LW_EleId_L1R', 'HLT_Ele15_LW_L1R', 'HLT_Ele15_SiStrip_L1R', 'HLT_Ele15_SC10_LW_L1R', 'HLT_Ele20_LW_L1R', 'HLT_DoubleEle5_SW_L1R',
                                                                'HLT_Photon15_L1R', 'HLT_Photon15_TrackIso_L1R', 'HLT_Photon15_LooseEcalIso_L1R', 'HLT_Photon20_L1R', 'HLT_DoublePhoton5_eeRes_L1R', 'HLT_DoublePhoton5_Jpsi_L1R', 'HLT_DoublePhoton5_Upsilon_L1R', 'HLT_DoublePhoton10_L1R',
                                                                'HLT_L1Mu14_L1SingleEG10'
                                                                ),
                                      StoreL3Objects = cms.untracked.bool(False)
                                      )

process.p = cms.Path(process.genEventScale + process.genEventWeight + process.genEventPdfInfo + process.photonIDSequence + process.patDefaultSequence + process.ePaxAnalysis)
