runOnData = False

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/antonius/test/photon_jets_48878A79-CCBC-DD11-85F3-0022199A2E95.root'
#'/store/user/pieta/test_311/test_1.root'
#'/store/mc/Fall08/PhotonJets200toInf-madgraph/GEN-SIM-RECO/IDEAL_V9_reco-v1/0026/00462FCA-C5FC-DD11-9B76-001A9227D3D1.root'
'file:/opt/scratch/608D435D-7B70-DE11-B093-001CC4782AF8.root'
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
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')
process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")


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


import os
cmsbase = os.environ.get('CMSSW_BASE')
execfile(cmsbase + "/src/MUSiCProject/Skimming/python/configurePAT_cff")

#process.pTHat = cms.EDFilter("PtHatFilter",
#         pt_hat_lower_bound = cms.double(200.),
#         pt_hat_upper_bound = cms.double(500.)
#)

# cut on the factorization scale e.g. suitable for cuts on inv. mass of resonances in Pythia!
#process.FacScale = cms.EDFilter("FacScaleFilter",
#         fac_scale_lower_bound = cms.double(200.),
#         fac_scale_upper_bound = cms.double(500.)
#)

process.Skimmer = cms.EDAnalyzer("MUSiCSkimmer",
                                      # label of file:
                                      FileName =  cms.untracked.string("test_run.pxlio"),
                                      # Debugging: 0 = off, 1 = human readable, 2 = insane
                                      debug = cms.untracked.int32(0),
                                      Process = cms.untracked.string("test_run"),
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
                                      reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                      reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                      barrelClusterCollection = cms.InputTag("correctedHybridSuperClusters","electronPixelSeeds"),
                                      endcapClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","electronPixelSeeds"),
                                      # Jet labels: used for Gen AND REC Jets , order of used algorithms must be identical , first entry is used for matching
                                      JetMCLabels = cms.vstring("sisCone5GenJets", "sisCone7GenJets", "iterativeCone5GenJets"),
                                      JetRecoLabels = cms.vstring( "SISC5", "SISC7", "IC5"),
                                      # MET
                                      METRecoLabel = cms.untracked.string("layer1METs"),

                                      #trigger menu: 1E31
                                      triggerProcess = cms.string( 'HLT' ),
                                      L1GlobalTriggerReadoutRecord = cms.InputTag("gtDigis", "", "HLT"),
                                      L1TriggerObjectMapTag = cms.InputTag("hltL1GtObjectMap", "", "HLT"),
                                      L1Triggers = cms.vstring( 'L1_SingleMu7', 'L1_DoubleMuOpen',
                                                                'L1_SingleEG8', 'L1_DoubleEG5',
                                                                'L1_SingleMu14', 'L1_SingleEG10', 'L1_Mu3QE8_EG5'
                                                                ),
                                      triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                                      triggerEvent = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                                      HLTriggers = cms.vstring( 'HLT_Mu9', 'HLT_DoubleMu0',
                                                                'HLT_Ele20_SW_L1R', 'HLT_DoubleEle10_SW_L1R',
                                                                'HLT_Photon25_L1R', 'HLT_DoublePhoton15_L1R',
                                                                'HLT_L1Mu14_L1SingleEG10', 'HLT_L2Mu5_Photon9_L1R'
                                                                ),

                                      #trigger menu: 8E29
                                      triggerProcess2 = cms.string( 'HLT8E29' ),
                                      L1GlobalTriggerReadoutRecord2 = cms.InputTag("gtDigis", "", "HLT8E29"),
                                      L1TriggerObjectMapTag2 = cms.InputTag("hltL1GtObjectMap", "", "HLT8E29"),
                                      L1Triggers2 = cms.vstring( 'L1_SingleMuOpen', 'L1_SingleMu0', 'L1_SingleMu3', 'L1_DoubleMuOpen',
                                                                 'L1_SingleEG5', 'L1_SingleEG8', 'L1_DoubleEG5',
                                                                 'L1_SingleMu14', 'L1_SingleEG10', 'L1_Mu3QE8_EG5'
                                                                 ),
                                      triggerResults2 = cms.InputTag("TriggerResults", "", "HLT8E29"),
                                      triggerEvent2 = cms.InputTag("hltTriggerSummaryAOD", "", "HLT8E29"),
                                      HLTriggers2 = cms.vstring( 'HLT_Mu3', 'HLT_DoubleMu0',
                                                                 'HLT_Ele10_LW_L1R', 'HLT_DoubleEle5_SW_L1R',
                                                                 'HLT_Photon15_L1R', 'HLT_DoublePhoton10_L1R'
                                                                 #no usefull cross-channel trigger in this menu
                                                                 ),
                                      
                                      
                                      CacheL1TriggerBits = cms.bool( True ),
                                      StoreL3Objects = cms.untracked.bool(False)
                                      )

process.p = cms.Path(process.genEventScale + process.genEventWeight + process.genEventPdfInfo + process.photonIDSequence + process.patDefaultSequence + process.Skimmer)
