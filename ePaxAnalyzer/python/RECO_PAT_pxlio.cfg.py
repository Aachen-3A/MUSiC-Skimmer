runOnData = False

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/photon_jets_48878A79-CCBC-DD11-85F3-0022199A2E95.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/photon_jets_B89E5A6F-CFBC-DD11-888C-00E0814002A9.root'
        'file:/tmp/Exotica_Zee_M200.root'	)
)

process.AdaptorConfig = cms.Service("AdaptorConfig",
    stats = cms.untracked.bool(True),
    enable = cms.untracked.bool(True),
    tempDir = cms.untracked.string("."),
    cacheHint = cms.untracked.string("lazy-download"),
    readHint = cms.untracked.string("auto-detect"))

process.load("Configuration/StandardSequences/GeometryPilot2_cff")

# for ClusterShape inspired by egamma hypernews
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

# Remove unneccessary stuff:
process.patLayer1.remove(process.layer1Hemispheres)
if runOnData:
    process.patHighLevelReco_withoutPFTau.remove( process.patJetFlavourId )
    process.patLayer0_withoutTrigMatch.remove( process.patMCTruth )
    process.allLayer1Muons.addGenMatch        = False
    process.allLayer1Jets.addGenPartonMatch   = False
    process.allLayer1Jets.addGenJetMatch      = False
    process.allLayer1Jets.getJetMCFlavour     = False
    process.allLayer1METs.addGenMET           = False



# this might be commented in in order to safe the edm root file containing the PAT Products
# Output module configuration
#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
#    # save only events passing the full path
#    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#    outputCommands = cms.untracked.vstring('drop *')
#)
#process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
#process.load("PhysicsTools.PatAlgos.patLayer1_EventContent_cff")
#process.out.outputCommands.extend(process.patLayer1EventContent.outputCommands)
#process.out.outputCommands.extend(["keep *_selectedLayer1Jets*_*_*"])

#insert following lines in order to get info about the HLT content (don't forget to include into path)
#import HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi
#process.hltAnalyzer = HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi.hltEventAnalyzerAOD.clone()

# add some other jet Collection (default one is: itCone5)
from PhysicsTools.PatAlgos.tools.jetTools import *
# which cleaner to use?!?
addJetCollection(process,'iterativeCone5CaloJets','IC5',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel=('IC5','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'sisCone5CaloJets','SISC5',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel=('SC5','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'sisCone7CaloJets','SISC7',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel=('SC7','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'kt4CaloJets','KT4',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel=('KT4','Calo'),doType1MET=True,doL1Counters=False)
addJetCollection(process,'kt6CaloJets','KT6',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel=('KT6','Calo'),doType1MET=True,doL1Counters=False)

#process.load("PhysicsTools.HepMCCandAlgos.genEventKTValue_cfi")

import os
cmsbase = os.environ.get('CMSSW_BASE')
execfile(cmsbase + "/src/ePaxDemo/ePaxAnalyzer/python/configurePAT_cff")

#process.pTHat = cms.EDFilter("PtHatFilter",
#         pt_hat_lower_bound = cms.double(200.),
#         pt_hat_upper_bound = cms.double(500.)
#)

# cut on the factorization scale e.g. suitable for cuts on inv. mass of resonances in Pythia!
process.FacScale = cms.EDFilter("FacScaleFilter",
         fac_scale_lower_bound = cms.double(200.),
         fac_scale_upper_bound = cms.double(500.)
)

process.ePaxAnalysis = cms.EDAnalyzer("ePaxAnalyzer",
         # label of file:
         FileName =  cms.untracked.string("testSIMphoton_temp.pxlio"),
         # Debugging: 0 = off, 1 = human readable, 2 = insane
         debug = cms.untracked.int32(0),
         Process = cms.untracked.string("Phythia8Photon35"),
         # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
	 GenOnly = cms.untracked.bool(False),
         # UseSIM true means to use SIM info for finding converted photons
         UseSIM = cms.untracked.bool(True),
         #labels of source
         genParticleCandidatesLabel = cms.untracked.string("genParticles"),
         METMCLabel = cms.untracked.string("genMet"),  # muon-correction needed ---> yes!
	 #FIXME make sure that this is the correct Collection! (BS = with beam spot constraints?)
         VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
         MuonRecoLabel = cms.untracked.string("selectedLayer1Muons"),
         ElectronRecoLabel = cms.untracked.string("selectedLayer1Electrons"),
         GammaRecoLabel = cms.untracked.string("selectedLayer1Photons"),
         # Jet labels: used for Gen AND REC Jets , order of used algorithms must be identical , first entry is used for matching
	 JetMCLabels = cms.vstring("sisCone5GenJets", "kt4GenJets","kt6GenJets", "sisCone7GenJets", "iterativeCone5GenJets"),
	 JetRecoLabels = cms.vstring( "SISC5" ,"KT4", "KT6", "SISC7", "IC5"),
	 L1GlobalTriggerReadoutRecord = cms.InputTag("hltGtDigis"),
	 #L1GlobalTriggerReadoutRecord = cms.InputTag("gtDigis"),
	 L1TriggerObjectMapTag = cms.InputTag("hltL1GtObjectMap"),
         # MET
         METRecoLabel = cms.untracked.string("selectedLayer1METs"),
	 reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    	 reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
     	 barrelClusterCollection = cms.InputTag("correctedHybridSuperClusters","electronPixelSeeds"),
    	 endcapClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","electronPixelSeeds"),
	 triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
	 triggerEvent = cms.InputTag("hltTriggerSummaryAOD", "", "HLT")

)

process.p = cms.Path(process.patLayer0 + process.patLayer1 + process.FacScale + process.ePaxAnalysis)

## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)


