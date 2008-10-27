import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/04419036-F385-DD11-B3A7-001617C3B6E8.root'
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/0A28F869-F285-DD11-AF3C-001617DBD5B2.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/162C4B5E-F585-DD11-872A-001617C3B64C.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/205E6CE3-F485-DD11-9D53-001617C3B76A.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/562BAFA1-F585-DD11-B931-001617DBD224.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/565AFE10-EF85-DD11-8353-000423D6B42C.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/5C66302A-F185-DD11-81D3-000423D98834.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/66F60641-F685-DD11-A493-000423D987FC.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/6E6A6E2D-F485-DD11-B707-001617DBD472.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/70191C8D-F485-DD11-8280-001617E30D06.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/9A204B65-F385-DD11-9CF1-000423D98B6C.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/A6E8BAB0-F085-DD11-9AB1-000423D986C4.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/B4D8FA14-F485-DD11-A41D-001617C3B76A.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/B6CA9FAB-F185-DD11-B66B-001617E30D0A.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/D8DAE3BF-F385-DD11-9C8E-001617C3B65A.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/EA02E08F-F285-DD11-8AF3-000423D9870C.root',
#'/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0001/6E9B44E2-0487-DD11-BFA7-001617C3B78C.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_04419036-F385-DD11-B3A7-001617C3B6E8_2.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_0A28F869-F285-DD11-AF3C-001617DBD5B2.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_162C4B5E-F585-DD11-872A-001617C3B64C.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_205E6CE3-F485-DD11-9D53-001617C3B76A.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_562BAFA1-F585-DD11-B931-001617DBD224.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_565AFE10-EF85-DD11-8353-000423D6B42C.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_5C66302A-F185-DD11-81D3-000423D98834.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_66F60641-F685-DD11-A493-000423D987FC.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_6E6A6E2D-F485-DD11-B707-001617DBD472.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_mm_2_1_9_relval_0A249693-FC85-DD11-AE0A-000423D99896.root'
#'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/16AAC418-218A-DD11-AC33-001F2908F0E4.root'

)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.Geometry_cff")

# for ClusterShape inspired by egamma hypernews
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

# Remove unneccessary stuff:
process.patLayer1.remove(process.layer1Hemispheres)


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
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel='Icone5',doType1MET=True,doL1Counters=False)
addJetCollection(process,'sisCone5CaloJets','SISC5',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel='Scone5',doType1MET=True,doL1Counters=False)
addJetCollection(process,'sisCone7CaloJets','SISC7',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel='Scone7',doType1MET=True,doL1Counters=False)
addJetCollection(process,'kt4CaloJets','KT4',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel='FKt4',doType1MET=True,doL1Counters=False)
addJetCollection(process,'kt6CaloJets','KT6',
                        runCleaner="CaloJet",doJTA=True,doBTagging=True,jetCorrLabel='FKt6',doType1MET=True,doL1Counters=False)

process.ePaxAnalysis = cms.EDAnalyzer("ePaxAnalyzer",
         # label of file:
         FileName =  cms.untracked.string("PAT_Zmumu.pxlio"),
         # Debugging: 0 = off, 1 = human readable, 2 = insane
         debug = cms.untracked.int32(0),
         Process = cms.untracked.string("Zmumu"),
         # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
	 GenOnly = cms.untracked.bool(False),
         #labels of source
         genParticleCandidatesLabel = cms.untracked.string("genParticles"),
         METMCLabel = cms.untracked.string("genMetNoNuBSM"),  # muon-correction needed?!??!
	 #FIXME make sure that this is the correct Collection! (BS = with beam spot constraints?)
         VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
         MuonRecoLabel = cms.untracked.string("selectedLayer1Muons"),
         ElectronRecoLabel = cms.untracked.string("selectedLayer1Electrons"),
         GammaRecoLabel = cms.untracked.string("selectedLayer1Photons"),
         # Jet labels: used for Gen AND REC Jets
	 JetMCLabels = cms.vstring("kt4GenJets","kt6GenJets","sisCone5GenJets", "sisCone7GenJets", "iterativeCone5GenJets"),
	 JetRecoLabels = cms.vstring("KT4", "KT6", "SISC5", "SISC7", "IC5"),
         # MET
         METRecoLabel = cms.untracked.string("selectedLayer1METs"),
	 reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    	 reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
     	 barrelClusterCollection = cms.InputTag("correctedHybridSuperClusters","electronPixelSeeds"),
    	 endcapClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","electronPixelSeeds"),
	 triggerResults = cms.InputTag("TriggerResults","","HLT"),
         triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")

)

process.p = cms.Path(process.patLayer0 + process.patLayer1 + process.ePaxAnalysis)


