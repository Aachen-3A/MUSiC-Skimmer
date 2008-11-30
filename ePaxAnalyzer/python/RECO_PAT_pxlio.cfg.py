import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/023A83CF-C1B1-DD11-9FCE-001F2908CF18.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/040AD322-BBB1-DD11-8173-001CC47C50D6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/0A29FD41-BBB1-DD11-A503-001F29087EFC.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/124D1436-CAB1-DD11-95CD-001F29089F2A.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/12C89773-CAB1-DD11-92F8-001CC47BCFDC.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/143E7A50-BBB1-DD11-A223-001CC4BDF336.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/1A79E058-CAB1-DD11-AE3E-001CC4445076.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/241843CB-C9B1-DD11-B8FA-001CC47BCDA6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/265891D1-C9B1-DD11-89F4-001E0B466C58.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/266AA095-BCB1-DD11-B868-001E0B469CAA.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/28956F89-BDB1-DD11-959B-001F2908CE58.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/2898C04B-CAB1-DD11-B6C4-001F29088E72.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/2A358C3B-CAB1-DD11-BF46-001F29078D4C.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/2ADB6368-CAB1-DD11-AE98-001F29087E0A.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/2CCAEF68-CAB1-DD11-A1C5-001CC47D43D4.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/2CFB9AB9-CAB1-DD11-B713-001CC47BCDA6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/2E2C87C3-C9B1-DD11-9639-001CC47B2ED2.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/320254D1-C9B1-DD11-A0BD-001E0B47B5E4.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/342C5ED2-C9B1-DD11-80F0-001E0B4746F0.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/3CAD2F2C-CAB1-DD11-9CB8-001F2908BE30.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/3E378EBC-DBB1-DD11-8A5C-001F2908CE58.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/42C8B16D-DBB1-DD11-9523-001F29097068.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/4448C2AA-CAB1-DD11-B03D-001CC4A7C0A4.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/46E3744F-BDB1-DD11-9336-001CC4BD4A46.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/569504D6-C1B1-DD11-B161-001E0BEC51FE.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/58E0BC4A-BBB1-DD11-98BA-001CC47D037C.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/5E420F54-BBB1-DD11-BC63-001CC4C1FECA.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/5E9C8F54-DBB1-DD11-A2E6-001F29086E74.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/6079AB1F-BBB1-DD11-85D2-001F2908CFBC.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/622D6069-CAB1-DD11-B2F6-001CC47D420E.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/6282F57E-CAB1-DD11-88D9-001CC416D644.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/64E71052-CAB1-DD11-A549-001CC47D8FBA.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/68A53854-BBB1-DD11-9AF2-001E0B46C9A0.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/6C357132-CAB1-DD11-8C2C-001CC47D8D40.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/6EDB86CA-C9B1-DD11-B367-001CC411CFC0.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/74CBA64F-BBB1-DD11-8A91-001E0BED1522.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/82A7AB5B-CAB1-DD11-BCAE-001CC4A6CD60.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/86594DBE-C9B1-DD11-830E-001CC47C0158.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/865E087B-CAB1-DD11-9008-001F2907EF7E.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/8AFA39C3-C9B1-DD11-AD62-001CC47D5DA2.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/8E54141C-BBB1-DD11-B9FD-001CC47B30A0.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/8E984BAA-CAB1-DD11-964B-001CC47D5DA2.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/964779B1-CAB1-DD11-B1A9-001F29079F98.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/989B793D-BBB1-DD11-A107-001CC47BF09C.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/9AACE531-BBB1-DD11-9713-001CC47C9040.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/9AD429CE-C1B1-DD11-B8FC-001CC4A68C80.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/A243B30C-DDB1-DD11-B9D7-001F2908AED8.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/A45680B0-CAB1-DD11-8217-001E0BEC51FE.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/A4913944-CAB1-DD11-8A20-001F2908AED8.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/A6329B49-CAB1-DD11-8363-001F2908BE72.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/AA63FAAF-CAB1-DD11-9C10-001F2907EE22.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/AC344AA1-CAB1-DD11-B251-001CC4A63C8E.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/AC68EEBC-DBB1-DD11-A1F4-001F2907EF28.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/AE8031CF-C1B1-DD11-8636-001F290860A6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/B0605D57-DBB1-DD11-AC33-001CC416D644.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/B45E77B0-CAB1-DD11-83B2-001F29082E9E.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/B622C83F-CAB1-DD11-8696-001F29082E8C.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/B6F432B2-CAB1-DD11-BE63-001F290789D6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/B8BD70DB-CAB1-DD11-A3D7-001F29086E20.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/BAC51C2B-BBB1-DD11-BAF5-001F29081EBA.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/BCBA9DD9-C1B1-DD11-AE63-001F29087E0A.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/BE16A1C7-DBB1-DD11-8C6D-001E0B486068.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/BEEB4DF4-DCB1-DD11-B8BB-001F2908AF72.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/C2B3A7A4-CAB1-DD11-8A2C-001F29081EB6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/C8260DC5-C1B1-DD11-9373-001F29085CF2.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/CAA69CD3-C9B1-DD11-BD9E-001E0B465712.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/CCC3F735-BBB1-DD11-946D-001CC47C90A6.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/D2F646B2-CAB1-DD11-814E-001F29087EFC.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/D6534D32-BBB1-DD11-9C4C-001F2908BEFA.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/D8884CC5-C9B1-DD11-95D2-001CC4BE3084.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/D8CF97A3-DBB1-DD11-94BF-001F29082E9E.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/DC8BEA1F-CAB1-DD11-A36E-001CC47D5BD2.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/DEE9AF63-CAB1-DD11-BFB6-001E0B462E8A.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/E04B3359-BBB1-DD11-9A1B-001CC47D037C.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/EAAB1B25-CAB1-DD11-89FE-001CC4A6CC90.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/EE032169-BBB1-DD11-A50D-001CC47D43D4.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/F2E565AF-DBB1-DD11-8B07-001F2908CE4A.root',
       '/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/FC0544C8-C9B1-DD11-B996-001CC47C90A6.root'
	)
)

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


import os
cmsbase = os.environ.get('CMSSW_BASE')
execfile(cmsbase + "/src/ePaxDemo/ePaxAnalyzer/python/configurePAT_cff")

process.ePaxAnalysis = cms.EDAnalyzer("ePaxAnalyzer",
         # label of file:
         FileName =  cms.untracked.string("PAT_LM3.pxlio"),
         # Debugging: 0 = off, 1 = human readable, 2 = insane
         debug = cms.untracked.int32(0),
         Process = cms.untracked.string("LM3"),
         # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
	 GenOnly = cms.untracked.bool(False),
         #labels of source
         genParticleCandidatesLabel = cms.untracked.string("genParticles"),
         METMCLabel = cms.untracked.string("genMetNoNuBSM"),  # muon-correction needed?!??! ---> no!
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

process.p = cms.Path(process.patLayer0 + process.patLayer1 + process.ePaxAnalysis)

## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)


