import FWCore.ParameterSet.Config as cms

process = cms.Process("Phythia8Photon35")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
<<<<<<< .mine
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/photon_jets_48878A79-CCBC-DD11-85F3-0022199A2E95.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/photon_jets_B89E5A6F-CFBC-DD11-888C-00E0814002A9.root'
#'/store/mc/Summer08/HerwigQCDPt300/GEN-SIM-RECO/IDEAL_V9_v1/0006/04189352-BFA5-DD11-A8A5-00D0680BF8C3.root'
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/16AAC418-218A-DD11-AC33-001F2908F0E4.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/1E19C1C2-EF89-DD11-A6AB-001E0B1C74DA.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/2AE099C8-1F8A-DD11-B30F-00144F2031D4.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/2CA4D6BE-EF89-DD11-AEFA-001F290860E4.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/32EE38A2-188A-DD11-827E-00144F283544.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/3AAE5B99-228A-DD11-8889-001CC47BF09C.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/3AB58659-178A-DD11-AFA1-00144F283E48.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/3CB730C2-0A8A-DD11-9CE2-001CC47D037C.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/3E58A651-178A-DD11-A33A-0015C5E9C186.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/4CD2835A-088A-DD11-8849-001CC47D037C.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/58D6B38D-EF89-DD11-B53B-001E0B475590.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/6AB2D1A1-2E8A-DD11-A65D-001CC445A5A8.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/722BA08C-EF89-DD11-8D2A-001CC4C0A4A4.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/72F144C2-EF89-DD11-9C58-001CC4C1FECA.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/74F721DD-EF89-DD11-B856-001E0BEC51FE.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/7679F11F-F089-DD11-95D8-001E0B469C96.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/7855DA8B-EF89-DD11-8F2B-001CC47C9040.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/7EA16907-1A8A-DD11-B0BD-001EC9DB2C1E.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/8466B531-248A-DD11-9B60-00E081402F3D.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/8CC7DE1F-218A-DD11-9CE3-001F29084160.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/9C95B62F-248A-DD11-AF7E-00E08140E719.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/BA331A99-188A-DD11-A11E-0015C5E9B30F.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/BC0A253A-078A-DD11-B22C-001CC47C7092.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/BEC43DE5-228A-DD11-9B39-001CC47BA216.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/C69E0D7B-EF89-DD11-A96C-001CC47D8FBA.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/CC176183-1E8A-DD11-B4EF-0015C5E9D0A6.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/DA45B4A7-2E8A-DD11-A743-001F2907EF7E.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/E232FAA0-188A-DD11-BB99-00144F2024D2.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/E28736D4-1C8A-DD11-944D-0015C5E9BEBD.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/E414E0A2-2E8A-DD11-9F95-001CC416D644.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/EA0B5A97-F089-DD11-9E54-001CC4A6CC32.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/F2A152B5-228A-DD11-9927-001F29072C1E.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/F65613C2-0A8A-DD11-8BE4-001CC47D43D4.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/F66338D3-1C8A-DD11-B242-0015C5E9BB60.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/FCC959CD-EF89-DD11-B167-001CC47D7FDE.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0004/FED1F3BF-EF89-DD11-9CBB-001E0B469778.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/10AFEFBB-268A-DD11-AA12-001CC47BCFDC.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/2C3F36ED-2F8A-DD11-96CD-001CC47C813C.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/3474D91E-288A-DD11-8845-001CC47D43D4.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/4238920E-358A-DD11-8054-001F29086E48.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/82BA9AFA-378A-DD11-9C9F-001E4F3D7654.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/A48BA2DA-2F8A-DD11-A97D-001CC445A5A8.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/BC5D0C20-288A-DD11-A4D9-001CC4A7D032.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/CED0ECE9-2F8A-DD11-925F-001F2907DA48.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/DC1CF760-278A-DD11-95B8-001CC47D2F90.root',
        #'/store/mc/Summer08/TauolaTTbar/GEN-SIM-RECO/IDEAL_V9_v1/0005/E0DA23B9-268A-DD11-9C2C-001CC47C813C.root'


=======
            #'/store/mc/Summer08/SUSY_LM4-sftsht/GEN-SIM-RECO/IDEAL_V9_v1/0000/FC0544C8-C9B1-DD11-B996-001CC47C90A6.root',
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/single_gamma_2_2_2_relval_1D08E7DEE-B4B9-DD11-B49E-001617E30D00.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Z_ee_2_1_9_relval_04419036-F385-DD11-B3A7-001617C3B6E8_2.root'
#'/store/mc/Summer08/HerwigQCDPt300/GEN-SIM-RECO/IDEAL_V9_v1/0006/04189352-BFA5-DD11-A8A5-00D0680BF8C3.root'
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/02FB229C-F9C9-DD11-98BB-001EC9AA92FC.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/04003EE7-04CA-DD11-8817-001C23C0B772.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/042187C8-04CA-DD11-A45E-001EC9AA9FE5.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/045527D4-33CA-DD11-9591-00304833457A.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/04BBF420-FFC9-DD11-B262-001D091C6771.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/061EC9A2-07CA-DD11-9860-0019B9D96FC5.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/064DE523-FFC9-DD11-B1C3-001EC9AA9FF4.root',
        '/store/mc/Fall08/WJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v1/0013/069222C9-04CA-DD11-AFF8-001EC9AAA058.root'
>>>>>>> .r634
	)
)

process.AdaptorConfig = cms.Service("AdaptorConfig",
    stats = cms.untracked.bool(True),
    enable = cms.untracked.bool(True),
    tempDir = cms.untracked.string(""),
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

process.p = cms.Path(process.patLayer0 + process.patLayer1 + process.ePaxAnalysis)

## Necessary fixes to run 2.2.X on 2.1.X data
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import run22XonSummer08AODSIM
run22XonSummer08AODSIM(process)


