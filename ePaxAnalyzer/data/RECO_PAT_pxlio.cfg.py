import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATLayer0Summary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default          = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    PATLayer0Summary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root'
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_1.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_2.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_3.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_4.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_5.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_6.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_7.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_8.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_9.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_10.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_11.root',
'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_12.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_13.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_14.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_15.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_16.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_17.root'
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/dcms/staschmitz/test/Zee_2_1_0_relval_18.root'
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.Geometry_cff")

#begin test for ClusterShape inspired by egamma hypernews
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#end test

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patLayer0_cff")
process.load("PhysicsTools.PatAlgos.patLayer1_cff")
#process.content = cms.EDAnalyzer("EventContentAnalyzer")
#process.p = cms.Path(
 #               process.patLayer0  
                #+ process.content # uncomment to get a dump of the output after layer 0
  #              + process.patLayer1  
   #         )

#process.p = cms.Path(
#                process.patLayer0  
                #+ process.content # uncomment to get a dump of the output after layer 0
#                + process.patLayer1  
#            )


# Output module configuration
#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
    # save only events passing the full path
#    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#    outputCommands = cms.untracked.vstring('drop *')
#)

process.ePaxAnalysis = cms.EDAnalyzer("ePaxAnalyzer",
# label of file:
         FileName =  cms.untracked.string("test_PAT_Zee_trig.pxlio"),
         # Debugging: 0 = off, 1 = human readable, 2 = insane
         debug = cms.untracked.int32(1),
         Process = cms.untracked.string("test_PAT_Zee"),
         # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
         #untracked bool GenOnly = false
         #labels of source
         #untracked string TruthVertexLabel = "trackingParticles"
         genParticleCandidatesLabel = cms.untracked.string("genParticles"),
         #untracked string KtJetMCLabel = "fastjet6GenJetsNoNuBSM"
         ItCone5JetMCLabel = cms.untracked.string("iterativeCone5GenJets"),
         METMCLabel = cms.untracked.string("genMetNoNuBSM"),  #here muon included thus no muon-correction needed
         #untracked string METMCLabel = "genMet"  
	 #FIXME make sure that this is the correct Collection! (BS = with beam spot constraints?)
         VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
         #untracked string SAMuonRecoLabel = "standAloneMuons"
	 # name of MuonCollection has changed!
         MuonRecoLabel = cms.untracked.string("selectedLayer1Muons"),
         ElectronRecoLabel = cms.untracked.string("selectedLayer1Electrons"),
         #untracked string PixelMatchElectronRecoLabel = "pixelMatchGsfElectrons"
         #string ElectronIDAssocProducer = "electronId" #this one is TRACKED!
         #string ElectronHcalIsolationProducer = "egammaTowerIsolation" #this one is TRACKED!
         #string ElectronEcalIsolationProducer = "egammaEcalIsolation" #this one is TRACKED!
         #string ElectronTrackIsolationProducer = "egammaElectronTkIsolation" #this one is TRACKED!
         #string ElectronTrackNumProducer = "egammaElectronTkNumIsolation" #this one is TRACKED!
        GammaRecoLabel = cms.untracked.string("selectedLayer1Photons"),
         #string GammaHcalIsolationProducer = "gammaTowerIsolation" #this one is TRACKED!
         #string GammaEcalIsolationProducer = "gammaEcalIsolation" #this one is TRACKED!
         #string GammaTrackIsolationProducer = "gammaPhotonTkIsolation" #this one is TRACKED!
         #string GammaTrackNumProducer = "gammaPhotonTkNumIsolation" #this one is TRACKED!
         #untracked string KtJetRecoLabel = "MCJetCorJetfastjet6"
        ItCone5JetRecoLabel = cms.untracked.string("selectedLayer1Jets"),
         #untracked string L2L3JESic5JetRecoLabel = "L3JetCorJetIcone5"
	 #"allLayer1METs" are apparently not written to root file!
         METRecoLabel = cms.untracked.string("selectedLayer1METs"),
         #untracked string METCorrRecoLabel = "corMetGlobalMuons" #save corrected MET as well
         #InputTag barrelClusterShapeAssociation = hybridSuperClusters:hybridShapeAssoc
         #InputTag endcapClusterShapeAssociation = islandBasicClusters:islandEndcapShapeAssoc
	 #endcapClusterShapeAssociation = cms.InputTag("hybridSuperClusters:islandEndcapShapeAssoc"), 
         #barrelClusterShapeAssociation = cms.InputTag("hybridSuperClusters:hybridShapeAssoc"),    
 	 #endcapClusterShapeAssociation = cms.InputTag("hybridSuperClusters","islandEndcapShapeAssoc"), 
         #barrelClusterShapeAssociation = cms.InputTag("hybridSuperClusters","hybridShapeAssoc")
	 reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    	 reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
     	 barrelClusterCollection = cms.InputTag("hybridSuperClusters:"),
    	 endcapClusterCollection = cms.InputTag("multi5x5BasicClusters:multi5x5EndcapBasicClusters"),
      	 TriggerResults = cms.untracked.string("hltTriggerSummaryAOD")
	 #untracked string fHBHELabel = "hbhereco"
         #untracked string fHBHEInstanceName = ""

)

process.p = cms.Path(
 	process.patLayer0  
  	+ process.patLayer1                  
	+ process.ePaxAnalysis
            )


