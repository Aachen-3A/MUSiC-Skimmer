runOnData = True
#run on Summer09 (ReReco'd or not)
runOnSummer09 = False
#run on ReReco'ed data or Summer09 MC
runOnReReco = False

if runOnData and runOnSummer09:
    print "runOnData and runOnSummer09 can't be true at the same time!"
    import sys
    sys.exit(1)

if runOnSummer09 and not runOnReReco:
    print 'Reco of CMSSW < 3.5.X is not supported anymore!'
    import sys
    sys.exit(1)


import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.limit = 100

import FWCore.Framework.test.cmsExceptionsFatalOption_cff
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    #open file in NOMERGE mode to avoid a memory leak
    fileMode = cms.untracked.string( 'NOMERGE' ),
    #stop processing on each and every thrown exception
    Rethrow = FWCore.Framework.test.cmsExceptionsFatalOption_cff.Rethrow
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
'/store/data/Commissioning10/MinimumBias/RAW-RECO/May6thPDSkim_Skim_logerror-v1/0130/0A06F34A-5A5D-DF11-8CCC-0018F3D095EC.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/May6thPDSkim_Skim_logerror-v1/0123/60525751-F35C-DF11-A397-0018F3D095EC.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/May6thPDSkim_Skim_logerror-v1/0122/CABC72FE-C95C-DF11-A944-001A92810AC8.root'
	)
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if runOnData:
    if runOnReReco:
        process.GlobalTag.globaltag = cms.string('GR_R_37X_V6::All')
    else:
        process.GlobalTag.globaltag = cms.string('GR_R_37X_V6::All')
else:
    process.GlobalTag.globaltag = cms.string('START37_V5::All')

process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")


import MUSiCProject.Skimming.Tools
MUSiCProject.Skimming.Tools.configurePAT( process, runOnData, runOnReReco, runOnSummer09 )


#filter on right BX in case of data
if runOnData:
    process.primaryVertexFilter = cms.EDFilter( "VertexSelector",
                                                src = cms.InputTag( "offlinePrimaryVertices" ),
                                                cut = cms.string( "!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2" ),
                                                filter = cms.bool( True ), # otherwise it won't filter the events, just produce an empty vertex collection.
                                                )                            

    process.scrapingFilter = cms.EDFilter( "FilterOutScraping",
                                           applyfilter = cms.untracked.bool( True ),
                                           debugOn = cms.untracked.bool( False ),
                                           numtrack = cms.untracked.uint32( 10 ),
                                           thresh = cms.untracked.double( 0.25 )
                                           )


    process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
    process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
    process.technicals = process.hltLevel1GTSeed.clone()
    process.technicals.L1TechTriggerSeeding = cms.bool(True)
    process.technicals.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ( (42 AND NOT 43) OR (43 AND NOT 42) )')

    process.incompleteECALFilter = cms.EDFilter( "recHitFilter" )

    process.p = cms.Path( process.primaryVertexFilter * process.scrapingFilter * process.technicals * process.incompleteECALFilter * process.patDefaultSequence )

else:
    process.p = cms.Path( process.patDefaultSequence )


process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )

if runOnSummer09:
    #add the high-lumi trigger
    process.Skimmer.triggers.HLT1E31 = cms.PSet(
        process = cms.string('HLT'),
        L1_result = cms.InputTag( "gtDigis" ),
        results = cms.string('TriggerResults'),
        event   = cms.string('hltTriggerSummaryAOD'),
        HLTriggers = cms.vstring(
            'HLT_Mu9', 'HLT_DoubleMu0',
            'HLT_Ele20_SW_L1R', 'HLT_DoubleEle10_SW_L1R',
            'HLT_Photon25_L1R', 'HLT_DoublePhoton15_L1R',
            'HLT_L1Mu14_L1SingleEG10', 'HLT_L2Mu5_Photon9_L1R'
        )
    )
    
    #use the re-digi trigger in re-reco
    process.Skimmer.triggers.HLT.process = 'REDIGI'


if not runOnData:
    MUSiCProject.Skimming.Tools.addFlavourMatching( process, process.Skimmer, process.p )

process.p += process.Skimmer
