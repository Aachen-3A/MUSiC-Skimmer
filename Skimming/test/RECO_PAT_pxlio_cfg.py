runOnData = True
#run on ReReco'ed data or Summer09 MC
runOnReReco = False
#run on NOT ReReco'ed Summer09
runOnSummer09 = False

if runOnReReco and runOnSummer09:
    print "runOnReReco and runOnSummer09 can't be true at the same time!"
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
'/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V3_test_production-v1/0000/962E4ECF-077C-DE11-871F-001F2907AF5C.root'
	)
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
if runOnData:
    if runOnReReco:
        process.GlobalTag.globaltag = cms.string('GR_R_35X_V7::All')
    else:
        process.GlobalTag.globaltag = cms.string('GR10_P_V5::All')
else:
    if runOnReReco:
        process.GlobalTag.globaltag = cms.string('START3X_V26::All')
    else:
        process.GlobalTag.globaltag = cms.string('MC_3XY_V26::All')

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

    process.p = cms.Path( process.primaryVertexFilter * process.scrapingFilter * process.technicals * process.incompleteECALFilter )


#add PAT to path
if runOnData:
    process.p += process.patDefaultSequence
else:
    process.p = cms.Path( process.patDefaultSequence )


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


#process.pTHat = cms.EDFilter("PtHatFilter",
#         pt_hat_lower_bound = cms.double(200.),
#         pt_hat_upper_bound = cms.double(500.)
#)

# cut on the factorization scale e.g. suitable for cuts on inv. mass of resonances in Pythia!
#process.FacScale = cms.EDFilter("FacScaleFilter",
#         fac_scale_lower_bound = cms.double(200.),
#         fac_scale_upper_bound = cms.double(500.)
#)


process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )

if runOnSummer09:
    #anti-kt 5 jets are called antikt5 in 3.1.X, but ak5 in later releases
    process.Skimmer.jets.AK5.MCLabel = 'antikt5GenJets'
    #add the high lumi trigger
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
    #and rename the low lumi trigger
    process.Skimmer.triggers.HLT8E29.process = 'HLT8E29'


if not runOnData:
    MUSiCProject.Skimming.Tools.addFlavourMatching( process, process.Skimmer, process.p )

process.p += process.Skimmer
