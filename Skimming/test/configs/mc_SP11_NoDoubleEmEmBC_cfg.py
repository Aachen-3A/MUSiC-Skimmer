runOnData = False
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# source
process.source = cms.Source("PoolSource",
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
'/store/mc/Spring11/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S1_START311_V1G1-v1/0015/B0298296-594F-E011-AA02-1CC1DE1CEFC8.root'
	)
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
if runOnData:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
else:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

import MUSiCProject.Skimming.Tools
MUSiCProject.Skimming.Tools.configurePAT( process, runOnData, runOnReReco, runOnSummer09 )
process.metJESCorAK5CaloJet.inputUncorMetLabel = 'metNoHF'

from PhysicsTools.PatAlgos.tools import pfTools
pfTools.usePF2PAT( process, runPF2PAT=True, jetAlgo='AK5', runOnMC= not runOnData, postfix="PFlow" )
process.patJetCorrFactorsPFlow.levels = cms.vstring( L1Offset, L2Relative, L3Absolute )

if runOnData:
    import PhysicsTools.PatAlgos.tools.coreTools
    PhysicsTools.PatAlgos.tools.coreTools.removeMCMatching( process, ['All'] )


#filter on right BX in case of data
if runOnData:
    process.scrapingFilter = cms.EDFilter( "FilterOutScraping",
                                           applyfilter = cms.untracked.bool( True ),
                                           debugOn = cms.untracked.bool( False ),
                                           numtrack = cms.untracked.uint32( 10 ),
                                           thresh = cms.untracked.double( 0.25 )
                                           )


    process.p = cms.Path( process.scrapingFilter * process.patDefaultSequence )

else:
    #this filter selects events containing two potential electron candidates
    process.doubleEMenrichingfilter = cms.EDFilter( "doubleEMEnrichingFilter",
                                                    filterAlgoPSet = cms.PSet( requireTrackMatch = cms.bool( False ),
                                                                               caloIsoMax = cms.double( 3.0 ),
                                                                               isoGenParConeSize = cms.double( 0.1 ),
                                                                               tkIsoMax = cms.double( 3.0 ),
                                                                               isoConeSize = cms.double( 0.2 ),
                                                                               isoGenParETMin = cms.double( 4.0 ),
                                                                               hOverEMax = cms.double( 0.5 ),
                                                                               clusterThreshold = cms.double( 4.0 ),
                                                                               seedThreshold = cms.double( 3.5 ),
                                                                               eTThreshold = cms.double( 3.0 ),
                                                                               genParSource = cms.InputTag( "genParticles" )
                                                                               )
                                                    )

    #this filter selects events containing one potential electron candidate
    process.emenrichingfilter = cms.EDFilter( "EMEnrichingFilter",
                                              filterAlgoPSet = cms.PSet( requireTrackMatch = cms.bool( False ),
                                                                         caloIsoMax = cms.double( 10.0 ),
                                                                         isoGenParConeSize = cms.double( 0.1 ),
                                                                         tkIsoMax = cms.double( 5.0 ),
                                                                         hOverEMax = cms.double( 0.5 ),
                                                                         isoGenParETMin = cms.double( 20.0 ),
                                                                         genParSource = cms.InputTag( "genParticles" ),
                                                                         isoConeSize = cms.double( 0.2 ),
                                                                         clusterThreshold = cms.double( 20.0 )
                                                                         )
                                              )

    #this filter selects events containing electrons that come from b or c hadrons
    process.bctoefilter = cms.EDFilter( "BCToEFilter",
                                        filterAlgoPSet = cms.PSet( genParSource = cms.InputTag( "genParticles" ),
                                                                   eTThreshold = cms.double(10)
                                                                   )
                                        )

    #but we don't want either of those events in this sample
    process.p = cms.Path( ~process.doubleEMenrichingfilter + ~process.emenrichingfilter + ~process.bctoefilter )
    if runOnSummer09:
        process.load("RecoJets.Configuration.GenJetParticles_cff")
        process.load("RecoJets.JetProducers.ak5GenJets_cfi")
        process.p += process.genParticlesForJets * process.ak5GenJets
    #add PAT
    process.p += process.patDefaultSequence


process.p += getattr(process,"patPF2PATSequencePFlow")
#store the result of the HCAL noise info
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.p += process.HBHENoiseFilterResultProducer


process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )

if not runOnData:
    MUSiCProject.Skimming.Tools.addFlavourMatching( process, process.Skimmer, process.p )

process.p += process.Skimmer
