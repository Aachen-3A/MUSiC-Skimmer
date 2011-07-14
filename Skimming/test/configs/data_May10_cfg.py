runOnData = True
#run on Summer09 (ReReco'd or not)
runOnSummer09 = False
#run on ReReco'ed data or Summer09 MC
runOnReReco = False
#run on GEN sample
runOnGen = False

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
    '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/7AF46394-8A7C-E011-9B29-00237DDC5B9E.root'
    )
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
if runOnData:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
else:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

print "INFO: Using global tag:", process.GlobalTag.globaltag
process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

import MUSiCProject.Skimming.Tools

if not runOnGen:
   MUSiCProject.Skimming.Tools.configurePAT( process, runOnData, runOnReReco, runOnSummer09 )
   process.metJESCorAK5CaloJet.inputUncorMetLabel = 'metNoHF'

   from PhysicsTools.PatAlgos.tools import pfTools
   pfTools.usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC= not runOnData, postfix="PFlow")
   process.patJetCorrFactorsPFlow.levels = cms.vstring( 'L1Offset', 'L2Relative', 'L3Absolute' )

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
   if runOnSummer09:
      process.load("RecoJets.Configuration.GenJetParticles_cff")
      process.load("RecoJets.JetProducers.ak5GenJets_cfi")
      process.p = cms.Path( process.genParticlesForJets * process.ak5GenJets * process.patDefaultSequence )
   elif runOnGen:
      process.p = cms.Path( process.patJetPartons )
   else:
      process.p = cms.Path( process.patDefaultSequence )

if not runOnGen:
   process.p += getattr(process,"patPF2PATSequencePFlow")
   #store the result of the HCAL noise info
   process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
   process.p += process.HBHENoiseFilterResultProducer

process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )

if runOnData:
   process.Skimmer.triggers.HLT.HLTriggers = cms.vstring( 'HLT_Mu15_v2',
                                                          'HLT_Mu20_v1',
                                                          'HLT_Mu24_v1',
                                                          'HLT_Mu24_v2',
                                                          'HLT_IsoMu12_v1',
                                                          'HLT_IsoMu15_v5',
                                                          'HLT_IsoMu17_v5',
                                                          'HLT_IsoMu17_v6',

                                                          'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                                                          'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
                                                          'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
                                                          'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                                                          'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',

                                                          'HLT_Photon20_CaloIdVL_IsoL_v1',
                                                          'HLT_Photon50_CaloIdVL_IsoL_v1',
                                                          'HLT_Photon75_CaloIdVL_IsoL_v1',
                                                          'HLT_Photon75_CaloIdVL_IsoL_v2',
                                                          'HLT_Photon75_CaloIdVL_IsoL_v3',
                                                          'HLT_Photon75_CaloIdVL_v1',
                                                          'HLT_Photon75_CaloIdVL_v2',
                                                          'HLT_Photon75_CaloIdVL_v3',

                                                          'HLT_Jet240_v1',
                                                          'HLT_Jet300_v1',
                                                          'HLT_Jet370_NoJetID_v1',
                                                          'HLT_Jet370_NoJetID_v2',
                                                          'HLT_Jet370_v1',
                                                          'HLT_Jet370_v2',

                                                          'HLT_MET200_v1',
                                                          'HLT_MET200_v2'
                                                         )
# HLTs for Summer11 MCs (using HLT config: /cdaq/physics/Run2011/5e32/v6.2/HLT/V1)
else:
   process.Skimmer.triggers.HLT.HLTriggers = cms.vstring( 'HLT_Mu15_v2',
                                                          'HLT_Mu20_v1',
                                                          'HLT_Mu24_v1',
                                                          'HLT_Mu30_v1',
                                                          'HLT_IsoMu12_v1',
                                                          'HLT_IsoMu15_v5',
                                                          'HLT_IsoMu17_v5',
                                                          'HLT_IsoMu24_v1',
                                                          'HLT_IsoMu30_v1',

                                                          'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
                                                          'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                                                          'HLT_Ele45_CaloIdVT_TrkIdT_v2',

                                                          'HLT_Photon50_CaloIdVL_IsoL_v1',
                                                          'HLT_Photon75_CaloIdVL_v2',
                                                          'HLT_Photon75_CaloIdVL_IsoL_v2',

                                                          'HLT_Jet240_v1',
                                                          'HLT_Jet370_v1',
                                                          'HLT_Jet370_NoJetID_v1',

                                                          'HLT_MET200_v1'
                                                          )

if not runOnData:
    MUSiCProject.Skimming.Tools.addFlavourMatching( process, process.Skimmer, process.p, runOnGen )

process.p += process.Skimmer
