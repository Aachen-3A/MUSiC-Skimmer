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
    #'/store/data/Run2011A/MET/AOD/PromptReco-v4/000/168/423/7C1E218D-1AA6-E011-A278-003048F11114.root'
    '/store/data/Run2011A/MET/AOD/PromptReco-v4/000/166/512/80F7C542-ED91-E011-99B6-001D09F24259.root'
    )
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
if runOnData:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
else:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.GlobalTag.globaltag = 'GR_R_42_V19::All'
print "INFO: Using global tag:", process.GlobalTag.globaltag

process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")

import MUSiCProject.Skimming.Tools

if not runOnGen:
   MUSiCProject.Skimming.Tools.configurePAT( process, runOnData, runOnReReco, runOnSummer09 )
   process.metJESCorAK5CaloJet.inputUncorMetLabel = 'metNoHF'

   process.load( 'JetMETCorrections.Configuration.DefaultJEC_cff' )
   process.load( 'RecoJets.Configuration.RecoPFJets_cff' )
   process.kt6PFJets.doRhoFastjet = True                          # Turn-on the FastJet density calculation
   process.ak5PFJets.doAreaFastjet = True                         # Turn-on the FastJet jet area calculation for ak5PFJets

   process.p = cms.Path( process.kt6PFJets * process.ak5PFJets )

   # for "PFnoPU" #
   from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

   # Create good primary vertices to be used for PF association
   process.goodOfflinePrimaryVertices = cms.EDFilter(
      'PrimaryVertexObjectFilter',
      filterParams = pvSelector.clone( minNdof = cms.double( 4.0 ), maxZ = cms.double( 24.0 ) ),
      src = cms.InputTag( 'offlinePrimaryVertices' )
      )
   ################

   from PhysicsTools.PatAlgos.tools import pfTools
   postfix = 'PFlow'
   pfTools.usePF2PAT( process, runPF2PAT = True, jetAlgo = 'AK5', runOnMC = not runOnData, postfix = postfix )
   #process.patJetCorrFactorsPFlow.levels = cms.vstring( 'L1Offset', 'L2Relative', 'L3Absolute' )
   if runOnData: process.patJetCorrFactorsPFlow.levels = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' )

   # for "PFnoPU" #
   process.pfPileUpPFlow.Enable = True
   process.pfPileUpPFlow.Vertices = 'goodOfflinePrimaryVertices'
   process.pfJetsPFlow.doAreaFastjet = True
   process.pfJetsPFlow.doRhoFastjet = False
   process.patJetCorrFactorsPFlow.rho = cms.InputTag( 'kt6PFJetsPFlow', 'rho' )
   process.pfPileUpPFlow.checkClosestZVertex = cms.bool( False )

   # Compute the mean pt per unit area ("rho") using KT6 Jets with the active areas method.
   from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
   process.kt6PFJetsPFlow = kt4PFJets.clone(
      rParam = cms.double(0.6),
      src = cms.InputTag( 'pfNoElectron' + postfix ),
      doAreaFastjet = cms.bool( True ),
      doRhoFastjet = cms.bool( True )
      )
   ################

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

    process.p += process.scrapingFilter * process.patDefaultSequence

else:
   if runOnSummer09:
      process.load("RecoJets.Configuration.GenJetParticles_cff")
      process.load("RecoJets.JetProducers.ak5GenJets_cfi")
      process.p += process.genParticlesForJets * process.ak5GenJets * process.patDefaultSequence
   elif runOnGen:
      process.p = cms.Path( process.patJetPartons )
   else:
      process.p += process.patDefaultSequence

if not runOnGen:
   # for "PFnoPU" #
   getattr( process, 'patPF2PATSequence' + postfix ).replace( getattr( process, 'pfNoElectron' + postfix ), getattr( process, 'pfNoElectron' + postfix ) * process.kt6PFJetsPFlow )
   process.patseq = cms.Sequence(
      process.goodOfflinePrimaryVertices*
      getattr( process, 'patPF2PATSequence' + postfix )
      )
   ################
   process.p += process.patseq

   # METnoPU (stolen from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/lucieg/METsWithPU/METsAnalyzer/python/pfMetNoPileUp_cff.py)
   process.pfMetNoPileUp       = getattr( process, 'pfMET' + postfix ).clone()
   process.pfMetNoPileUp.alias = 'pfMetNoPileUp'
   process.pfMetNoPileUp.src   = cms.InputTag( 'pfNoPileUp' + postfix )

   process.p += process.pfMetNoPileUp

   patMETsPFlowNoPU = 'patMETs' + postfix + 'NoPU'
   setattr( process, patMETsPFlowNoPU, getattr( process, 'patMETs' + postfix ).clone() )
   getattr( process, patMETsPFlowNoPU ).metSource = cms.InputTag( 'pfMetNoPileUp' )

   process.p += getattr( process, patMETsPFlowNoPU )


   #store the result of the HCAL noise info
   process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
   process.p += process.HBHENoiseFilterResultProducer

process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )

if runOnData:
   process.Skimmer.triggers.HLT.HLTriggers = cms.vstring( 'HLT_Mu30_v3',
                                                          'HLT_Mu40_v1',
                                                          'HLT_Mu40_v2',
                                                          'HLT_Mu40_v3',
                                                          'HLT_IsoMu17_v8',
                                                          'HLT_IsoMu24_v4',
                                                          'HLT_IsoMu24_v5',
                                                          'HLT_IsoMu24_v6',
                                                          'HLT_IsoMu24_v7',

                                                          'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
                                                          'HLT_Ele52_CaloIdVT_TrkIdT_v1',
                                                          'HLT_Ele52_CaloIdVT_TrkIdT_v2',
                                                          'HLT_Ele52_CaloIdVT_TrkIdT_v3',

                                                          'HLT_Photon75_CaloIdVL_IsoL_v4',
                                                          'HLT_Photon90_CaloIdVL_IsoL_v1',
                                                          'HLT_Photon90_CaloIdVL_IsoL_v2',
                                                          'HLT_Photon90_CaloIdVL_IsoL_v3',

                                                          'HLT_Jet300_v2',
                                                          'HLT_Jet300_v3',
                                                          'HLT_Jet300_v4',
                                                          'HLT_Jet300_v5',

                                                          'HLT_MET200_v3',
                                                          'HLT_MET200_v4',
                                                          'HLT_MET200_v5',
                                                          'HLT_MET200_v6'
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
