runOnData = True
#run on GEN sample
runOnGen = False

if runOnGen and runOnData :
    print "runOnData and runOnGen can't be true at the same time!"
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
    #'/store/data/Run2011A/MET/AOD/PromptReco-v4/000/166/512/80F7C542-ED91-E011-99B6-001D09F24259.root'
    '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/D4C112A8-697C-E011-9625-0018F3D0960A.root'
    #'/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S4_START42_V11-v1/0000/C45712AA-A6A8-E011-8679-0024E86E8D4C.root'
    )
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
if runOnData:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
else:
    process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

process.load( 'Configuration.StandardSequences.GeometryPilot2_cff' )
process.load( 'Configuration.StandardSequences.MagneticField_38T_cff' )

process.content = cms.EDAnalyzer("EventContentAnalyzer")

import MUSiCProject.Skimming.Tools

# create an empty path because modules will be added by calling the functions below
process.p = cms.Path()

if not runOnGen:
    # Keep the following functions in the right order as they will add modules to the path!
    if runOnData: MUSiCProject.Skimming.Tools.addScrapingFilter( process )
    MUSiCProject.Skimming.Tools.configureJEC( process, runOnData )
    MUSiCProject.Skimming.Tools.configurePAT( process, runOnData )
    process.metJESCorAK5CaloJet.inputUncorMetLabel = 'metNoHF'

    postfix = 'PFlow'
    MUSiCProject.Skimming.Tools.configurePF( process, runOnData, postfix )
    MUSiCProject.Skimming.Tools.configurePFnoPU( process, postfix )

    import PhysicsTools.PatAlgos.tools.coreTools
    PhysicsTools.PatAlgos.tools.coreTools.removeMCMatching( process, ['All'] )

    # store the result of the HCAL noise info
    process.load( 'CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi' )
    process.p += process.HBHENoiseFilterResultProducer

if not runOnData:
   process.p += process.patJetPartons

#rerun PFTau reco
process.load( "RecoTauTag.Configuration.RecoPFTauTag_cff" )
process.p += process.PFTau

process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )

if runOnData:
    process.Skimmer.triggers.HLT.HLTriggers = cms.vstring( 'HLT_Mu30_v3',
                                                           'HLT_Mu40_v1',
                                                           'HLT_Mu40_v2',
                                                           'HLT_Mu40_v3',
                                                           'HLT_Mu40_v5',
                                                           'HLT_IsoMu17_v8',
                                                           'HLT_IsoMu24_v4',
                                                           'HLT_IsoMu24_v5',
                                                           'HLT_IsoMu24_v6',
                                                           'HLT_IsoMu24_v7',
                                                           'HLT_IsoMu24_v8',

                                                           'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
                                                           'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                                                           'HLT_Ele52_CaloIdVT_TrkIdT_v1',
                                                           'HLT_Ele52_CaloIdVT_TrkIdT_v2',
                                                           'HLT_Ele52_CaloIdVT_TrkIdT_v3',
                                                           'HLT_Ele65_CaloIdVT_TrkIdT_v3',

                                                           'HLT_Photon75_CaloIdVL_IsoL_v4',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v1',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v2',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v3',
                                                           'HLT_Photon135_v1',

                                                           'HLT_Jet300_v2',
                                                           'HLT_Jet300_v3',
                                                           'HLT_Jet300_v4',
                                                           'HLT_Jet300_v5',

                                                           'HLT_MET200_v3',
                                                           'HLT_MET200_v4',
                                                           'HLT_MET200_v5',
                                                           'HLT_MET200_v6',
                                                           'HLT_MET200_v7'
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

                                                           'HLT_MET200_v1',

                                                            #tau trigger
                                                           'HLT_IsoPFTau35_Trk20_MET45_v2',
                                                           'HLT_DoubleIsoPFTau20_Trk5_v2',
                                                           'HLT_Mu15_LooseIsoPFTau20_v2',
                                                           'HLT_IsoMu12_LooseIsoPFTau10_v2',
                                                           'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2'
                                                           )

if not runOnData:
    MUSiCProject.Skimming.Tools.addFlavourMatching( process, process.Skimmer, process.p, runOnGen )

process.p += process.Skimmer

print 'INFO: Using global tag:', process.GlobalTag.globaltag
