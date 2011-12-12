runOnData = True
#run on GEN sample
runOnGen = False

if runOnGen and runOnData :
    print "runOnData and runOnGen can't be true at the same time!"
    import sys
    sys.exit(1)

import MUSiCProject.Skimming.Tools as Tools

process = Tools.prepare( runOnGen, runOnData )

import FWCore.ParameterSet.Config as cms

# source
process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    fileNames = cms.untracked.vstring(
        '/store/data/Run2011A/METBTag/AOD/May10ReReco-v1/0000/D4C112A8-697C-E011-9625-0018F3D0960A.root'
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )

# Take care of HLTs in data.
# HLTs for Summer11 MCs configured in MUSiCSkimmer_cfi.py
# (using HLT config: /cdaq/physics/Run2011/5e32/v6.2/HLT/V1)
#
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

print 'INFO: Using global tag:', process.GlobalTag.globaltag
