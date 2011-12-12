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
    process.Skimmer.triggers.HLT.HLTriggers = cms.vstring( 'HLT_Mu30_v3',
                                                           'HLT_Mu40_v1',
                                                           'HLT_Mu40_v2',
                                                           'HLT_Mu40_v3',
                                                           'HLT_Mu40_v5',
                                                           'HLT_Mu40_eta2p1_v1',
                                                           'HLT_Mu40_eta2p1_v4',
                                                           'HLT_Mu40_eta2p1_v5',
                                                           'HLT_IsoMu17_v8',
                                                           'HLT_IsoMu24_v4',
                                                           'HLT_IsoMu24_v5',
                                                           'HLT_IsoMu24_v6',
                                                           'HLT_IsoMu24_v7',
                                                           'HLT_IsoMu24_v8',
                                                           'HLT_IsoMu30_eta2p1_v3',
                                                           'HLT_IsoMu30_eta2p1_v6',
                                                           'HLT_IsoMu30_eta2p1_v7',

                                                           'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
                                                           'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                                                           'HLT_Ele52_CaloIdVT_TrkIdT_v1',
                                                           'HLT_Ele52_CaloIdVT_TrkIdT_v2',
                                                           'HLT_Ele52_CaloIdVT_TrkIdT_v3',
                                                           'HLT_Ele65_CaloIdVT_TrkIdT_v3',
                                                           'HLT_Ele65_CaloIdVT_TrkIdT_v4',
                                                           'HLT_Ele80_CaloIdVT_TrkIdT_v2',
                                                           'HLT_Ele80_CaloIdVT_TrkIdT_v3',

                                                           'HLT_Photon75_CaloIdVL_IsoL_v4',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v1',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v2',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v3',
                                                           'HLT_Photon90_CaloIdVL_IsoL_v5',
                                                           'HLT_Photon135_v1',
                                                           'HLT_Photon135_v2',

                                                           'HLT_Jet300_v2',
                                                           'HLT_Jet300_v3',
                                                           'HLT_Jet300_v4',
                                                           'HLT_Jet300_v5',
                                                           'HLT_Jet300_v6',
                                                           'HLT_Jet370_v10',

                                                           'HLT_MET200_v3',
                                                           'HLT_MET200_v4',
                                                           'HLT_MET200_v5',
                                                           'HLT_MET200_v6',
                                                           'HLT_MET200_v7'
                                                           )

print 'INFO: Using global tag:', process.GlobalTag.globaltag
