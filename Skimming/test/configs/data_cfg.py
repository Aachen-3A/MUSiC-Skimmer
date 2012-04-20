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
    process.Skimmer.triggers.HLT.HLTriggers = cms.vstring(
    'HLT_Mu15_v2',
    'HLT_Mu20_v1',
    'HLT_Mu24_v1',
    'HLT_Mu24_v2',
    'HLT_Mu30_v1',
    'HLT_Mu30_v2',
    'HLT_Mu30_v3',
    'HLT_Mu40_v1',
    'HLT_Mu40_v2',
    'HLT_Mu40_v3',
    'HLT_Mu40_v5',
    'HLT_Mu40_eta2p1_v1',
    'HLT_Mu40_eta2p1_v4',
    'HLT_Mu40_eta2p1_v5',

    'HLT_IsoMu12_v1',
    'HLT_IsoMu15_v5',
    'HLT_IsoMu17_v5',
    'HLT_IsoMu17_v6',
    'HLT_IsoMu17_v8',
    'HLT_IsoMu20_eta2p1_v1',
    'HLT_IsoMu24_v1',
    'HLT_IsoMu24_v2',
    'HLT_IsoMu24_v4',
    'HLT_IsoMu24_v5',
    'HLT_IsoMu24_v6',
    'HLT_IsoMu24_v7',
    'HLT_IsoMu24_v8',
    'HLT_IsoMu30_eta2p1_v3',
    'HLT_IsoMu30_eta2p1_v6',
    'HLT_IsoMu30_eta2p1_v7',

    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
    'HLT_Ele45_CaloIdVT_TrkIdT_v1',
    'HLT_Ele45_CaloIdVT_TrkIdT_v2',
    'HLT_Ele45_CaloIdVT_TrkIdT_v3',
    'HLT_Ele52_CaloIdVT_TrkIdT_v1',
    'HLT_Ele52_CaloIdVT_TrkIdT_v2',
    'HLT_Ele52_CaloIdVT_TrkIdT_v3',
    'HLT_Ele65_CaloIdVT_TrkIdT_v1',
    'HLT_Ele65_CaloIdVT_TrkIdT_v2',
    'HLT_Ele65_CaloIdVT_TrkIdT_v3',
    'HLT_Ele65_CaloIdVT_TrkIdT_v4',
    'HLT_Ele80_CaloIdVT_TrkIdT_v2',
    'HLT_Ele80_CaloIdVT_TrkIdT_v3',

    'HLT_Photon20_CaloIdVL_IsoL_v1',
    'HLT_Photon50_CaloIdVL_IsoL_v1',
    'HLT_Photon75_CaloIdVL_IsoL_v1',
    'HLT_Photon75_CaloIdVL_IsoL_v2',
    'HLT_Photon75_CaloIdVL_IsoL_v3',
    'HLT_Photon75_CaloIdVL_IsoL_v4',
    'HLT_Photon90_CaloIdVL_IsoL_v1',
    'HLT_Photon90_CaloIdVL_IsoL_v2',
    'HLT_Photon90_CaloIdVL_IsoL_v3',
    'HLT_Photon90_CaloIdVL_IsoL_v5',
    'HLT_Photon135_v1',
    'HLT_Photon135_v2',

    'HLT_Jet240_v1',
    'HLT_Jet300_v1',
    'HLT_Jet300_v2',
    'HLT_Jet300_v3',
    'HLT_Jet300_v4',
    'HLT_Jet300_v5',
    'HLT_Jet300_v6',
    'HLT_Jet370_v6',
    'HLT_Jet370_v7',
    'HLT_Jet370_v10',

    'HLT_MET100_v1',
    'HLT_MET120_v1',
    'HLT_MET120_v2',
    'HLT_MET120_v3',
    'HLT_MET120_v4',
    'HLT_MET120_v5',
    'HLT_MET120_v6',
    'HLT_MET120_v7',
    'HLT_MET200_v1',
    'HLT_MET200_v2',
    'HLT_MET200_v3',
    'HLT_MET200_v4',
    'HLT_MET200_v5',
    'HLT_MET200_v6',
    'HLT_MET200_v7',

    'HLT_Mu13_Mu8_v2',
    'HLT_Mu13_Mu8_v3',
    'HLT_Mu13_Mu8_v4',
    'HLT_Mu13_Mu8_v6',
    'HLT_Mu13_Mu8_v7',
    'HLT_Mu17_Mu8_v7',
    'HLT_Mu17_Mu8_v10',
    'HLT_Mu17_Mu8_v11',
    'HLT_Mu17_TkMu8_v3',
    'HLT_Mu17_TkMu8_v4',

    'HLT_DoubleMu6_v1',
    'HLT_DoubleMu7_v1',
    'HLT_DoubleMu7_v2',

    'HLT_IsoMu12_LooseIsoPFTau10_v1',
    'HLT_IsoMu12_LooseIsoPFTau10_v2',
    'HLT_IsoMu12_LooseIsoPFTau10_v4',
    'HLT_IsoMu15_LooseIsoPFTau15_v2',
    'HLT_IsoMu15_LooseIsoPFTau15_v4',
    'HLT_IsoMu15_LooseIsoPFTau15_v5',
    'HLT_IsoMu15_LooseIsoPFTau15_v6',
    'HLT_IsoMu15_LooseIsoPFTau15_v8',
    'HLT_IsoMu15_TightIsoPFTau20_v2',
    'HLT_IsoMu15_TightIsoPFTau20_v3',
    'HLT_IsoMu15_TightIsoPFTau20_v4',
    'HLT_IsoMu15_TightIsoPFTau20_v6',
    'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1',
    'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v5',
    'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v6',
    'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v1',
    'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v5',
    'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v6',
    'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v1',
    'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v5',
    'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v6',

    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v1',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v1',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v6',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v8',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v9',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v3',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1',
    'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1',
    'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v5',
    'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v6',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2',

    'HLT_IsoPFTau35_Trk20_MET45_v1',
    'HLT_IsoPFTau35_Trk20_MET45_v2',
    'HLT_IsoPFTau35_Trk20_MET45_v4',
    'HLT_IsoPFTau35_Trk20_MET45_v6',
    'HLT_IsoPFTau35_Trk20_MET60_v2',
    'HLT_IsoPFTau35_Trk20_MET60_v3',
    'HLT_IsoPFTau35_Trk20_MET60_v4',
    'HLT_IsoPFTau35_Trk20_MET60_v6',

    'HLT_MediumIsoPFTau35_Trk20_MET60_v1',
    'HLT_MediumIsoPFTau35_Trk20_MET60_v5',
    'HLT_MediumIsoPFTau35_Trk20_MET60_v6',
    )

print 'INFO: Using global tag:', process.GlobalTag.globaltag
