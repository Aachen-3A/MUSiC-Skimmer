runOnData = True
#run on GEN sample
runOnGen = False
# Choose the type of effective area correction you want to use.
# Possible values:
#     NoCorr
#     Data2011
#     Data2012
#     Summer11MC
#     Fall11MC
eleEffAreaTarget = cms.untracked.string( 'Data2011' )

if runOnGen and runOnData :
    print "runOnData and runOnGen can't be true at the same time!"
    import sys
    sys.exit(1)

import MUSiCProject.Skimming.Tools as Tools

process = Tools.prepare( runOnGen, runOnData, eleEffAreaTarget )

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
    process.Skimmer.triggers.HLT.HLTriggers.extend( [

    # Single Mu:
    #
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
    'HLT_IsoMu34_eta2p1_v1',

    # Single Ele:
    #
    'HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
    'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
    'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3',
    'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',

    'HLT_Ele52_CaloIdVT_TrkIdT_v1',
    'HLT_Ele52_CaloIdVT_TrkIdT_v2',
    'HLT_Ele52_CaloIdVT_TrkIdT_v3',
    'HLT_Ele80_CaloIdVT_TrkIdT_v2',
    'HLT_Ele80_CaloIdVT_TrkIdT_v3',

    # Tau
    #
    # Single Tau
    'HLT_MediumIsoPFTau35_Trk20_v1',
    'HLT_MediumIsoPFTau35_Trk20_v5',
    'HLT_MediumIsoPFTau35_Trk20_v6',

    # Double Tau
    'HLT_IsoPFTau40_IsoPFTau30_Trk5_eta2p1_v2',

    'HLT_DoubleIsoPFTau20_Trk5_v1',
    'HLT_DoubleIsoPFTau20_Trk5_v2',
    'HLT_DoubleIsoPFTau20_Trk5_v4',
    'HLT_DoubleIsoPFTau25_Trk5_eta2p1_v2',
    'HLT_DoubleIsoPFTau35_Trk5_eta2p1_v2',
    'HLT_DoubleIsoPFTau35_Trk5_eta2p1_v3',
    'HLT_DoubleIsoPFTau35_Trk5_eta2p1_v4',
    'HLT_DoubleIsoPFTau40_Trk5_eta2p1_v2',
    'HLT_DoubleIsoPFTau40_Trk5_eta2p1_v3',
    'HLT_DoubleIsoPFTau40_Trk5_eta2p1_v4',
    'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v2',
    'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v3',
    'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v7',
    'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v8',
    'HLT_DoubleIsoPFTau55_Trk5_eta2p1_v4',
    'HLT_DoubleIsoPFTau55_Trk5_eta2p1_v5',

    # Tau + MET
    'HLT_IsoPFTau35_Trk20_MET45_v1',
    'HLT_IsoPFTau35_Trk20_MET45_v2',
    'HLT_IsoPFTau35_Trk20_MET45_v4',
    'HLT_IsoPFTau35_Trk20_MET45_v6',
    'HLT_IsoPFTau35_Trk20_MET60_v2',
    'HLT_IsoPFTau35_Trk20_MET60_v3',
    'HLT_IsoPFTau35_Trk20_MET60_v4',
    'HLT_IsoPFTau35_Trk20_MET60_v6',
    'HLT_IsoPFTau35_Trk20_MET70_v2',
    'HLT_IsoPFTau45_Trk20_MET60_v2',
    'HLT_IsoPFTau45_Trk20_MET60_v3',
    'HLT_IsoPFTau45_Trk20_MET60_v4',

    'HLT_MediumIsoPFTau35_Trk20_MET60_v1',
    'HLT_MediumIsoPFTau35_Trk20_MET60_v5',
    'HLT_MediumIsoPFTau35_Trk20_MET60_v6',
    'HLT_MediumIsoPFTau35_Trk20_MET70_v1',

    # Mu + Tau
    'HLT_Mu15_LooseIsoPFTau20_v1',
    'HLT_Mu15_LooseIsoPFTau20_v2',

    'HLT_IsoMu12_LooseIsoPFTau10_v1',
    'HLT_IsoMu12_LooseIsoPFTau10_v2',
    'HLT_IsoMu12_LooseIsoPFTau10_v4',
    'HLT_IsoMu15_LooseIsoPFTau15_v2',
    'HLT_IsoMu15_LooseIsoPFTau15_v4',
    'HLT_IsoMu15_LooseIsoPFTau15_v5',
    'HLT_IsoMu15_LooseIsoPFTau15_v6',
    'HLT_IsoMu15_LooseIsoPFTau15_v8',
    'HLT_IsoMu15_LooseIsoPFTau20_v2',
    'HLT_IsoMu15_LooseIsoPFTau20_v3',
    'HLT_IsoMu15_LooseIsoPFTau20_v4',
    'HLT_IsoMu15_LooseIsoPFTau20_v6',
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

    # Ele + Tau
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v1',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v1',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v6',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v8',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v9',
    'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v3',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1',
    'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2',
    'HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v1',
    'HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v5',
    'HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v6',
    'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1',
    'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v5',
    'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v6',

    # Single Photon:
    #
    'HLT_Photon20_CaloIdVL_IsoL_v1',
    'HLT_Photon30_CaloIdVL_IsoL_v8',
    'HLT_Photon30_CaloIdVL_v7',
    'HLT_Photon50_CaloIdVL_IsoL_v1',
    'HLT_Photon75_CaloIdVL_IsoL_v1',
    'HLT_Photon75_CaloIdVL_IsoL_v2',
    'HLT_Photon75_CaloIdVL_IsoL_v3',
    'HLT_Photon75_CaloIdVL_IsoL_v4',
    'HLT_Photon75_CaloIdVL_v1',
    'HLT_Photon75_CaloIdVL_v2',
    'HLT_Photon75_CaloIdVL_v3',
    'HLT_Photon90_CaloIdVL_IsoL_v1',
    'HLT_Photon90_CaloIdVL_IsoL_v2',
    'HLT_Photon90_CaloIdVL_IsoL_v3',
    'HLT_Photon90_CaloIdVL_IsoL_v5',
    'HLT_Photon135_v1',
    'HLT_Photon135_v2',

    # Single Jet:
    #
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

    # MET:
    #
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

    # Double Mu:
    #
    'HLT_DoubleMu6_v1',
    'HLT_DoubleMu7_v1',
    'HLT_DoubleMu7_v2',

    'HLT_Mu13_Mu8_v2',
    'HLT_Mu13_Mu8_v3',
    'HLT_Mu13_Mu8_v4',
    'HLT_Mu13_Mu8_v6',
    'HLT_Mu13_Mu8_v7',
    'HLT_Mu17_Mu8_v10',
    'HLT_Mu17_Mu8_v11',
    'HLT_Mu17_Mu8_v2',
    'HLT_Mu17_Mu8_v3',
    'HLT_Mu17_Mu8_v4',
    'HLT_Mu17_Mu8_v6',
    'HLT_Mu17_Mu8_v7',
    'HLT_Mu17_TkMu8_v3',
    'HLT_Mu17_TkMu8_v4',

    # Double Ele:
    #
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1',
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2',
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3',
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4',
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5',
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8',
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5',

    # Double Photon:
    #
    'HLT_Photon26_CaloIdL_IsoVL_Photon18_v1',
    'HLT_Photon26_CaloIdL_IsoVL_Photon18_v2',
    'HLT_Photon26_CaloIdL_IsoVL_Photon18_v3',
    'HLT_Photon26_CaloIdL_IsoVL_Photon18_v4',
    'HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v1',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v1',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v2',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v3',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v4',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v5',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v6',
    'HLT_Photon26_IsoVL_Photon18_IsoVL_v7',
    'HLT_Photon26_IsoVL_Photon18_v1',
    'HLT_Photon26_IsoVL_Photon18_v2',
    'HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1',
    'HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2',
    'HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v1',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v2',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v3',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v4',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v6',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v7',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v1',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v2',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_v3',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_v1',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_v2',
    'HLT_Photon36_CaloIdL_IsoVL_Photon22_v3',
    'HLT_Photon36_CaloIdL_Photon22_CaloIdL_v1',
    'HLT_Photon36_CaloIdL_Photon22_CaloIdL_v2',
    'HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3',
    'HLT_Photon36_CaloIdL_Photon22_CaloIdL_v4',
    'HLT_Photon36_IsoVL_Photon22_v1',
    'HLT_Photon36_IsoVL_Photon22_v2',
    'HLT_Photon36_IsoVL_Photon22_v3',

    # Single Mu + Single Ele:
    #
    'HLT_Mu8_Ele17_CaloIdL_v1',
    'HLT_Mu8_Ele17_CaloIdL_v2',
    'HLT_Mu8_Ele17_CaloIdL_v3',
    'HLT_Mu8_Ele17_CaloIdL_v4',
    'HLT_Mu8_Ele17_CaloIdL_v5',
    'HLT_Mu8_Ele17_CaloIdL_v6',
    'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1',
    'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3',
    'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4',
    'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7',
    'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8',
    'HLT_Mu10_Ele10_CaloIdL_v2',
    'HLT_Mu10_Ele10_CaloIdL_v3',
    'HLT_Mu10_Ele10_CaloIdL_v4',
    'HLT_Mu17_Ele8_CaloIdL_v1',
    'HLT_Mu17_Ele8_CaloIdL_v2',
    'HLT_Mu17_Ele8_CaloIdL_v3',
    'HLT_Mu17_Ele8_CaloIdL_v4',
    'HLT_Mu17_Ele8_CaloIdL_v5',
    'HLT_Mu17_Ele8_CaloIdL_v6',
    'HLT_Mu17_Ele8_CaloIdL_v8',
    'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1',
    'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3',
    'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4',
    'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7',
    'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8',
    ] )


print 'INFO: Using global tag:', process.GlobalTag.globaltag
