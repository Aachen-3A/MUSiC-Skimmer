runOnData = False
#run on GEN sample
runOnGen = False

import FWCore.ParameterSet.Config as cms

# Choose the type of effective area correction you want to use.
# Possible values:
#     NoCorr
#     Data2011
#     Data2012
#     Summer11MC
#     Fall11MC
eleEffAreaTarget = cms.untracked.string( 'Fall11MC' )

# Verbosity: 0 = normal messaging, 1 = human readable, 2 = insane, 3 = INFO from all modules
verbosityLvl = 0

if runOnGen and runOnData :
    print "runOnData and runOnGen can't be true at the same time!"
    import sys
    sys.exit( 1 )

import MUSiCProject.Skimming.Tools as Tools

process = Tools.prepare( runOnGen, runOnData, eleEffAreaTarget, verbosityLvl )

# source
process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    fileNames = cms.untracked.vstring(
        #'/store/mc/Fall11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/16DC5C61-BF07-E111-88BB-0030487D5EB3.root'
        '/store/mc/Fall11/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/AODSIM/PU_S6-START44_V5-v1/0001/4031B29A-0306-E111-9E34-0030487D8661.root'
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )

# This is used for QCD samples only, and the filters only write a flag in the
# event if they fired or not. The actual selection must happen in the
# classification, i.e. you have to set in the config file which flag you want to
# consider.
#
if not runOnGen:
    if not runOnData:
        # remove events with electrons that come from b or c hadrons
        Tools.addBCtoEFilter( process )

        # remove events containing b quarks
        Tools.addBFilter( process )

        # remove events with at least one potential electron candidate
        Tools.addEMFilter( process )

        # remove events with a muon of more than 5 GeV
        Tools.addMuGenFilter( process, pt = 5 )

        # remove events with a muon of more than 10 GeV
        Tools.addMuGenFilter( process, pt = 10 )

        # remove events with a muon of more than 15 GeV
        Tools.addMuGenFilter( process, pt = 15 )

        process.Skimmer.triggers.HLT.HLTriggers = cms.vstring(
            # MC HLTs for Fall11 samples:
            # ---------------------------
            # There are different production revisions that were processed
            # with different CMSSW releases. Each of these used a different
            # HLT config, so different triggers are available.
            #
            # Always using the HLT config found in 'CMSSW/HLTrigger/Configuration/python/HLT_GRun_cff.py'
            # for the respective release.
            #
            # Fall11_R1 = CMSSW_4_2_9_HLT1_patch1 (/online/collisions/2011/3e33/v1.1/HLT/V10)
            #
            # Single Mu:
            #
            'HLT_Mu40_eta2p1_v1',
            'HLT_IsoMu30_eta2p1_v3',

            # Single Ele:
            #
            'HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',

            # Tau:
            #
            # Single Tau:
            'HLT_MediumIsoPFTau35_Trk20_v1',

            # Tau + MET:
            'HLT_MediumIsoPFTau35_Trk20_MET70_v1',

            # Double Tau:
            'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v3',

            # Mu + Tau:
            'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1',
            'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v1',
            'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v1',

            # Ele + Tau:
            'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1',
            'HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v1',

            # Single Photon:
            #
            'HLT_Photon30_CaloIdVL_v7',

            # Single Jet:
            #
            'HLT_Jet300_v5',

            # MET:
            #
            'HLT_MET200_v7',

            # Double Mu:
            #
            'HLT_Mu13_Mu8_v7',

            # Double Ele:
            #
            'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4',

            # Double Photon:
            #
            'HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v1',

            # Single Mu + Single Ele:
            #
            'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4',


            # Fall11_R2 = CMSSW_4_2_8_patch7 # (/online/collisions/2011/5e32/v6.2/HLT/V4)
            #
            # Single Mu:
            #
            'HLT_Mu30_v1',
            'HLT_IsoMu12_v1',

            # Single Ele:
            #
            'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',

            # Tau:
            #
            # Double Tau:
            'HLT_DoubleIsoPFTau20_Trk5_v2',

            # Tau + MET:
            'HLT_IsoPFTau35_Trk20_MET45_v2',

            # Mu + Tau:
            'HLT_IsoMu12_LooseIsoPFTau10_v2',

            # Ele + Tau:
            'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2',

            # Single Photon:
            #
            'HLT_Photon50_CaloIdVL_IsoL_v1',

            # Single Jet:
            #
            'HLT_Jet240_v1',

            # MET:
            #
            'HLT_MET200_v1',

            # Double Mu:
            #
            'HLT_DoubleMu6_v1',

            # Double Ele:
            #
            'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2',

            # Double Photon:
            #
            'HLT_Photon26_IsoVL_Photon18_v2',

            # Single Mu + Single Ele:
            #
            'HLT_Mu10_Ele10_CaloIdL_v3',
            'HLT_Mu17_Ele8_CaloIdL_v2',
            'HLT_Mu8_Ele17_CaloIdL_v2',


            # Fall11_Chamonix = CMSSW_4_4_2_patch9 (/online/collisions/2011/5e33/v3.0/HLT/V6)
            # Fall11_R4 = CMSSW_4_4_2_patch8 (/online/collisions/2011/5e33/v3.0/HLT/V6)
            # Fall11_TSG = CMSSW_4_4_2_patch6 (/online/collisions/2011/5e33/v3.0/HLT/V4)
            #
            # Single Mu:
            #
            'HLT_Mu40_eta2p1_v5',
            'HLT_IsoMu30_eta2p1_v7',

            # Single Ele:
            #
            'HLT_Ele80_CaloIdVT_TrkIdT_v3',

            # Tau:
            #
            # Single Tau:
            'HLT_MediumIsoPFTau35_Trk20_v6',

            # Double Tau:
            'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v8',

            # Tau + MET:
            'HLT_MediumIsoPFTau35_Trk20_MET60_v6',

            # Mu + Tau:
            'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v6',
            'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v6',
            'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v6',

            # Ele + Tau:
            'HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v6',

            # Single Photon:
            #
            'HLT_Photon135_v2',

            # Single Jet:
            #
            'HLT_Jet370_v10',

            # MET:
            #
            'HLT_MET200_v7',

            # Double Mu:
            #
            'HLT_Mu17_Mu8_v11',
            'HLT_Mu17_TkMu8_v4',

            # Double Ele:
            #
            'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10',

            # Double Photon:
            #
            'HLT_Photon36_CaloIdL_IsoVL_Photon22_CaloIdL_IsoVL_v7',

            # Single Mu + Single Ele:
            #
            'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8',
            'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8',


            # Fall11_R3 = CMSSW_4_4_0_patch3 (/dev/CMSSW_4_2_0/GRun/V226)
            #
            # Single Mu:
            #
            'HLT_IsoMu30_eta2p1_v3',
            'HLT_Mu40_eta2p1_v1',

            # Single Ele:
            #
            'HLT_Ele25_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v5',

            # Tau:
            #
            # Single Tau:
            'HLT_MediumIsoPFTau35_Trk20_v1',

            # Double Tau:
            'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v3',

            # Tau + MET:
            'HLT_MediumIsoPFTau35_Trk20_MET70_v1',

            # Mu + Tau:
            'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1',
            'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v1',
            'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v1',

            # Ele + Tau:
            'HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v1',

            # Single Photon:
            #
            'HLT_Photon90_CaloIdVL_IsoL_v5',

            # Single Jet:
            #
            'HLT_Jet300_v5',

            # MET:
            #
            'HLT_MET200_v7',

            # Double Mu:
            #
            'HLT_Mu13_Mu8_v7',

            # Double Ele:
            #
            'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8',

            # Double Photon:
            #
            'HLT_Photon26_CaloIdXL_IsoXL_Photon18_CaloIdXL_IsoXL_v1',

            # Single Mu + Single Ele:
            #
            'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4',
            'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4',
            )

print 'INFO: Using global tag:', process.GlobalTag.globaltag
