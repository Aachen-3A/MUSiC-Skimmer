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

import PxlSkimmer.Skimming.Tools as Tools

process = Tools.prepare( runOnGen, runOnData, eleEffAreaTarget, verbosityLvl )

# source
process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12_DR53X/T_t-channel_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v3/00000/54834283-CE0F-E211-AFB8-003048678FAE.root'
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )


if not runOnGen:
    if not runOnData:
        # This is used for QCD samples only, and the filters only write a flag in the
        # event if they fired or not. The actual selection must happen in the
        # classification, i.e. you have to set in the config file which flag you want to
        # consider.
        #
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


print 'INFO: Using global tag:', process.GlobalTag.globaltag
