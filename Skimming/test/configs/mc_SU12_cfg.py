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
        '/store/mc/Summer12/TTJets_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S8_START52_V9-v1/0000/8EF198B4-FBC1-E111-8882-003048FFD728.root'
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
