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

import MUSiCProject.Skimming.Tools_miniAOD as Tools_miniAOD

process = Tools_miniAOD.prepare( runOnGen, runOnData, eleEffAreaTarget, verbosityLvl )

# source
process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    duplicateCheckMode = cms.untracked.string( "noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #'/store/mc/Spring14dr/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/124EBB03-F1E6-E311-9837-002590A8DC50.root'
        #'/store/cmst3/user/gpetrucc/miniAOD/v1/DYJetsToLL_M-50_13TeV-madgraph-pythia8_Flat20to50_PAT.root'
        #'file://WprimeTauMiniAOD.root'
        # '/store/mc/Spring14miniaod/WprimeToMuNu_M_5800_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/E4AD1244-2809-E411-919F-0025904B1452.root'
        'file:///disk1/erdweg/MINIAOD_files/0603D444-2D70-E411-AF03-002618943922.root'
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

###FIXME for miniAOD!!!
#if not runOnGen:
    #if not runOnData:
        # This is used for QCD samples only, and the filters only write a flag in the
        # event if they fired or not. The actual selection must happen in the
        # classification, i.e. you have to set in the config file which flag you want to
        # consider.
        #
        # remove events with electrons that come from b or c hadrons
        #Tools_miniAOD.addBCtoEFilter( process )

        # remove events containing b quarks
        #Tools_miniAOD.addBFilter( process )

        # remove events with at least one potential electron candidate
        #Tools_miniAOD.addEMFilter( process )

        # remove events with a muon of more than 5 GeV
        #Tools_miniAOD.addMuGenFilter( process, pt = 5 )

        # remove events with a muon of more than 10 GeV
        #Tools_miniAOD.addMuGenFilter( process, pt = 10 )

        # remove events with a muon of more than 15 GeV
        #Tools_miniAOD.addMuGenFilter( process, pt = 15 )


print 'INFO: Using global tag:', process.GlobalTag.globaltag
