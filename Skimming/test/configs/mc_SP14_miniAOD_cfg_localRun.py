runOnData = False
#run on GEN sample
runOnGen = False

import FWCore.ParameterSet.Config as cms
import sys
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

import PxlSkimmer.Skimming.Tools_miniAOD as Tools_miniAOD

print sys.argv
name="test"
datasetpath="dummy"
globalTag="MCRUN2_74_V9"
for option in sys.argv:
    splitoption=option.split('=')
    if "name" in option and len(splitoption) > 1:
        name = splitoption[1]
    if "datasetpath" in option and len(splitoption) > 1:
        datasetpath = splitoption[1]
    if "globalTag" in option and len(splitoption) > 1:
        globalTag = splitoption[1]

print "music_crab3 name %s"%name
print "music_crab3 datasetpath %s"%datasetpath
print "music_crab3 globalTag %s"%globalTag

process = Tools_miniAOD.prepare( runOnGen, runOnData, eleEffAreaTarget, name, datasetpath, globalTag, verbosityLvl )

# source
process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    duplicateCheckMode = cms.untracked.string( "noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
        #'/store/mc/Spring14dr/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/124EBB03-F1E6-E311-9837-002590A8DC50.root'
        #'/store/cmst3/user/gpetrucc/miniAOD/v1/DYJetsToLL_M-50_13TeV-madgraph-pythia8_Flat20to50_PAT.root'
        #'file://WprimeTauMiniAOD.root'
        '/store/mc/RunIISpring15DR74/WprimeToMuNu_M-1600_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/A449BF1F-DCFC-E411-9350-00259073E442.root'
        # 'file:///disk1/erdweg/MINIAOD_files/0603D444-2D70-E411-AF03-002618943922.root'
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
