runOnData = False
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
        '/store/mc/Summer11/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/AODSIM/PU_S4_START42_V11-v1/0000/C45712AA-A6A8-E011-8679-0024E86E8D4C.root'
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )

print 'INFO: Using global tag:', process.GlobalTag.globaltag
