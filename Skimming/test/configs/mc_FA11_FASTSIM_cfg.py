runOnData = False
# Running on a GEN sample ?
runOnGen = False
# Running on a FASTSIM sample?
runOnFast = True

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

process = Tools.prepare( runOnGen, runOnData, eleEffAreaTarget, verbosityLvl, runOnFast=runOnFast )

# source
process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    fileNames = cms.untracked.vstring(
        'file:////user/papacz/production/run/Hadronizer_MgmMatchTuneZ2star_7TeV_madgraph_tauola_cff_py_GEN_FASTSIM_HLT_PU.root'
        )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 10 ) )

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

    if runOnFast:
        # In FASTSIM, the triggers are not assigned to any datastreams.
        # We leave them empty and thus all (unprescaled) triggers are written
        # into the pxlio file.
        # This is not very space efficient, but we do not plan do many samples
        # with FASTSIM.
        process.Skimmer.triggers.HLT.datastreams = cms.vstring()

print 'INFO: Using global tag:', process.GlobalTag.globaltag
