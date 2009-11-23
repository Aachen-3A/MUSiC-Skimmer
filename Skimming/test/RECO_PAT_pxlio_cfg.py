runOnData = False

import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    fileMode = cms.untracked.string( 'NOMERGE' )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# source
process.source = cms.Source("PoolSource", 
     skipEvents = cms.untracked.uint32(0),
     fileNames = cms.untracked.vstring(
#'dcap://grid-dcache.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/antonius/test/photon_jets_48878A79-CCBC-DD11-85F3-0022199A2E95.root'
'/store/mc/Summer09/TTbar/GEN-SIM-RECO/MC_31X_V3_test_production-v1/0000/962E4ECF-077C-DE11-871F-001F2907AF5C.root'
#'file:/opt/scratch/608D435D-7B70-DE11-B093-001CC4782AF8.root'
	)
)


process.load("Configuration/StandardSequences/GeometryPilot2_cff")

# rerun photon ID
process.load("RecoEgamma.PhotonIdentification.photonId_cff")
# for ClusterShape inspired by egamma hypernews
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')
process.load("Configuration/StandardSequences/MagneticField_38T_cff")

# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.content = cms.EDAnalyzer("EventContentAnalyzer")


# Remove unneccessary stuff:
if runOnData:
    process.patHighLevelReco_withoutPFTau.remove( process.patJetFlavourId )
    process.patDefaultSequence_withoutTrigMatch.remove( process.patMCTruth )
    process.allLayer1Muons.addGenMatch        = False
    process.allLayer1Jets.addGenPartonMatch   = False
    process.allLayer1Jets.addGenJetMatch      = False
    process.allLayer1Jets.getJetMCFlavour     = False
    process.allLayer1METs.addGenMET           = False



# this might be commented in in order to safe the edm root file containing the PAT Products
# Output module configuration
# process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
#    # save only events passing the full path
#    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#    outputCommands = cms.untracked.vstring('drop *')
# )
# process.outpath = cms.EndPath(process.out)
# save PAT Layer 1 output
# from PhysicsTools.PatAlgos.patEventContent_cff import *
# process.out.outputCommands += patEventContent
# process.out.outputCommands.extend(["keep *_selectedLayer1Jets*_*_*"])

# insert following lines in order to get info about the HLT content (don't forget to include into path)
# import HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi
# process.hltAnalyzer = HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi.hltEventAnalyzerAOD.clone()


import os
cmsbase = os.environ.get('CMSSW_BASE')
execfile(cmsbase + "/src/MUSiCProject/Skimming/python/configurePAT_cff")

#process.pTHat = cms.EDFilter("PtHatFilter",
#         pt_hat_lower_bound = cms.double(200.),
#         pt_hat_upper_bound = cms.double(500.)
#)

# cut on the factorization scale e.g. suitable for cuts on inv. mass of resonances in Pythia!
#process.FacScale = cms.EDFilter("FacScaleFilter",
#         fac_scale_lower_bound = cms.double(200.),
#         fac_scale_upper_bound = cms.double(500.)
#)


process.load( "MUSiCProject.Skimming.MUSiCSkimmer_cfi" )


process.p = cms.Path( process.photonIDSequence + process.patDefaultSequence + process.Skimmer )
