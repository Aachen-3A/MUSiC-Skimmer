import FWCore.ParameterSet.Config as cms

process = cms.Process("Dump")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
     '/store/relval/CMSSW_2_1_9/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v2/0000/04419036-F385-DD11-B3A7-001617C3B6E8.root'
     )
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(process.content)
