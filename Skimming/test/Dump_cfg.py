import FWCore.ParameterSet.Config as cms

process = cms.Process("Dump")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
     '/store/mc/Summer08/TTJets-madgraph/GEN-SIM-RECO/IDEAL_V9_v2/0002/00437D8F-04A9-DD11-8D99-0015C5E9B2AB.root'
     )
)

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(process.content)
