import FWCore.ParameterSet.Config as cms
import PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi as electronMatch
import PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi as muonMatch
import PhysicsTools.PatAlgos.mcMatchLayer0.photonMatch_cfi as photonMatch
import PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi as jetMatch

#configure PAT matching
electronMatch.checkCharge = cms.bool(False),
electronMatch.resolveByMatchQuality = cms.bool(True),
electronMatch.maxDeltaR = cms.double(0.2),  
electronMatch.maxDPtRel = cms.double(1000000.0),
muonMatch.checkCharge = cms.bool(False),
muonMatch.resolveByMatchQuality = cms.bool(True),
muonMatch.maxDeltaR = cms.double(0.2),  
muonMatch.maxDPtRel = cms.double(1000000.0),
photonMatch.checkCharge = cms.bool(False),
photonMatch.resolveByMatchQuality = cms.bool(True),
photonMatch.maxDeltaR = cms.double(0.2),  
photonMatch.maxDPtRel = cms.double(1000000.0), 
jetMatch.checkCharge = cms.bool(False),
jetMatch.resolveByMatchQuality = cms.bool(True),
jetMatch.maxDeltaR = cms.double(0.4),  
jetMatch.maxDPtRel = cms.double(1000000.0), 
