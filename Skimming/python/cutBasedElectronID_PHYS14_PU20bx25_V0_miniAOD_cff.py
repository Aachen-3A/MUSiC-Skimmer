import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry

ebCutOff = 1.479

#
# This is the PRELIMINARY version of PHYS14 cuts, it is a duplicate of the corresponding
# CSA14 scenario, as explained below. The only change: the class for full 5x5 sigma_ieta_ieta
# used here is different since this variable can now be accessed directly from the
# electron object.
#
# The ID cuts below are optimized IDs for CSA14 Scenario 2 (25ns, PU 20)
# The cut values are taken from the twiki:
#       https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
#       (where they may not stay, if a newer version of cuts becomes available for these
#        conditions)
# See also the presentation explaining these working points (this will not change):
#        https://indico.cern.ch/event/292929/contribution/3/material/slides/0.pdf

cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_veto = cms.PSet(
    idName = cms.string("cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-veto"),
    cutFlow = cms.VPSet(
        cms.PSet( cutName = cms.string("MinPtCut"),
                  minPt = cms.double(5.0),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)                ),
        cms.PSet( cutName = cms.string("GsfEleSCEtaMultiRangeCut"),
                  useAbsEta = cms.bool(True),
                  allowedEtaRanges = cms.VPSet(
                                  cms.PSet( minEta = cms.double(0.0),
                                            maxEta = cms.double(ebCutOff) ),
                                  cms.PSet( minEta = cms.double(ebCutOff),
                                            maxEta = cms.double(2.5) )
                                  ),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDEtaInCut'),
                  dEtaInCutValueEB = cms.double(0.0200),
                  dEtaInCutValueEE = cms.double(0.0141),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDPhiInCut'),
                  dPhiInCutValueEB = cms.double(0.258),
                  dPhiInCutValueEE = cms.double(0.259),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut72X'),
                  full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0125),
                  full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0371),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleHadronicOverEMCut'),
                  hadronicOverEMCutValueEB = cms.double(0.256),
                  hadronicOverEMCutValueEE = cms.double(0.134),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDxyCut'),
                  dxyCutValueEB = cms.double(0.025),
                  dxyCutValueEE = cms.double(0.223),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDzCut'),
                  dzCutValueEB = cms.double(0.586),
                  dzCutValueEE = cms.double(0.951),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                  eInverseMinusPInverseCutValueEB = cms.double(0.151),
                  eInverseMinusPInverseCutValueEE = cms.double(0.154),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDeltaBetaIsoCutStandalone'),
                  isoCutEBLowPt = cms.double(0.331),
                  isoCutEBHighPt = cms.double(0.331),
                  isoCutEELowPt = cms.double(0.382),
                  isoCutEEHighPt = cms.double(0.382),
                  isRelativeIso = cms.bool(True),
                  deltaBetaConstant = cms.double(0.5),
                  ptCutOff = cms.double(20.0),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleConversionVetoCut'),
                  conversionSrc = cms.InputTag('reducedEgamma:reducedConversions'),
                  beamspotSrc = cms.InputTag('offlineBeamSpot'),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleMissingHitsCut'),
                  maxMissingHitsEB = cms.uint32(2),
                  maxMissingHitsEE = cms.uint32(3),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False))
    )
)

cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_loose = cms.PSet(
    idName = cms.string("cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-loose"),
    cutFlow = cms.VPSet(
        cms.PSet( cutName = cms.string("MinPtCut"),
                  minPt = cms.double(5.0),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)                ),
        cms.PSet( cutName = cms.string("GsfEleSCEtaMultiRangeCut"),
                  useAbsEta = cms.bool(True),
                  allowedEtaRanges = cms.VPSet(
                                  cms.PSet( minEta = cms.double(0.0),
                                            maxEta = cms.double(ebCutOff) ),
                                  cms.PSet( minEta = cms.double(ebCutOff),
                                            maxEta = cms.double(2.5) )
                                  ),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDEtaInCut'),
                  dEtaInCutValueEB = cms.double(0.0181),
                  dEtaInCutValueEE = cms.double(0.0124),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDPhiInCut'),
                  dPhiInCutValueEB = cms.double(0.0936),
                  dPhiInCutValueEE = cms.double(0.0642),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut72X'),
                  full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0123),
                  full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0350),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleHadronicOverEMCut'),
                  hadronicOverEMCutValueEB = cms.double(0.141),
                  hadronicOverEMCutValueEE = cms.double(0.112),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDxyCut'),
                  dxyCutValueEB = cms.double(0.0166),
                  dxyCutValueEE = cms.double(0.098),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDzCut'),
                  dzCutValueEB = cms.double(0.543),
                  dzCutValueEE = cms.double(0.919),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                  eInverseMinusPInverseCutValueEB = cms.double(0.135),
                  eInverseMinusPInverseCutValueEE = cms.double(0.144),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDeltaBetaIsoCutStandalone'),
                  isoCutEBLowPt = cms.double(0.24),
                  isoCutEBHighPt = cms.double(0.24),
                  isoCutEELowPt = cms.double(0.353),
                  isoCutEEHighPt = cms.double(0.353),
                  isRelativeIso = cms.bool(True),
                  deltaBetaConstant = cms.double(0.5),
                  ptCutOff = cms.double(20.0),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleConversionVetoCut'),
                  conversionSrc = cms.InputTag('reducedEgamma:reducedConversions'),
                  beamspotSrc = cms.InputTag('offlineBeamSpot'),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleMissingHitsCut'),
                  maxMissingHitsEB = cms.uint32(1),
                  maxMissingHitsEE = cms.uint32(1),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False))
    )
)

cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_medium = cms.PSet(
    idName = cms.string("cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-medium"),
    cutFlow = cms.VPSet(
        cms.PSet( cutName = cms.string("MinPtCut"),
                  minPt = cms.double(5.0),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)                ),
        cms.PSet( cutName = cms.string("GsfEleSCEtaMultiRangeCut"),
                  useAbsEta = cms.bool(True),
                  allowedEtaRanges = cms.VPSet(
                                  cms.PSet( minEta = cms.double(0.0),
                                            maxEta = cms.double(ebCutOff) ),
                                  cms.PSet( minEta = cms.double(ebCutOff),
                                            maxEta = cms.double(2.5) )
                                  ),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDEtaInCut'),
                  dEtaInCutValueEB = cms.double(0.0106),
                  dEtaInCutValueEE = cms.double(0.0108),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDPhiInCut'),
                  dPhiInCutValueEB = cms.double(0.0323),
                  dPhiInCutValueEE = cms.double(0.0455),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut72X'),
                  full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0107),
                  full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0318),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleHadronicOverEMCut'),
                  hadronicOverEMCutValueEB = cms.double(0.067),
                  hadronicOverEMCutValueEE = cms.double(0.097),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDxyCut'),
                  dxyCutValueEB = cms.double(0.0131),
                  dxyCutValueEE = cms.double(0.0845),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDzCut'),
                  dzCutValueEB = cms.double(0.223),
                  dzCutValueEE = cms.double(0.752),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                  eInverseMinusPInverseCutValueEB = cms.double(0.104),
                  eInverseMinusPInverseCutValueEE = cms.double(0.120),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDeltaBetaIsoCutStandalone'),
                  isoCutEBLowPt = cms.double(0.218),
                  isoCutEBHighPt = cms.double(0.218),
                  isoCutEELowPt = cms.double(0.254),
                  isoCutEEHighPt = cms.double(0.254),
                  isRelativeIso = cms.bool(True),
                  deltaBetaConstant = cms.double(0.5),
                  ptCutOff = cms.double(20.0),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleConversionVetoCut'),
                  conversionSrc = cms.InputTag('reducedEgamma:reducedConversions'),
                  beamspotSrc = cms.InputTag('offlineBeamSpot'),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleMissingHitsCut'),
                  maxMissingHitsEB = cms.uint32(1),
                  maxMissingHitsEE = cms.uint32(1),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False))
    )
)

cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_tight = cms.PSet(
    idName = cms.string("cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-tight"),
    cutFlow = cms.VPSet(
        cms.PSet( cutName = cms.string("MinPtCut"),
                  minPt = cms.double(5.0),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)                ),
        cms.PSet( cutName = cms.string("GsfEleSCEtaMultiRangeCut"),
                  useAbsEta = cms.bool(True),
                  allowedEtaRanges = cms.VPSet(
                                  cms.PSet( minEta = cms.double(0.0),
                                            maxEta = cms.double(ebCutOff) ),
                                  cms.PSet( minEta = cms.double(ebCutOff),
                                            maxEta = cms.double(2.5) )
                                  ),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDEtaInCut'),
                  dEtaInCutValueEB = cms.double(0.0091),
                  dEtaInCutValueEE = cms.double(0.0106),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDPhiInCut'),
                  dPhiInCutValueEB = cms.double(0.0310),
                  dPhiInCutValueEE = cms.double(0.0359),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut72X'),
                  full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0106),
                  full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0305),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleHadronicOverEMCut'),
                  hadronicOverEMCutValueEB = cms.double(0.0532),
                  hadronicOverEMCutValueEE = cms.double(0.0835),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDxyCut'),
                  dxyCutValueEB = cms.double(0.0126),
                  dxyCutValueEE = cms.double(0.0163),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDzCut'),
                  dzCutValueEB = cms.double(0.0116),
                  dzCutValueEE = cms.double(0.60),
                  vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                  eInverseMinusPInverseCutValueEB = cms.double(0.0609),
                  eInverseMinusPInverseCutValueEE = cms.double(0.113),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleDeltaBetaIsoCutStandalone'),
                  isoCutEBLowPt = cms.double(0.165),
                  isoCutEBHighPt = cms.double(0.165),
                  isoCutEELowPt = cms.double(0.208),
                  isoCutEEHighPt = cms.double(0.208),
                  isRelativeIso = cms.bool(True),
                  deltaBetaConstant = cms.double(0.5),
                  ptCutOff = cms.double(20.0),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleConversionVetoCut'),
                  conversionSrc = cms.InputTag('reducedEgamma:reducedConversions'),
                  beamspotSrc = cms.InputTag('offlineBeamSpot'),
                  needsAdditionalProducts = cms.bool(True),
                  isIgnored = cms.bool(False)),
        cms.PSet( cutName = cms.string('GsfEleMissingHitsCut'),
                  maxMissingHitsEB = cms.uint32(1),
                  maxMissingHitsEE = cms.uint32(1),
                  barrelCutOff = cms.double(ebCutOff),
                  needsAdditionalProducts = cms.bool(False),
                  isIgnored = cms.bool(False))
    )
)

central_id_registry.register(cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_veto.idName,
                             'ee73e19af1ec4bf4f71f46013791e98b')
central_id_registry.register(cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_loose.idName,
                             'd6fc42108cefce6e0177c80b40860382')
central_id_registry.register(cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_medium.idName,
                             '7491d788175df497b06f5922579d2ae2')
central_id_registry.register(cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_standalone_tight.idName,
                             '2fec01235da5321c4116ddef9a8bd0ff')