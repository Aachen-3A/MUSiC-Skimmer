import FWCore.ParameterSet.Config as cms

Skimmer = cms.EDAnalyzer(
    "MUSiCSkimmer",
    #output file name
    FileName =  cms.untracked.string("test_run.pxlio"),
    #symbolic name of the processed data
    Process = cms.untracked.string("test_run"),
    # Debugging: 0 = off, 1 = human readable, 2 = insane
    debug = cms.untracked.int32( 0 ),
    # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
    GenOnly = cms.untracked.bool( False ),
    # UseSIM true means to use SIM info for finding converted photons
    UseSIM = cms.untracked.bool( True ),
    #labels of source
    genParticleCandidatesLabel = cms.untracked.string( "genParticles" ),
    METMCLabel = cms.untracked.string( "genMetCalo" ),  # muon-correction needed ---> yes!
    #vertices with beam spot constraint
    VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
    #the following is all PAT
    MuonRecoLabel = cms.untracked.string("cleanPatMuons"),
    ElectronRecoLabel = cms.untracked.string("cleanPatElectrons"),
    GammaRecoLabel = cms.untracked.string("cleanPatPhotons"),
    METRecoLabel = cms.untracked.string("patMETs"),
    #ECAL RecHits for supercluster information
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    #HCAL noise
    HCALNoiseInfo = cms.InputTag( 'hcalnoise' ),
    
    jets = cms.PSet(
        # REMARK: The names of the following PSets will be used as the names for the PXL particles that are the jets
        AK5 = cms.PSet(
            MCLabel = cms.InputTag( "ak5GenJets" ),
            RecoLabel = cms.InputTag( "cleanPatJets" )
            ),
        # the following vector must hold the names of the IDs in the same sequence
        # as the qualities in PhysicsTools/PatUtils/interface/JetIDSelectionFunctor.h
        IDs = cms.vstring( 'minimal', 'loose_aod', 'loose', 'tight' )
        ),
    
    triggers = cms.PSet(
        #REMARK: The names of the following PSets will be used as the trigger identifier in the PXL output
        # Trigger menu: 8e29
        HLT = cms.PSet(
            process = cms.string('HLT'),
            L1_result = cms.InputTag( "gtDigis" ),
            results = cms.string('TriggerResults'),
            event   = cms.string('hltTriggerSummaryAOD'),
            HLTriggers = cms.vstring(
                'HLT_L1Mu20','HLT_L2Mu9','HLT_L2Mu11','HLT_L1Mu14_L1SingleEG10','HLT_L1Mu14_L1SingleJet6U','HLT_L1Mu14_L1ETM30','HLT_L2DoubleMu0','HLT_L1DoubleMuOpen','HLT_DoubleMu0','HLT_DoubleMu3','HLT_Mu3','HLT_Mu5','HLT_Mu9','HLT_IsoMu3','HLT_Mu0_L1MuOpen','HLT_Mu0_Track0_Jpsi','HLT_Mu3_L1MuOpen','HLT_Mu3_Track0_Jpsi','HLT_Mu5_L1MuOpen','HLT_Mu5_Track0_Jpsi','HLT_Mu0_L2Mu0','HLT_Mu3_L2Mu0','HLT_Mu5_L2Mu0',
                'HLT_Photon15_L1R','HLT_Photon15_LooseEcalIso_L1R','HLT_Photon20_L1R','HLT_Photon30_L1R_8E29','HLT_DoublePhoton4_Jpsi_L1R','HLT_DoublePhoton4_Upsilon_L1R','HLT_DoublePhoton5_Jpsi_L1R','HLT_DoublePhoton5_Upsilon_L1R','HLT_DoublePhoton5_L1R','HLT_DoublePhoton10_L1R','HLT_DoubleEle5_SW_L1R','HLT_Ele20_LW_L1R','HLT_Ele15_SiStrip_L1R','HLT_Ele15_SC10_LW_L1R','HLT_Ele15_LW_L1R','HLT_Ele10_LW_EleId_L1R','HLT_Ele10_LW_L1R','HLT_Photon15_TrackIso_L1R',
                'HLT_FwdJet20U','HLT_Jet30U','HLT_Jet50U','HLT_DiJetAve30U_8E29','HLT_QuadJet15U','HLT_MET45','HLT_MET100','HLT_HT100U','HLT_SingleLooseIsoTau20','HLT_DoubleLooseIsoTau15','HLT_DoubleJet15U_ForwardBackward','HLT_BTagMu_Jet10U','HLT_BTagIP_Jet50U','HLT_StoppedHSCP_8E29'
                )
            ),
        StoreL3Objects = cms.untracked.bool(False)
        ),

    cuts = cms.PSet(
        min_muon_pt = cms.double( 5 ),
        min_ele_pt = cms.double( 5 ),
        min_gamma_pt = cms.double( 5 ),
        min_jet_pt = cms.double( 20 ),
        min_met = cms.double( 20 ),
        max_eta = cms.double( 3 ),
        vertex_minNDOF = cms.double( 3 ),
        vertex_maxZ = cms.double( 20 ),
        vertex_maxR = cms.double( 4 ),
        PV_minNDOF = cms.double( 5 ),
        PV_maxZ = cms.double( 15 ),
        PV_maxR = cms.double( 2 )
        )
    )
