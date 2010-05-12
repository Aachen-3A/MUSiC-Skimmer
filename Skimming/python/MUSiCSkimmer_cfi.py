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
                'HLT_Mu3', 'HLT_DoubleMu0',
                'HLT_Ele10_LW_L1R', 'HLT_DoubleEle5_SW_L1R',
                'HLT_Photon15_L1R', 'HLT_DoublePhoton10_L1R'
                #no usefull cross-channel trigger in this menu
                )
            ),
        StoreL3Objects = cms.untracked.bool(False)
        ),

    cuts = cms.PSet(
        min_muon_pt = cms.double( 10 ),
        min_ele_pt = cms.double( 10 ),
        min_gamma_pt = cms.double( 10 ),
        min_jet_pt = cms.double( 30 ),
        min_met = cms.double( 30 ),
        max_eta = cms.double( 3 ),
        vertex_minNDOF = cms.double( 3 ),
        vertex_maxZ = cms.double( 20 ),
        vertex_maxR = cms.double( 4 ),
        PV_minNDOF = cms.double( 5 ),
        PV_maxZ = cms.double( 15 ),
        PV_maxR = cms.double( 2 )
        )
    )
