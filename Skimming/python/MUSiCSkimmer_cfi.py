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
    #vertices with beam spot constraint
    VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
    #the following is all PAT
    MuonRecoLabel = cms.untracked.string("cleanPatMuons"),
    ElectronRecoLabel = cms.untracked.string("cleanPatElectrons"),
    GammaRecoLabel = cms.untracked.string("cleanPatPhotons"),
    #ECAL RecHits for supercluster information
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    #HCAL noise
    HCALNoise = cms.InputTag( 'HBHENoiseFilterResultProducer', 'HBHENoiseFilterResult' ),

    METs = cms.PSet(
        # REMARK: The names of the following PSets will be used as the names for the PXL particles that are the METs
        caloMET = cms.PSet(
                MCLabel = cms.InputTag( "genMetCalo" ),
                RecoLabel = cms.InputTag("patMETs")
                ),
        pfMET = cms.PSet(
                MCLabel = cms.InputTag( "genMetCalo" ),
                RecoLabel = cms.InputTag("patMETsPFlow") #patMETsPFlow")
                )
        ),
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
                'HLT_Mu5','HLT_Mu9',
                'HLT_Photon15_L1R','HLT_Photon15_Cleaned_L1R','HLT_Photon20_L1R','HLT_Photon20_Cleaned_L1R','HLT_Photon30_L1R','HLT_Photon30_L1R_8E29',
                'HLT_Ele15_SW_L1R','HLT_Ele15_LW_L1R','HLT_Ele20_SW_L1R','HLT_Ele20_LW_L1R',
                'HLT_Jet50U','HLT_BTagMu_Jet10U','HLT_BTagIP_Jet50U',
                'HLT_MET45','HLT_MET100'
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
        min_rechit_energy = cms.double( 20 ),
        min_rechit_swiss_cross = cms.double( 0.8 ),
        min_rechit_R19 = cms.double( 0.8 ),
        vertex_minNDOF = cms.double( 3 ),
        vertex_maxZ = cms.double( 30 ),
        vertex_maxR = cms.double( 3 ),
        PV_minNDOF = cms.double( 5 ),
        PV_maxZ = cms.double( 24 ),
        PV_maxR = cms.double( 2 )
        )
    )
