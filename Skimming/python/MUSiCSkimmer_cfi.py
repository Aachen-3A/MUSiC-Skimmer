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
    UseSIM = cms.untracked.bool( False ),
    # name of the LHgrid for pdf weights
    #LHgridName = cms.untracked.string("cteq61.LHgrid"),
    LHgridName = cms.untracked.string("cteq61.LHgrid"),
    # number of pdf error sets in the LHgrid for pdf weights
    NumLHgridErrorSets = cms.untracked.int32(40),
    #labels of source
    genParticleCandidatesLabel = cms.untracked.string( "genParticles" ),
    #vertices with beam spot constraint
    VertexRecoLabel = cms.untracked.string("offlinePrimaryVerticesWithBS"),
    #the following is all PAT
    TauRecoLabel = cms.untracked.string("hpsPFTauProducer"),
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
                ),
        pfMETNoPU = cms.PSet(
                MCLabel = cms.InputTag( 'genMetCalo' ),
                RecoLabel = cms.InputTag( 'patMETsPFlowNoPU' )
                )
        ),
    jets = cms.PSet(
        # REMARK: The names of the following PSets will be used as the names for the PXL particles that are the jets
        AK5 = cms.PSet(
            MCLabel = cms.InputTag( "ak5GenJets" ),
            RecoLabel = cms.InputTag( "cleanPatJets" ),
            isPF = cms.bool(False),
            # the following vector must hold the names of the IDs in the same sequence
            # as the qualities in PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h
            IDs = cms.vstring( 'MINIMAL', 'LOOSE_AOD', 'LOOSE', 'TIGHT' )
            ),
        pfAK5 = cms.PSet(
            MCLabel = cms.InputTag( "ak5GenJets" ),
            RecoLabel = cms.InputTag( "selectedPatJetsPFlow" ),
            isPF = cms.bool(True),
            # the following vector must hold the names of the IDs in the same sequence
            # as the qualities in PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
            IDs = cms.vstring( 'LOOSE', 'TIGHT' )
            ),
        ),
    
    triggers = cms.PSet(
        #REMARK: The names of the following PSets will be used as the trigger identifier in the PXL output
        HLT = cms.PSet(
            process = cms.string( 'auto' ),
            L1_result = cms.InputTag( "gtDigis" ),
            results = cms.string('TriggerResults'),
            event   = cms.string('hltTriggerSummaryAOD'),
            HLTriggers = cms.vstring(
                # MC HLTs for Summer11 MCs
                # (using HLT config: /cdaq/physics/Run2011/5e32/v6.2/HLT/V1)
                #
                'HLT_Mu15_v2',
                'HLT_Mu20_v1',
                'HLT_Mu24_v1',
                'HLT_Mu30_v1',
                'HLT_IsoMu12_v1',
                'HLT_IsoMu15_v5',
                'HLT_IsoMu17_v5',
                'HLT_IsoMu24_v1',
                'HLT_IsoMu30_v1',

                'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',
                'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                'HLT_Ele45_CaloIdVT_TrkIdT_v2',

                'HLT_Photon50_CaloIdVL_IsoL_v1',
                'HLT_Photon75_CaloIdVL_v2',
                'HLT_Photon75_CaloIdVL_IsoL_v2',

                'HLT_Jet240_v1',
                'HLT_Jet370_v1',
                'HLT_Jet370_NoJetID_v1',

                'HLT_MET200_v1',

                'HLT_IsoPFTau35_Trk20_MET45_v2',
                'HLT_Mu15_LooseIsoPFTau20_v2',
                'HLT_IsoMu12_LooseIsoPFTau10_v2',
                'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2'

                # Default data triggers. The actual once are defined in data*cfg.py
                #
                'HLT_Mu15_v2', 'HLT_Mu20_v1', 'HLT_IsoMu12_v1', 'HLT_IsoMu15_v5',
                'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2', 'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1',
                'HLT_Photon50_CaloIdVL_IsoL_v1', 'HLT_Photon75_CaloIdVL_IsoL_v2',
                'HLT_Jet240_v1', 'HLT_Jet370_v1',
                'HLT_MET200_v1'
                )
            ),
        StoreL3Objects = cms.untracked.bool(False)
        ),

    # This is used to access the results of all filters that ran.
    #
    filters = cms.PSet(
        AllFilters = cms.PSet(
            process = cms.string( 'PAT' ),
            results = cms.string( 'TriggerResults' ),
            paths = cms.vstring()
        )
    ),

    cuts = cms.PSet(
        min_tau_pt  = cms.double( 10 ),
        min_muon_pt = cms.double( 5 ),
        min_ele_pt = cms.double( 5 ),
        min_gamma_pt = cms.double( 5 ),
        min_jet_pt = cms.double( 20 ),
        min_met = cms.double( 0 ),
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
