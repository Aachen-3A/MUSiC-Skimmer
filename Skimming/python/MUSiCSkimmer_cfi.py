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
    TauRecoLabel = cms.untracked.string( "patTausPFlow" ),
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

            # Data triggers must be defined in the actual data*cfg.py.
            # Examples in $CMSSW_BASE/src/MUSiCProject/Skimming/test/configs
            #
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
                'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2',

                # MC HLTs for Fall11 samples:
                # ---------------------------
                # There are different production revisions that were processed
                # with different CMSSW releases. Each of these used a different
                # HLT config, so different triggers are available.
                #
                # Always using the HLT config found in 'CMSSW/HLTrigger/Configuration/python/HLT_GRun_cff.py'
                # for the respective release.
                #
                # Fall11_R1 + Fall11_R3, CMSSW_4_2_9_HLT1_patch1 + CMSSW_4_4_0_patch3,
                # /online/collisions/2011/3e33/v1.1/HLT/V10 + /dev/CMSSW_4_2_0/GRun/V226 (equivalent to /online/collisions/2011/3e33/v2.0/HLT/V11):
                #
                'HLT_Mu40_eta2p1_v1',

                'HLT_IsoMu30_eta2p1_v3',

                'HLT_Ele65_CaloIdVT_TrkIdT_v4',

                'HLT_Photon90_CaloIdVL_IsoL_v5',

                'HLT_Jet300_v5',

                'HLT_MET200_v7',

                'HLT_Mu13_Mu8_v7',

                'HLT_MediumIsoPFTau35_Trk20_MET60_v1',
                'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v3',

                'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1',
                'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v1',
                'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v1',

                'HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1',

                # Fall11_R2, CMSSW_4_2_8_patch4,
                # /online/collisions/2011/5e32/v6.2/HLT/V4:
                #
                'HLT_Mu15_v2',

                'HLT_IsoMu12_v1',

                'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2',

                'HLT_Photon50_CaloIdVL_IsoL_v1',

                'HLT_Jet240_v1',

                'HLT_MET200_v1',

                'HLT_DoubleMu6_v1',

                'HLT_DoubleIsoPFTau20_Trk5_v2',
                'HLT_IsoPFTau35_Trk20_MET45_v2',

                'HLT_Mu15_LooseIsoPFTau20_v2',
                'HLT_IsoMu12_LooseIsoPFTau10_v2',

                'HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2',

                # Fall11_R4 + Fall11_Chamonix + Fall11_TSG, CMSSW_4_4_2_patch8 + CMSSW_4_4_2_patch8 + CMSSW_4_4_2_patch6,
                # /online/collisions/2011/5e33/v3.0/HLT/V6 + /online/collisions/2011/5e33/v3.0/HLT/V4:
                #
                'HLT_Mu50_eta2p1_v2',

                'HLT_IsoMu30_eta2p1_v7',

                'HLT_Ele80_CaloIdVT_TrkIdT_v3',

                'HLT_Photon135_v2',

                'HLT_Jet370_v10',

                'HLT_MET200_v7',

                'HLT_Mu17_Mu8_v11',

                'HLT_DoubleIsoPFTau45_Trk5_eta2p1_v8',
                'HLT_MediumIsoPFTau35_Trk20_MET60_v6',

                'HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v6',
                'HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v6',
                'HLT_IsoMu15_eta2p1_TightIsoPFTau20_v6',

                'HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v6',
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
