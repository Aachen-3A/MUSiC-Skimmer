import FWCore.ParameterSet.Config as cms

Skimmer = cms.EDAnalyzer(
    "PxlSkimmer_miniAOD",
    #output file name
    FileName =  cms.untracked.string("test_run.pxlio"),
    #symbolic name of the processed data
    Process = cms.untracked.string("test_run"),
    # Dataset of the processed data or MC.
    Dataset = cms.untracked.string( 'test_run' ),
    # GenOnly true mean no Rec-info in event, check for GenJets and GenMET
    GenOnly = cms.untracked.bool( False ),
    # Are we running on a FASTSIM sample?
    FastSim = cms.bool( False ),
    # UseSIM true means to use SIM info for finding converted photons
    UseSIM = cms.untracked.bool( False ),
    # name of the LHgrid for pdf weights
    LHgridName = cms.untracked.string("cteq61.LHgrid"),
    # number of pdf error sets in the LHgrid for pdf weights
    NumLHgridErrorSets = cms.untracked.int32(40),
    #labels of source
    genParticleCandidatesLabel  = cms.InputTag( "prunedGenParticles" ),
    genFinalParticlesLabel      = cms.InputTag( "packedGenParticles" ),
    #vertices with beam spot constraint
    VertexRecoLabel     = cms.untracked.string("offlineSlimmedPrimaryVertices"),
    #the following is all PAT
    patMuonLabel       = cms.InputTag("slimmedMuons"),
    patElectronLabel   = cms.InputTag("slimmedElectrons"),
    patTauTag          = cms.InputTag( 'slimmedTaus' ),
    patGammaLabel      = cms.InputTag("slimmedPhotons"),
    patMETTag          = cms.InputTag( 'slimmedMETs' ),
    PUPPIMETTag        = cms.InputTag( 'slimmedMETsPuppi' ),
    noHFMETTag        = cms.InputTag( 'slimmedMETsNoHF' ),
    newUncertMETTag        = cms.InputTag( 'slimmedMETsnewUncert' ),
    patPFCandiates     = cms.InputTag( 'packedPFCandidates' ),

    rhos = cms.VInputTag( cms.InputTag( 'fixedGridRhoAll' ),
                       cms.InputTag( 'fixedGridRhoFastjetAll' ),
                       cms.InputTag( 'fixedGridRhoFastjetAllCalo' ),
                       cms.InputTag( 'fixedGridRhoFastjetCentralCalo' ),
                       cms.InputTag( 'fixedGridRhoFastjetCentralChargedPileUp' ),
                       cms.InputTag( 'fixedGridRhoFastjetCentralNeutral' )
                       ),

    eleIDs = cms.VInputTag(
                            #have the 50ns just in case!!
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-veto"),
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-loose"),
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-medium"),
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V1-standalone-tight"),

                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                            cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                            cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                            cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                            cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60")
                           ),

    gammaIDs = cms.VInputTag(   cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                                cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
                                cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
                                cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
                                cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"),
                                cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight"),
                             ),

    bits = cms.InputTag("TriggerResults","","HLT"),
    objects = cms.InputTag("selectedPatTrigger"),

    #ECAL RecHits for supercluster information
    reducedSuperClusterCollection   = cms.InputTag("reducedEgamma","reducedESClusters"),
    reducedEBClusterCollection      = cms.InputTag("reducedEgamma","reducedEBEEClusters"),
    #HCAL noise
    HCALNoise                      = cms.InputTag( "HBHENoiseFilterResultProducer", "HBHENoiseFilterResult" ),
    METFilterTag                   = cms.InputTag("TriggerResults","","PAT"),

    conversionsTag                  = cms.InputTag( "reducedEgamma","reducedConversions" ),
    conversionsSingleLegTag         = cms.InputTag("reducedEgamma","reducedSingleLegConversions" ),

    jets = cms.PSet(
        # REMARK: The names of the following PSets will be used as the names for the PXL particles that are the jets
        AK4Corr = cms.PSet(
            MCLabel = cms.InputTag( "slimmedGenJets" ),
            RecoLabel = cms.InputTag( "patJetsUpdated" ),
            isPF = cms.bool(True),
            # the following vector must hold the names of the IDs in the same sequence
            # as the qualities in PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h
            IDs = cms.vstring( 'LOOSE', 'TIGHT' )
            ),
        AK4 = cms.PSet(
            MCLabel = cms.InputTag( "slimmedGenJets" ),
            RecoLabel = cms.InputTag( "slimmedJets" ),
            isPF = cms.bool(True),
            # the following vector must hold the names of the IDs in the same sequence
            # as the qualities in PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h
            IDs = cms.vstring( 'LOOSE', 'TIGHT' )
            ),
        AK8 = cms.PSet(
            MCLabel = cms.InputTag( "" ),
            RecoLabel = cms.InputTag( "slimmedJetsAK8" ),
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

            # A list of triggers can be defined in the actual *cfg.py.
            # Otherwise all unprescaled triggers from the HLT config will
            # be used for each run.
            HLTriggers = cms.vstring(),
            # Only triggers from datastreams whose name is given in the following list
            # will be considered. Make sure to update this list regularly.
            #datastreams = cms.vstring(
                                        #"InitialPD",
                                        #"Templates"
            #),
             datastreams = cms.vstring(
                                         "BJetPlusX",
                                         "BTag",
                                         #"Commissioning",
                                         #"Cosmics",
                                         "DoubleElectron",
                                         "DoubleMu",
                                         "DoubleMuParked",
                                         "DoublePhoton",
                                         "DoublePhotonHighPt",
                                         "ElectronHad",
                                         #"FEDMonitor",
                                         #"HLTPhysicsParked",
                                         "HTMHT",
                                         "HTMHTParked",
                                         #"HcalHPDNoise",
                                         #"HcalNZS",
                                         "JetHT",
                                         "JetMon",
                                         #"LogMonitor",
                                         "MET",
                                         "METParked",
                                         "MinimumBias",
                                         "EGamma",
                                         "MuHad",
                                         #"MuOnia",
                                         #"MuOniaParked",
                                         "MultiJet",
                                         "MultiJet1Parked",
                                         "NoBPTX",
                                         "PPMuon",
                                         "PPPhoton",
                                         "PhotonHad",
                                         "SingleElectron",
                                         "SingleMu",
                                         "SinglePhoton",
                                         "SinglePhotonParked",
                                         "Tau",
                                         "TauParked",
                                         "TauPlusX",
                                         "VBF1Parked",
                                         #"ZeroBiasParked",
             ),
        ),
        StoreL3Objects = cms.bool(False)
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


    # "borrowed" from IB
    # for miniIsolation PF-weighted isolation
    # https://github.com/SUSYDileptonAC/SuSyAachen/blob/master/TagAndProbeTreeWriter/python/isolationFunctor_cfi.py
    # this version:
    # https://github.com/SUSYDileptonAC/SuSyAachen/commit/e9f0eaf2fd3828eabe18545d88fa83a9a8eb72f5
    isolationDefinitions = cms.PSet(
        rhoSource = cms.InputTag("fixedGridRhoFastjetAll"),
        candSource = cms.InputTag("packedPFCandidates"),
    ),

    cuts = cms.PSet(
        min_tau_pt  = cms.double( 10 ),
        min_muon_pt = cms.double( 5 ),
        min_ele_pt = cms.double( 5 ),
        min_gamma_pt = cms.double( 5 ),
        min_jet_pt = cms.double( 20 ),
        min_met = cms.double( 10 ),
        max_eta = cms.double( 3 ),
        min_rechit_energy = cms.double( 20 ),
        min_rechit_swiss_cross = cms.double( 0.8 ),
        min_rechit_R19 = cms.double( 0.8 ),
        vertex_minNDOF = cms.double( 3 ),
        vertex_maxZ = cms.double( 30 ),
        vertex_maxR = cms.double( 3 ),
        # These cuts come from:
        # https://twiki.cern.ch/twiki/bin/viewauth/CMSPublic/WorkBookChapter8?rev=27
        # See also:
        # CMS PAS TRK-10-005
        PV_minNDOF = cms.double( 4 ),
        PV_maxZ = cms.double( 24 ),
        PV_maxR = cms.double( 2 )
    )
)
