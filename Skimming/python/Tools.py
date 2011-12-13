import FWCore.ParameterSet.Config as cms

def prepare( runOnGen, runOnData ):
    process = cms.Process( 'PAT' )

    # Initialize MessageLogger and output report.
    #
    process.load( 'FWCore.MessageLogger.MessageLogger_cfi' )
    process.MessageLogger.cerr.FwkReport.limit = 100

    import FWCore.Framework.test.cmsExceptionsFatalOption_cff
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool( True ),
        # Open file in NOMERGE mode to avoid a memory leak.
        #
        fileMode = cms.untracked.string( 'NOMERGE' ),
        # Stop processing on each and every thrown exception.
        #
        Rethrow = FWCore.Framework.test.cmsExceptionsFatalOption_cff.Rethrow
        )

    # The global tag is retrieved automatically but can be changed by the
    # configureJEC function.
    #
    process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )
    from Configuration.PyReleaseValidation.autoCond import autoCond
    if runOnData:
        process.GlobalTag.globaltag = cms.string( autoCond[ 'com10' ] )
    else:
        process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

    process.load( 'Configuration.StandardSequences.GeometryPilot2_cff' )
    process.load( 'Configuration.StandardSequences.MagneticField_38T_cff' )

    # do we need this ?
    #process.content = cms.EDAnalyzer( 'EventContentAnalyzer' )

    # Create an empty path because modules will be added by calling the
    # functions below.
    #
    process.p = cms.Path()

    process.load( 'MUSiCProject.Skimming.MUSiCSkimmer_cfi' )

    # Several filters are used while running over data or MC.
    # In order to be flexible, events *not* passing these filtes we do not want
    # to throw these events away but rather write the outcome as a bool into the
    # event.
    # To do this, each filter runs in an own path. These paths are stored in the
    # filterlist. This list is later used to access the value of the filter with
    # help of edm::TriggerResult.
    #
    process.Skimmer.filterlist = cms.vstring()

    if not runOnGen:
        addScrapingFilter( process )

        # Keep the following functions in the right order as they will add modules to the path!
        #
        configureJEC( process, runOnData )
        configurePAT( process, runOnData )
        process.metJESCorAK5CaloJet.inputUncorMetLabel = 'metNoHF'

        postfix = 'PFlow'
        configurePF( process, runOnData, postfix )
        configurePFnoPU( process, postfix )

        import PhysicsTools.PatAlgos.tools.coreTools
        PhysicsTools.PatAlgos.tools.coreTools.removeMCMatching( process, [ 'All' ] )

        addHCALnoiseFilter( process )

    if not runOnData:
       process.p += process.patJetPartons

    if not runOnData:
        addFlavourMatching( process, process.Skimmer, process.p, runOnGen )

    configureTaus( process )

    process.Skimmer.filters.AllFilters.paths = process.Skimmer.filterlist
    process.Skimmer.filters.AllFilters.process = process.name_()

    # The skimmer is in the endpath because then the results of all preceding paths
    # are available. This is used to access the outcome of filters that ran.
    #
    process.e = cms.EndPath( process.Skimmer )

    return process


def configurePAT( process, runOnData ):
    # PAT Layer 0+1
    process.load( 'PhysicsTools.PatAlgos.patSequences_cff' )

    #do not store TagInfos, as they are not in AOD
    process.patJets.addTagInfos = False
    #do not embed muon tracks, as it breaks the TeV-refit
    process.patMuons.embedCombinedMuon = False
    process.patMuons.embedStandAloneMuon = False

    if runOnData:
        #configure PAT matching
        process.electronMatch.checkCharge = False
        process.electronMatch.resolveByMatchQuality = True
        process.electronMatch.maxDeltaR = 0.2
        process.electronMatch.maxDPtRel = 1000000.0

        process.muonMatch.checkCharge = False
        process.muonMatch.resolveByMatchQuality = True
        process.muonMatch.maxDeltaR = 0.2
        process.muonMatch.maxDPtRel = 1000000.0

        process.photonMatch.checkCharge = False
        process.photonMatch.resolveByMatchQuality = True
        process.photonMatch.maxDeltaR = 0.2
        process.photonMatch.maxDPtRel = 1000000.0

        process.patJetPartonMatch.checkCharge = False
        process.patJetPartonMatch.resolveByMatchQuality = True
        process.patJetPartonMatch.maxDeltaR = 0.4
        process.patJetPartonMatch.maxDPtRel = 1000000.0

        process.patJetGenJetMatch.checkCharge = False
        process.patJetGenJetMatch.resolveByMatchQuality = True
        process.patJetGenJetMatch.maxDeltaR = 0.4
        process.patJetGenJetMatch.maxDPtRel = 1000000.0


        # save a reference to the gen-object instead of a copy.
        # this of course only works if the gen collection is still in the event
        # if we run PAT ourself on GEN-SIM-RECO, it is, so everything is fine
        process.patElectrons.embedGenMatch = False
        process.patMuons.embedGenMatch = False
        process.patPhotons.embedGenMatch = False
        process.patJets.embedGenJetMatch = False
        process.patJets.embedGenPartonMatch = False

        process.patJetCorrFactors.levels = cms.vstring( 'L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual' )
    else:
        process.patJetCorrFactors.levels = cms.vstring( 'L1Offset', 'L2Relative', 'L3Absolute' )

    process.p += process.patDefaultSequence


#adds flavour information for all Gen and Rec-Jets used in skimmer
def addFlavourMatching( process, skimmer, path, runOnGen ):
    for jet_name,jet_def in skimmer.jets.parameters_().items():
        if isinstance( jet_def, cms.PSet ):

            setattr( process,
                     jet_name+'GenJetPartonAssociation',
                     cms.EDProducer( 'JetPartonMatcher',
                                   jets = jet_def.MCLabel,
                                   partons = cms.InputTag( 'patJetPartons' ),
                                   coneSizeToAssociate = cms.double( 0.3 )
                                   )
                     )
            path += getattr( process, jet_name+'GenJetPartonAssociation' )

            setattr( process,
                     jet_name+'GenJetFlavourAlgo',
                     cms.EDProducer( 'JetFlavourIdentifier',
                                   srcByReference = cms.InputTag( jet_name+'GenJetPartonAssociation' ),
                                   physicsDefinition = cms.bool( False )
                                   )
                     )
            path += getattr( process, jet_name+'GenJetFlavourAlgo' )

            setattr( process,
                     jet_name+'GenJetFlavourPhysics',
                     cms.EDProducer( 'JetFlavourIdentifier',
                                   srcByReference = cms.InputTag( jet_name+'GenJetPartonAssociation' ),
                                   physicsDefinition = cms.bool( True )
                                   )
                     )
            path += getattr( process, jet_name+'GenJetFlavourPhysics' )

            if not runOnGen:
               setattr( process,
                        jet_name+'RecoJetPartonAssociation',
                        cms.EDProducer( 'JetPartonMatcher',
                                    jets = jet_def.RecoLabel,
                                    partons = cms.InputTag( 'patJetPartons' ),
                                    coneSizeToAssociate = cms.double( 0.3 )
                                    )
                        )
               path += getattr( process, jet_name+'RecoJetPartonAssociation' )

               setattr( process,
                        jet_name+'RecoJetFlavourPhysics',
                        cms.EDProducer( 'JetFlavourIdentifier',
                                    srcByReference = cms.InputTag( jet_name+'RecoJetPartonAssociation' ),
                                    physicsDefinition = cms.bool( True )
                                    )
                        )
               path += getattr( process, jet_name+'RecoJetFlavourPhysics' )


def configureJEC( process, runOnData ):
    if runOnData:
        jecGlobalTag = cms.string( 'GR_R_42_V19::All' )
    else:
        jecGlobalTag = cms.string( 'START42_V13::All' )
    GlobalTag = process.GlobalTag.globaltag
    if process.GlobalTag.globaltag != jecGlobalTag:
        process.GlobalTag.globaltag = jecGlobalTag
        print "INFO: GlobalTag was '%s' and was changed by configureJEC() to: '%s'" % (GlobalTag, jecGlobalTag)

    process.load( 'JetMETCorrections.Configuration.DefaultJEC_cff' )
    process.load( 'RecoJets.Configuration.RecoPFJets_cff' )
    process.kt6PFJets.doRhoFastjet = True                          # Turn-on the FastJet density calculation
    process.ak5PFJets.doAreaFastjet = True                         # Turn-on the FastJet jet area calculation for ak5PFJets

    process.p += process.kt6PFJets * process.ak5PFJets


def configurePF( process, runOnData, postfix ):
    defaultPostfix = 'PFlow'
    if not postfix:
        postfix = defaultPostfix
        print "WARNING: No postfix provided, setting to: '%s'" %postfix

    from PhysicsTools.PatAlgos.tools import pfTools
    pfTools.usePF2PAT( process,
                       runPF2PAT = True,
                       jetAlgo = 'AK5',
                       jetCorrections = ( 'AK5PFchs', process.patJetCorrFactors.levels ),
                       runOnMC = not runOnData,
                       postfix = postfix
                       )


def configurePFnoPU( process, postfix ):
    defaultPostfix = 'PFlow'
    if not postfix:
        postfix = defaultPostfix
        print "WARNING: No postfix provided, setting to: '%s'" %postfix
    from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

    # Create good primary vertices to be used for PF association
    process.goodOfflinePrimaryVertices = cms.EDFilter(
        'PrimaryVertexObjectFilter',
        filterParams = pvSelector.clone( minNdof = cms.double( 4.0 ), maxZ = cms.double( 24.0 ) ),
        src = cms.InputTag( 'offlinePrimaryVertices' )
        )

    process.pfPileUpPFlow.Enable = True
    process.pfPileUpPFlow.Vertices = 'goodOfflinePrimaryVertices'
    process.pfJetsPFlow.doAreaFastjet = True
    process.pfJetsPFlow.doRhoFastjet = False
    process.patJetCorrFactorsPFlow.rho = cms.InputTag( 'kt6PFJetsPFlow', 'rho' )
    process.pfPileUpPFlow.checkClosestZVertex = cms.bool( False )

    # Compute the mean pt per unit area ("rho") using KT6 Jets with the active areas method.
    from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
    process.kt6PFJetsPFlow = kt4PFJets.clone(
        rParam = cms.double(0.6),
        src = cms.InputTag( 'pfNoElectron' + postfix ),
        doAreaFastjet = cms.bool( True ),
        doRhoFastjet = cms.bool( True )
        )

    getattr( process, 'patPF2PATSequence' + postfix ).replace( getattr( process, 'pfNoElectron' + postfix ), getattr( process, 'pfNoElectron' + postfix ) * process.kt6PFJetsPFlow )
    process.patseq = cms.Sequence(
       process.goodOfflinePrimaryVertices*
       getattr( process, 'patPF2PATSequence' + postfix )
       )
    process.p += process.patseq

    # METnoPU (stolen from: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/lucieg/METsWithPU/METsAnalyzer/python/pfMetNoPileUp_cff.py)
    process.pfMetNoPileUp       = getattr( process, 'pfMET' + postfix ).clone()
    process.pfMetNoPileUp.alias = 'pfMetNoPileUp'
    process.pfMetNoPileUp.src   = cms.InputTag( 'pfNoPileUp' + postfix )

    process.p += process.pfMetNoPileUp

    patMETsPFlowNoPU = 'patMETs' + postfix + 'NoPU'
    setattr( process, patMETsPFlowNoPU, getattr( process, 'patMETs' + postfix ).clone() )
    getattr( process, patMETsPFlowNoPU ).metSource = cms.InputTag( 'pfMetNoPileUp' )

    process.p += getattr( process, patMETsPFlowNoPU )


def addScrapingFilter( process ):
    process.scrapingFilter = cms.EDFilter( 'FilterOutScraping',
                                           applyfilter = cms.untracked.bool( True ),
                                           debugOn = cms.untracked.bool( False ),
                                           numtrack = cms.untracked.uint32( 10 ),
                                           thresh = cms.untracked.double( 0.25 )
                                           )

    process.p_scrapingFilter = cms.Path( process.scrapingFilter )
    process.Skimmer.filterlist.append( 'p_scrapingFilter' )


def addMuGenFilter( process, pt ):
    mugenfilterName = 'mugenfilter' + str( pt )
    # this filter selects events containing muons with pt > pt GeV ...
    mugenfilter = cms.EDFilter( 'MCSmartSingleGenParticleFilter',
                                MaxDecayRadius = cms.untracked.vdouble( 2000.0, 2000.0 ),
                                Status = cms.untracked.vint32( 1, 1 ),
                                MinPt = cms.untracked.vdouble( float( pt ), float( pt ) ),
                                ParticleID = cms.untracked.vint32( 13, -13 ),
                                MaxEta = cms.untracked.vdouble( 2.5, 2.5 ),
                                MinEta = cms.untracked.vdouble( -2.5, -2.5 ),
                                MaxDecayZ = cms.untracked.vdouble( 4000.0, 4000.0 ),
                                MinDecayZ = cms.untracked.vdouble( -4000.0, -4000.0 ),
                                genParSource = cms.InputTag( 'genParticles' )
                                )
    setattr( process, mugenfilterName, mugenfilter.clone() )

    # ... but we don't want these events
    setattr( process, 'p_' + mugenfilterName, cms.Path( ~getattr( process, mugenfilterName ) ) )
    process.Skimmer.filterlist.append( 'p_' + mugenfilterName )


def addEMFilter( process ):
    # this filter selects events containing at least one potential electron candidate ...
    process.emenrichingfilter = cms.EDFilter( 'EMEnrichingFilter',
                                              filterAlgoPSet = cms.PSet( requireTrackMatch = cms.bool( False ),
                                                                         caloIsoMax = cms.double( 10.0 ),
                                                                         isoGenParConeSize = cms.double( 0.1 ),
                                                                         tkIsoMax = cms.double( 5.0 ),
                                                                         hOverEMax = cms.double( 0.5 ),
                                                                         isoGenParETMin = cms.double( 20.0 ),
                                                                         genParSource = cms.InputTag( 'genParticles' ),
                                                                         isoConeSize = cms.double( 0.2 ),
                                                                         clusterThreshold = cms.double( 20.0 )
                                                                         )
                                              )

    # ... but we don't want these events
    process.p_emenrichingfilter = cms.Path( ~process.emenrichingfilter )
    process.Skimmer.filterlist.append( 'p_emenrichingfilter' )


def addBCtoEFilter( process ):
    # this filter selects events containing electrons that come from b or c hadrons ...
    process.bctoefilter = cms.EDFilter( 'BCToEFilter',
                                        filterAlgoPSet = cms.PSet( genParSource = cms.InputTag( 'genParticles' ),
                                                                   eTThreshold = cms.double( 10 )
                                                                   )
                                        )

    # ... but we don't want these events
    process.p_bctoefilter = cms.Path( ~process.bctoefilter )
    process.Skimmer.filterlist.append( 'p_bctoefilter' )


def addBFilter( process ):
    # this filter selects events containing b quarks
    process.bbfilter = cms.EDFilter( 'MCSingleGenParticleFilter',
                                     genParSource = cms.InputTag( 'genParticles' ),
                                     ParticleID = cms.untracked.vint32( 5, -5 ),
                                     Status = cms.untracked.vint32( 2, 2 )
                                     )

    # ... but we don't want these events
    process.p_bbfilter = cms.Path( ~process.bbfilter )
    process.Skimmer.filterlist.append( 'p_bbfilter' )


def addHCALnoiseFilter( process ):
    # Store the result of the HCAL noise info.
    #
    process.HBHENoiseFilterResultProducer = cms.EDProducer(
        'HBHENoiseFilterResultProducer',
        label = cms.InputTag( 'hcalnoise', '', 'RECO' ),
        minRatio = cms.double( -999 ),
        maxRatio = cms.double( 999 ),
        minHPDHits = cms.int32( 17 ),
        minRBXHits = cms.int32( 999 ),
        minHPDNoOtherHits = cms.int32( 10 ),
        minZeros = cms.int32( 10 ),
        minHighEHitTime = cms.double( -9999.0 ),
        maxHighEHitTime = cms.double( 9999.0 ),
        maxRBXEMF = cms.double( -999.0 ),
        minNumIsolatedNoiseChannels = cms.int32( 9999 ),
        minIsolatedNoiseSumE = cms.double( 9999 ),
        minIsolatedNoiseSumEt = cms.double( 9999 ),
        useTS4TS5 = cms.bool( True )
        )

    process.p += process.HBHENoiseFilterResultProducer


def configureTaus( process ):
    # rerun PFTau reco
    process.load( 'RecoTauTag.Configuration.RecoPFTauTag_cff' )
    process.p += process.PFTau
