import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *


def prepare( runOnGen, runOnData, eleEffAreaTarget,name ,datasetpath ,globalTag , verbosity=0, runOnFast=False ):
    process = cms.Process( 'PATprepare' )

    configureMessenger( process, verbosity )

    import FWCore.Framework.test.cmsExceptionsFatalOption_cff
    process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool( True ),
        # Open file in NOMERGE mode to avoid a memory leak.
        #
        #fileMode = cms.untracked.string( 'NOMERGE' ),
        # Stop processing on each and every thrown exception.
        #
        #Rethrow = FWCore.Framework.test.cmsExceptionsFatalOption_cff.Rethrow
        )
    #process.Tracer = cms.Service("Tracer")
    # The global tag is retrieved automatically but can be changed by the
    # configureJEC function.
    #
    # process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_cff' )

    # Transient track builder is used for the muon vertex refit in the ADD analysis
    process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

    # The global tag is set in pset file or overidden by the calling
    # script (e.g. music_crab3.py9
    process.load( 'Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff' )
    process.GlobalTag.globaltag = globalTag

    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
    process.load('Configuration.StandardSequences.MagneticField_cff')

    # do we need this ?
    #process.content = cms.EDAnalyzer( 'EventContentAnalyzer' )

    # Create an empty path because modules will be added by calling the
    # functions below.

    process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
    process.load( 'PxlSkimmer.Skimming.PxlSkimmer_miniAOD_cfi' )

    process.Skimmer.FastSim = runOnFast

    # Several filters are used while running over data or MC.
    # In order to be flexible, events *not* passing these filtes we do not want
    # to throw these events away but rather write the outcome as a bool into the
    # event.
    # To do this, each filter runs in an own path. These paths are stored in the
    # filterlist. This list is later used to access the value of the filter with
    # help of edm::TriggerResult.
    #
    process.Skimmer.filterlist = cms.vstring()

    # Set effective electron area for corrections.
    #
    #process.Skimmer.EleEffAreaTargetLabel = eleEffAreaTarget

    addElectronIDs( process )
    addGammaIDs( process )
    addNoHFMET( process, runOnData )

    process.Skimmer.FileName = name+'.pxlio'
    process.Skimmer.Process = name
    process.Skimmer.Dataset = datasetpath

    if runOnData:
        process.Skimmer.triggers.HLT.datastreams=cms.vstring(
                                            #"AlCaLumiPixels",
                                            #"AlCaP0",
                                            #"AlCaPhiSym",
                                            "BTagCSV",
                                            "BTagMu",
                                            #"Charmonium",
                                            #"Commissioning",
                                            #"DisplacedJet",
                                            "DoubleEG",
                                            "DoubleMuon",
                                            "DoubleMuonLowMass",
                                            "EcalLaser",
                                            #"EventDisplay",
                                            #"ExpressPhysics",
                                            #"FullTrack",
                                            #"HINCaloJetsOther",
                                            #"HINMuon",
                                            #"HINPFJetsOther",
                                            #"HINPhoton",
                                            #"HLTPhysics",
                                            "HTMHT",
                                            #"HcalHPDNoise",
                                            #"HcalNZS",
                                            #"HighMultiplicity",
                                            "JetHT",
                                            #"L1Accept",
                                            #"LookAreaPD",
                                            "MET",
                                            #"MuOnia",
                                            "MuonEG",
                                            #"NoBPTX",
                                            #"OnlineMonitor",
                                            #"RPCMonitor",
                                            "SingleElectron",
                                            "SingleMuon",
                                            "SinglePhoton",
                                            "Tau",
                                            #"TestEnablesEcalHcal",
                                            #"TestEnablesEcalHcalDQM",
                                            #"ZeroBias",
                                        )
        if not "July17" in datasetpath:
            process.Skimmer.METFilterTag=cms.InputTag("TriggerResults","","RECO")

    if not runOnGen:

        #postfix = ''
        #configurePFJet( process, runOnData, postfix )
        addHCALnoiseFilter( process )
        addHCALLaserEventFilter( process )

    process.Skimmer.filters.AllFilters.paths = process.Skimmer.filterlist
    process.Skimmer.filters.AllFilters.process = process.name_()

    # The skimmer is in the endpath because then the results of all preceding paths
    # are available. This is used to access the outcome of filters that ran.
    #
    process.e = cms.EndPath( process.Skimmer )

    return process

def addGammaIDs( process ):

    # based on https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos/PhotonNtupler/test/runPhotons_VID_CutBased_PHYS14_demo.py
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

    # define which IDs we want to produce
    my_id_modules = [
    'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff',
    'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff',
    ]
    for idmod in my_id_modules:
         setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

    process.p = cms.Path(process.egmPhotonIDSequence )

def addElectronIDs( process ):
    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

    # define which IDs we want to produce
    my_id_modules = [
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                     ]
    #Add them to the VID producer
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

    process.p = cms.Path( process.egmGsfElectronIDSequence )

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

        process.patJetCorrFactors.levels = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' )
    else:
        process.patJetCorrFactors.levels = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute' )

    #process.p += process.patDefaultSequence
    process.patDefaultSequence


#adds flavour information for all Gen and Rec-Jets used in skimmer
def addFlavourMatching( process, skimmer, runOnGen ):
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

# See also:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?rev=116#JetEnCor2012Summer13
# https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC?rev=59
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?rev=449#Winter13_2012_A_B_C_D_datasets_r


#def configureJEC( process, runOnData ):
#    if runOnData:
#        # Newest GT for the ReReco-22Jan2013 data.
#        jecGlobalTag = cms.string( 'FT_53_V21_AN6::All' )
#    else:
#        # Newest GT for CMSSW >= CMSSW_5_3_8_patch3 MC.
#        jecGlobalTag = cms.string( 'START53_V27::All' )
#
#    GlobalTag = process.GlobalTag.globaltag
#
#    process.GlobalTag.globaltag = jecGlobalTag
#    print "INFO: GlobalTag was '%s' and was changed by configureJEC() to: '%s'" % (GlobalTag, jecGlobalTag)


# PF2PAT configuration for jets.
#
def configurePFJet( process, runOnData, postfix ):
    from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
    from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets

    #process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

    process.ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True) # no idea while doArea is false by default, but it's True in RECO so we have to set it
    process.ak5PFJetsCHS = ak5PFJets.clone(src = 'chs', doAreaFastjet = True) # no idea while doArea is false by default, but it's True in RECO so we have to set it
    process.ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')


    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    jetCorrFactors = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute' )
    if runOnData: jetCorrFactors = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' )
    addJetCollection(
       process,
       postfix   = "",
       labelName = 'AK5PFCHS',
       jetSource = cms.InputTag('ak5PFJetsCHS'),
       trackSource = cms.InputTag('unpackedTracksAndVertices'),
       pvSource = cms.InputTag('unpackedTracksAndVertices'),
       jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']),'type-1'),
       btagDiscriminators = [      'combinedSecondaryVertexBJetTags'     ]
       )
    #process.patJetPartonMatchPatJetsAK5PFCHS.matched = "prunedGenParticles"
    process.patJetPartons.src = "prunedGenParticles"
    process.patJetPartons.skipFirstN = cms.uint32(0) # do not skip first 6 particles, we already pruned some!
    process.patJetPartons.acceptNoDaughters = cms.bool(True) # as we drop intermediate stuff, we need to accept quarks with no siblings
    #process.patJetCorrFactorsPatJetsAK5PFCHS.primaryVertices = "offlineSlimmedPrimaryVertices"

    #recreate tracks and pv for btagging
    process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
    process.combinedSecondaryVertex.trackMultiplicityMin = 1 #silly sv, uses un filtered tracks.. i.e. any pt


    #defaultPostfix = 'PFlow'
    #if not postfix:
        #postfix = defaultPostfix
        #print "WARNING: No postfix provided, setting to: '%s'" %postfix

    ## PFnoPU (jets with charged hadron subtraction), see also:
    ## https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?rev=116#JetEnCorPFnoPU2012
    ##
    ## 1. Create good primary vertices to be used for PF association.
    ##
    #from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
    #process.goodOfflinePrimaryVertices = cms.EDFilter(
        #'PrimaryVertexObjectFilter',
        #filterParams = pvSelector.clone( minNdof = cms.double( 4.0 ), maxZ = cms.double( 24.0 ) ),
        #src = cms.InputTag( 'offlinePrimaryVertices' )
        #)

    ## 2. Create the "top-down projection" for the PF2PAT sequence.
    ##
    #from PhysicsTools.PatAlgos.tools import pfTools

    #process.goodOfflinePrimaryVertices

    #jetCorrFactors = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute' )
    #if runOnData: jetCorrFactors = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual' )
    #pfTools.usePF2PAT( process,
                       #runPF2PAT = True,
                       #jetAlgo = 'AK5',
                       ##jetCorrections = ( 'AK5PFchs', jetCorrFactors ),
                       #runOnMC = not runOnData,
                       #postfix = postfix,
                       #pvCollection = cms.InputTag( 'goodOfflinePrimaryVertices' ),
                       #outputModules = []
                       #)
    #process.pfPileUpPFlow.checkClosestZVertex = False

    # 3. Add them all to the sequence.
    #
    #process.patseq = cms.Sequence(
    #   process.goodOfflinePrimaryVertices*
    #   getattr( process, 'patPF2PATSequence' + postfix )
    #   )


# PF2PAT and MET corrections build up on each other, so everything is defined in this one function.

def configurePFMET( process, runOnData ):

   process.load( 'JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff' )

   if runOnData:
      process.corrPfMetType1.jetCorrLabel = cms.string( 'ak5PFL1FastL2L3Residual' )
   else:
      process.corrPfMetType1.jetCorrLabel = cms.string( 'ak5PFL1FastL2L3' )

   process.load( 'JetMETCorrections.Type1MET.correctionTermsPfMetType0PFCandidate_cff' )

   process.load( 'JetMETCorrections.Type1MET.correctionTermsPfMetType0RecoTrack_cff' )

   process.load( 'JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff' )

   if runOnData:
      process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
   else:
      process.corrPfMetShiftXY.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc

   process.load( 'JetMETCorrections.Type1MET.correctedMet_cff' )

   process.correctionTermsPfMetType1Type2
   process.correctionTermsPfMetType0RecoTrack
   process.correctionTermsPfMetType0PFCandidate
   process.correctionTermsPfMetShiftXY

   # pfMET + Type-0(track)
   process.pfMetT0rt
   # PFMET + Type-0(track) + Type-I
   process.pfMetT0rtT1
   # PFMET + Type-0(track) + Type-II
   process.pfMetT0rtT2
   # PFMET + Type-0(track) + xy-shift
   process.pfMetT0rtTxy
   # PFMET + Type-0(track) + Type-I + Type-II
   process.pfMetT0rtT1T2
   # PFMET + Type-0(track) + Type-I + xy-shift
   process.pfMetT0rtT1Txy
   # PFMET + Type-0(track) + Type-II + xy-shift
   process.pfMetT0rtT2Txy
   # PFMET + Type-0(track) + Type-I + Type-II + xy-shift
   process.pfMetT0rtT1T2Txy

   # PFMET + Type-0(pfcand)
   process.pfMetT0pc
   # PFMET + Type-0(pfcand) + Type-I
   process.pfMetT0pcT1
   # PFMET + Type-0(pfcand) + xy-shift
   process.pfMetT0pcTxy
   # PFMET + Type-0(pfcand) + Type-I + xy-shift
   process.pfMetT0pcT1Txy

   # PFMET + Type-I
   process.pfMetT1
   # PFMET + Type-I + Type-II
   process.pfMetT1T2
   # PFMET + Type-I + xy-shift
   process.pfMetT1Txy
   # PFMET + Type-I + Type-II + xy-shift
   process.pfMetT1T2Txy


# Following the "recipe":
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/CommonTools/ParticleFlow/test/pfIsolation_cfg.py?revision=1.2&view=markup
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/CommonTools/ParticleFlow/test/PFIsoReaderDemo.cc?view=markup
#
def configurePFIso( process ):
    from CommonTools.ParticleFlow.Tools import pfIsolation

    process.eleIsoSequence = pfIsolation.setupPFElectronIso( process, 'cleanPatElectrons' )
    process.phoIsoSequence = pfIsolation.setupPFPhotonIso( process, 'cleanPatPhotons' )

    process.pfParticleSelectionSequence
    process.eleIsoSequence
    process.phoIsoSequence


# Median jet pt per area for each event.
# See also:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Rho_for_2011_Effective_Areas
# https://twiki.cern.ch/twiki/bin/view/CMS/Vgamma2011PhotonID#Recommended_cuts
#
def addRhoVariable( process ):
    from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
    process.kt6PFJets50 = kt6PFJets.clone( doRhoFastjet = True )
    process.kt6PFJets50.Rho_EtaMax = cms.double( 5.0 )

    process.kt6PFJets25 = kt6PFJets.clone( doRhoFastjet = True )
    process.kt6PFJets25.Rho_EtaMax = cms.double( 2.5 )

    process.kt6PFJets44 = kt6PFJets.clone( doRhoFastjet = True )

    process.fjSequence = cms.Sequence( process.kt6PFJets25 +
                                       process.kt6PFJets50 +
                                       process.kt6PFJets44
                                       )
    #process.p += process.fjSequence


def addScrapingFilter( process ):
    process.scrapingFilter = cms.EDFilter( 'FilterOutScraping',
                                           applyfilter = cms.untracked.bool( True ),
                                           debugOn = cms.untracked.bool( False ),
                                           numtrack = cms.untracked.uint32( 10 ),
                                           thresh = cms.untracked.double( 0.25 )
                                           )

    process.p_scrapingFilter = cms.Path( process.scrapingFilter )
    process.Skimmer.filterlist.append( 'p_scrapingFilter' )

def addNoHFMET( process, runOnData ):
    import os

    jecUncertaintyFile="PxlSkimmer/Skimming/data/Summer15_25nsV5_DATA_Uncertainty_AK4PFchs.txt"
    process.noHFCands = cms.EDFilter("CandPtrSelector",
                                     src=cms.InputTag("packedPFCandidates"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )

    if True:
        from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
        dBFile =  "Summer15_25nsV5_DATA.db"
        print "If the file "+dBFile+" is not found copy them to your running dir!"
        process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                               connect = cms.string( "sqlite_file:"+dBFile ),
                               toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_Summer15_25nsV5_DATA_AK4PF"),
                label= cms.untracked.string("AK4PF")
                ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_Summer15_25nsV5_DATA_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
                ),
            )
        )
        process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')


    process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
          src = cms.InputTag("slimmedJets"),
          levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
          payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.patJetsReapplyJEC = process.patJetsUpdated.clone(
          jetSource = cms.InputTag("slimmedJets"),
          jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
          )
    process.p += cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )


    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    #default configuration for miniAOD reprocessing, change the isData flag to run on data
    #for a full met computation, remove the pfCandColl input
    runMetCorAndUncFromMiniAOD(process,
                               isData=runOnData,
                               jecUncFile=jecUncertaintyFile,
                               postfix="newUncert"
                               )

    #if not useHFCandidates:
    runMetCorAndUncFromMiniAOD(process,
                               isData=runOnData,
                               pfCandColl=cms.InputTag("noHFCands"),
                               jecUncFile=jecUncertaintyFile,
                               postfix="NoHFhomeMade"
                               )




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
    # (HCAL DPG recommended baseline filter.)
    ##___________________________HCAL_Noise_Filter________________________________||
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
       inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
       reverseDecision = cms.bool(False)
    )

    process.p += process.HBHENoiseFilterResultProducer



def addKinematicsFilter( process ):
    process.load( 'GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi' )

    process.p_kinematicsfilter = cms.Path( process.totalKinematicsFilter )
    process.Skimmer.filterlist.append( 'p_kinematicsfilter' )


def addCSCHaloFilter( process ):
    process.load( 'RecoMET.METAnalyzers.CSCHaloFilter_cfi' )

    process.p_cschalofilter = cms.Path( process.CSCTightHaloFilter )
    process.Skimmer.filterlist.append( 'p_cschalofilter' )


def addHCALLaserEventFilter( process ):
    process.load( 'RecoMET.METFilters.hcalLaserEventFilter_cfi' )
    process.hcalLaserEventFilter.vetoByRunEventNumber = cms.untracked.bool( False )
    process.hcalLaserEventFilter.vetoByHBHEOccupancy  = cms.untracked.bool( True )

    process.p_hcallasereventfilter = cms.Path( process.hcalLaserEventFilter )
    process.Skimmer.filterlist.append( 'p_hcallasereventfilter' )


def addECALDeadCellFilter( process ):
    process.load( 'RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi' )
    process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag( 'ecalTPSkimNA' )

    process.load( 'RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi' )
    process.EcalDeadCellBoundaryEnergyFilter.taggingMode                    = cms.bool( False )
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB      = cms.untracked.double( 10 )
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE      = cms.untracked.double( 10 )
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB            = cms.untracked.double( 100 )
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE            = cms.untracked.double( 100 )
    process.EcalDeadCellBoundaryEnergyFilter.enableGap                      = cms.untracked.bool( False )
    process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32( 12, 14 )
    process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32( 12, 14 )

    # Use BE+TP filter
    process.p_ecaldeadcellfilter = cms.Path( process.EcalDeadCellTriggerPrimitiveFilter * process.EcalDeadCellBoundaryEnergyFilter )
    process.Skimmer.filterlist.append( 'p_ecaldeadcellfilter' )


def addTrackingFailureFilter( process ):
    process.goodVertices = cms.EDFilter(
        'VertexSelector',
        filter = cms.bool( False ),
        src = cms.InputTag( 'offlinePrimaryVertices' ),
        cut = cms.string( '!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2' )
        )

    process.load( 'RecoMET.METFilters.trackingFailureFilter_cfi' )

    process.p_trackingfailurefilter = cms.Path( process.goodVertices * process.trackingFailureFilter )
    process.Skimmer.filterlist.append( 'p_trackingfailurefilter' )


def addEEBadSCFilter( process ):
    process.load('RecoMET.METFilters.eeBadScFilter_cfi')

    process.p_eebadscfilter = cms.Path( process.eeBadScFilter )
    process.Skimmer.filterlist.append( 'p_eebadscfilter' )


def addMuonPFCandidateFilter( process ):
    process.load( 'RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi' )
    process.load( 'RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi' )

    process.p_muonpfcandidatefilter = cms.Path( process.greedyMuonPFCandidateFilter * process.inconsistentMuonPFCandidateFilter )
    process.Skimmer.filterlist.append( 'p_muonpfcandidatefilter' )


# Make sure you check out this user code:
#   cvs co -r Seema11Apr12_52X_V1 -d SandBox/Skims UserCode/seema/SandBox/Skims
# For info & code see this Twiki page:
#   https://twiki.cern.ch/twiki/bin/view/CMS/SusyRA2NJetsInData2011#EB_or_EE_Xtals_with_large_laser
#
def addECALLaserCorrFilter( process ):
    process.load( 'RecoMET.METFilters.ecalLaserCorrFilter_cfi' )

    process.p_ecallasercorrfilter = cms.Path( process.ecalLaserCorrFilter )
    process.Skimmer.filterlist.append( 'p_ecallasercorrfilter' )


def configureTaus( process ):
    # rerun PFTau reco
    process.load( 'RecoTauTag.Configuration.RecoPFTauTag_cff' )
    #process.p += process.PFTau


# Initialize MessageLogger and output report.
#
def configureMessenger( process, verbosity = 0 ):
    process.load( 'FWCore.MessageLogger.MessageLogger_cfi' )
    process.MessageLogger.cerr.threshold = 'INFO'
    process.MessageLogger.cerr.default.limit = -1
    process.MessageLogger.cerr.FwkReport.limit = 100
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000


    process.MessageLogger.categories.append( 'TRIGGERINFO_PXLSKIMMER' )
    process.MessageLogger.categories.append( 'PDFINFO_PXLSKIMMER' )

    if verbosity > 0:
        process.MessageLogger.categories.append( 'EventInfo' )
        process.MessageLogger.categories.append( 'FilterInfo' )
        process.MessageLogger.categories.append( 'TriggerInfo' )
        process.MessageLogger.categories.append( 'PDFInfo' )

    if verbosity > 1:
        process.MessageLogger.categories.append( 'PxlSkimmer' )

    if verbosity > 2:
        process.MessageLogger.cerr.INFO = cms.untracked.PSet( limit = cms.untracked.int32( -1 ) )
