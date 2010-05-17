import FWCore.ParameterSet.Config as cms

def configurePAT( process, runOnData, runOnReReco, runOnSummer09 ):
    if runOnData:
        import PhysicsTools.PatAlgos.tools.coreTools
        PhysicsTools.PatAlgos.tools.coreTools.removeMCMatching( process, ['All'] )
    else:
        import PhysicsTools.PatAlgos.tools.cmsswVersionTools
        if runOnSummer09:
            #in ReReco of Summer09 there are no ak5GenJets, so add them
            PhysicsTools.PatAlgos.tools.cmsswVersionTools.run36xOn35xInput( process, genJets='ak5GenJets' )
        else:
            #change b-tag algos for samples reco'ed with CMSSW < 3.6.X
            PhysicsTools.PatAlgos.tools.cmsswVersionTools.run36xOn35xInput( process )


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



#adds flavour information for all Gen and Rec-Jets used in skimmer
def addFlavourMatching( process, skimmer, path):
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
