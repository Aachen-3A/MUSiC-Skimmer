import FWCore.ParameterSet.Config as cms

def configurePAT( process, runOnData, runOnReReco, runOnSummer09 ):
    import PhysicsTools.PatAlgos.tools.jetTools
    #stay consistent wth older samples
    if runOnSummer09:
        PhysicsTools.PatAlgos.tools.jetTools.switchJECSet( process, "Summer09_7TeV" )
    else:
        PhysicsTools.PatAlgos.tools.jetTools.switchJECSet( process, "Summer09_7TeV_ReReco332" )

    if runOnData:
        import PhysicsTools.PatAlgos.tools.coreTools
        PhysicsTools.PatAlgos.tools.coreTools.removeMCMatching( process, ['All'] )
    else:
        if runOnSummer09:
            #in 33x, anti-kt jets are called ak*, however in the 31x they are called antikt*.
            import PhysicsTools.PatAlgos.tools.cmsswVersionTools
            PhysicsTools.PatAlgos.tools.cmsswVersionTools.run33xOn31xMC( process )
        elif runOnReReco:
            #in ReReco of Summer09 there are no ak5GenJets, so add them
            import PhysicsTools.PatAlgos.tools.cmsswVersionTools
            PhysicsTools.PatAlgos.tools.cmsswVersionTools.run33xOnReRecoMC( process )


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
                     cms.EDFilter( 'JetPartonMatcher',
                                   jets = jet_def.MCLabel,
                                   partons = cms.InputTag( 'patJetPartons' ),
                                   coneSizeToAssociate = cms.double( 0.3 )
                                   )
                     )
            path += getattr( process, jet_name+'GenJetPartonAssociation' )

            setattr( process,
                     jet_name+'GenJetFlavourAlgo',
                     cms.EDFilter( 'JetFlavourIdentifier',
                                   srcByReference = cms.InputTag( jet_name+'GenJetPartonAssociation' ),
                                   physicsDefinition = cms.bool( False )
                                   )
                     )
            path += getattr( process, jet_name+'GenJetFlavourAlgo' )

            setattr( process,
                     jet_name+'GenJetFlavourPhysics',
                     cms.EDFilter( 'JetFlavourIdentifier',
                                   srcByReference = cms.InputTag( jet_name+'GenJetPartonAssociation' ),
                                   physicsDefinition = cms.bool( True )
                                   )
                     )
            path += getattr( process, jet_name+'GenJetFlavourPhysics' )

            setattr( process,
                     jet_name+'RecoJetPartonAssociation',
                     cms.EDFilter( 'JetPartonMatcher',
                                   jets = jet_def.RecoLabel,
                                   partons = cms.InputTag( 'patJetPartons' ),
                                   coneSizeToAssociate = cms.double( 0.3 )
                                   )
                     )
            path += getattr( process, jet_name+'RecoJetPartonAssociation' )

            setattr( process,
                     jet_name+'RecoJetFlavourPhysics',
                     cms.EDFilter( 'JetFlavourIdentifier',
                                   srcByReference = cms.InputTag( jet_name+'RecoJetPartonAssociation' ),
                                   physicsDefinition = cms.bool( True )
                                   )
                     )
            path += getattr( process, jet_name+'RecoJetFlavourPhysics' )
