def configurePAT( process, runOnData ):
    if runOnData:
        import PhysicsTools.PatAlgos.tools.coreTools
        PhysicsTools.PatAlgos.tools.coreTools.removeMCMatching( process, 'All' )
    else:
        #in 33x, anti-kt jets are called ak*, however in the 31x they are called antikt*.
        #the following function call will fix this
        import PhysicsTools.PatAlgos.tools.cmsswVersionTools
        PhysicsTools.PatAlgos.tools.cmsswVersionTools.run33xOn31xMC( process )
        
        
        #use 7 TeV JEC
        process.jetCorrFactors.corrSample = 'Summer09_7TeV'

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
        
        process.jetPartonMatch.checkCharge = False
        process.jetPartonMatch.resolveByMatchQuality = True
        process.jetPartonMatch.maxDeltaR = 0.4
        process.jetPartonMatch.maxDPtRel = 1000000.0
        
        process.jetGenJetMatch.checkCharge = False
        process.jetGenJetMatch.resolveByMatchQuality = True
        process.jetGenJetMatch.maxDeltaR = 0.4
        process.jetGenJetMatch.maxDPtRel = 1000000.0
        
        
        # save a reference to the gen-object instead of a copy.
        # this of course only works if the gen collection is still in the event
        # if we run PAT ourself on GEN-SIM-RECO, it is, so everything is fine
        process.allLayer1Electrons.embedGenMatch = False
        process.allLayer1Muons.embedGenMatch = False
        process.allLayer1Photons.embedGenMatch = False
        process.allLayer1Jets.embedGenJetMatch = False
        process.allLayer1Jets.embedGenPartonMatch = False
