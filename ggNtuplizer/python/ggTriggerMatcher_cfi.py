import FWCore.ParameterSet.Config as cms


# matches to HLT_IsoMu24_eta2p1_v*
muonTriggerMatchHLTIsoMu24eta2p1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                  src     = cms.InputTag( "cleanPatMuons" ),
                                                  matched = cms.InputTag( "patTrigger" ),
                                                  matchedCuts = cms.string( 'path( "HLT_IsoMu24_eta2p1_v*",1,0)' ),
                                                  maxDPtRel = cms.double( 0.5 ),
                                                  maxDeltaR = cms.double( 0.3 ),
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

# matches to HLT_IsoMu24_v*
muonTriggerMatchHLTIsoMu24 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                            src     = cms.InputTag( "cleanPatMuons" ),
                                            matched = cms.InputTag( "patTrigger" ),
                                            matchedCuts = cms.string( 'path( "HLT_IsoMu24_v*",1,0)' ),
                                            maxDPtRel = cms.double( 0.5 ),
                                            maxDeltaR = cms.double( 0.3 ),
                                            resolveAmbiguities    = cms.bool( True ),
                                            resolveByMatchQuality = cms.bool( True )
                                            )

# matches to HLT_Mu17_Mu8
muonTriggerMatchHLTMu17Mu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                            src     = cms.InputTag( "cleanPatMuons" ),
                                            matched = cms.InputTag( "patTrigger" ),
                                            matchedCuts = cms.string( 'path( "HLT_Mu17_Mu8_v*",1,0)' ),
                                            maxDPtRel = cms.double( 0.5 ),
                                            maxDeltaR = cms.double( 0.3 ),
                                            resolveAmbiguities    = cms.bool( True ),
                                            resolveByMatchQuality = cms.bool( True )
                                            )

# matches to Mu17 of HLT_Mu17_Mu8_v*
muonTriggerMatchHLTMu17forMu17Mu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( "cleanPatMuons" ),
        matched = cms.InputTag( "patTrigger" ),
        matchedCuts = cms.string( 'filter("hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17") || filter("hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17")' ),
        maxDPtRel = cms.double( 0.5 ),
        maxDeltaR = cms.double( 0.3 ),
        resolveAmbiguities    = cms.bool( True ),
        resolveByMatchQuality = cms.bool( True )
        )

# matches to Mu8 of HLT_Mu17_Mu8_v*
muonTriggerMatchHLTMu8forMu17Mu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( "cleanPatMuons" ),
        matched = cms.InputTag( "patTrigger" ),
        matchedCuts = cms.string( 'filter("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8") || filter("hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8")' ),
        maxDPtRel = cms.double( 0.5 ),
        maxDeltaR = cms.double( 0.3 ),
        resolveAmbiguities    = cms.bool( True ),
        resolveByMatchQuality = cms.bool( True )
        )

# matches to Mu17 of HLT_Mu17_TkMu8_v*
muonTriggerMatchHLTMu17forMu17TkMu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( "cleanPatMuons" ),
        matched = cms.InputTag( "patTrigger" ),
        matchedCuts = cms.string( 'filter("hltL3fL1sMu10MuOpenL1f0L2f10L3Filtered17") || filter("hltL3fL1sMu10MuOpenOR3p5L1f0L2f10L3Filtered17")' ),
        maxDPtRel = cms.double( 0.5 ),
        maxDeltaR = cms.double( 0.3 ),
        resolveAmbiguities    = cms.bool( True ),
        resolveByMatchQuality = cms.bool( True )
        )

# matches to Mu8 of HLT_Mu17_TkMu8_v*
muonTriggerMatchHLTMu8forMu17TkMu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
        src     = cms.InputTag( "cleanPatMuons" ),
        matched = cms.InputTag( "patTrigger" ),
        matchedCuts = cms.string( 'filter("hltDiMuonGlb17Trk8DzFiltered0p2") || filter("hltDiMuonGlb17Trk8DzFiltered0p2")' ),
        maxDPtRel = cms.double( 0.5 ),
        maxDeltaR = cms.double( 0.3 ),
        resolveAmbiguities    = cms.bool( True ),
        resolveByMatchQuality = cms.bool( True )
        )

# matches to HLT_Mu17_TkMu8
muonTriggerMatchHLTMu17TkMu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                              src     = cms.InputTag( "cleanPatMuons" ),
                                              matched = cms.InputTag( "patTrigger" ),
                                              matchedCuts = cms.string( 'path( "HLT_Mu17_TkMu8_v*",1,0)' ),
                                              maxDPtRel = cms.double( 0.5 ),
                                              maxDeltaR = cms.double( 0.3 ),
                                              resolveAmbiguities    = cms.bool( True ),
                                              resolveByMatchQuality = cms.bool( True )
                                              )

# matches to HLT_Mu22_TkMu8
muonTriggerMatchHLTMu22TkMu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                              src     = cms.InputTag( "cleanPatMuons" ),
                                              matched = cms.InputTag( "patTrigger" ),
                                              matchedCuts = cms.string( 'path( "HLT_Mu22_TkMu8_v*",1,0)' ),
                                              maxDPtRel = cms.double( 0.5 ),
                                              maxDeltaR = cms.double( 0.3 ),
                                              resolveAmbiguities    = cms.bool( True ),
                                              resolveByMatchQuality = cms.bool( True )
                                              )

# matches to HLT_Mu22_Mu8
muonTriggerMatchHLTMu22Mu8 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                            src     = cms.InputTag( "cleanPatMuons" ),
                                            matched = cms.InputTag( "patTrigger" ),
                                            matchedCuts = cms.string( 'path( "HLT_Mu22_Mu8_v*",1,0)' ),
                                            maxDPtRel = cms.double( 0.5 ),
                                            maxDeltaR = cms.double( 0.3 ),
                                            resolveAmbiguities    = cms.bool( True ),
                                            resolveByMatchQuality = cms.bool( True )
                                            )

# matches to hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter (Single Electron that Electron Et > 17)
electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                                                              #src            = cms.InputTag("cleanPatElectrons"),
                                                                                              src            = cms.InputTag("calibratedPatElectrons"),
                                                                                              matched        = cms.InputTag("patTrigger"),
                                                                                              matchedCuts    = cms.string('filter("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter")'),
                                                                                              maxDPtRel      = cms.double( 0.5 ) ,
                                                                                              maxDeltaR      = cms.double( 0.1 ) ,
                                                                                              resolveAmbiguities    = cms.bool( True ),
                                                                                              resolveByMatchQuality = cms.bool( True )
                                                                                              )

# matches to hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter (Double Electron that Electron Et > 17 for leading and Et > 8 for trailing)
electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter = cms.EDProducer("PATTriggerMatcherDRLessByR" ,
                                                                                                    #src            = cms.InputTag("cleanPatElectrons"),
                                                                                                    src            = cms.InputTag("calibratedPatElectrons"),
                                                                                                    matched        = cms.InputTag( "patTrigger" ),
                                                                                                    matchedCuts    = cms.string( 'filter("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter")'),
                                                                                                    maxDPtRel      = cms.double( 0.5 ) ,
                                                                                                    maxDeltaR      = cms.double( 0.1 ) ,
                                                                                                    resolveAmbiguities    = cms.bool( True ),
                                                                                                    resolveByMatchQuality = cms.bool( True )
                                                                                                    )

# matches to hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ // 8 GeV leg of electron
electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ = cms.EDProducer("PATTriggerMatcherDRLessByR" ,
                                                                                          #src            = cms.InputTag("cleanPatElectrons"),
                                                                                          src            = cms.InputTag("calibratedPatElectrons"),
                                                                                          matched        = cms.InputTag( "patTrigger" ),
                                                                                          matchedCuts    = cms.string( 'filter("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ")'),
                                                                                          maxDPtRel      = cms.double( 0.5 ) ,
                                                                                          maxDeltaR      = cms.double( 0.1 ) ,
                                                                                          resolveAmbiguities    = cms.bool( True ),
                                                                                          resolveByMatchQuality = cms.bool( True )
                                                                                          )

# matches to HLT_Ele27_WP80
electronTriggerMatchHLTEle27WP80 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                  #src         = cms.InputTag( "cleanPatElectrons" ),
                                                  src         = cms.InputTag("calibratedPatElectrons"),
                                                  matched     = cms.InputTag( "patTrigger" ),
                                                  matchedCuts = cms.string( 'path( "HLT_Ele27_WP80_v*",1,0)' ),
                                                  maxDPtRel   = cms.double( 0.5 ),
                                                  maxDeltaR   = cms.double( 0.3 ),
                                                  resolveAmbiguities    = cms.bool( True ),
                                                  resolveByMatchQuality = cms.bool( True )
                                                  )

# matches to HLT_PAPhoton15_TightCaloIdVL for pA analysis
electronTriggerMatchHLTPAPhoton15TightCaloIdVL = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                                src         = cms.InputTag("calibratedPatElectrons"),
                                                                matched     = cms.InputTag( "patTrigger" ),
                                                                matchedCuts = cms.string( 'path( "HLT_PAPhoton15_TightCaloIdVL_v*",1,0)' ),
                                                                maxDPtRel   = cms.double( 0.5 ),
                                                                maxDeltaR   = cms.double( 0.3 ),
                                                                resolveAmbiguities    = cms.bool( True ),
                                                                resolveByMatchQuality = cms.bool( True )
                                                                )

# matches to HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v
electronTriggerMatchEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                                                              #src     = cms.InputTag( "cleanPatElectrons" ),
                                                                                                              src     = cms.InputTag("calibratedPatElectrons"),
                                                                                                              matched = cms.InputTag( "patTrigger" ),
                                                                                                              matchedCuts = cms.string( 'path( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",1,0)' ),
                                                                                                              maxDPtRel = cms.double( 0.5 ),
                                                                                                              maxDeltaR = cms.double( 0.3 ),
                                                                                                              resolveAmbiguities    = cms.bool( True ),
                                                                                                              resolveByMatchQuality = cms.bool( True )
                                                                                                              )

# matches to HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
electronTriggerMatchEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                                                              #src     = cms.InputTag( "cleanPatElectrons" ),
                                                                                                              src     = cms.InputTag("calibratedPatElectrons"),
                                                                                                              matched = cms.InputTag( "patTrigger" ),
                                                                                                              matchedCuts = cms.string( 'path( "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_**",1,0)' ),
                                                                                                              maxDPtRel = cms.double( 0.5 ),
                                                                                                              maxDeltaR = cms.double( 0.3 ),
                                                                                                              resolveAmbiguities    = cms.bool( True ),
                                                                                                              resolveByMatchQuality = cms.bool( True )
                                                                                                              )

# matches to HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*
electronTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                                       #src     = cms.InputTag( "cleanPatElectrons" ),
                                                                                       src     = cms.InputTag("calibratedPatElectrons"),
                                                                                       matched = cms.InputTag( "patTrigger" ),
                                                                                       matchedCuts = cms.string( 'path( "HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*",1,0)' ),
                                                                                       maxDPtRel = cms.double( 0.5 ),
                                                                                       maxDeltaR = cms.double( 0.1 ),
                                                                                       resolveAmbiguities    = cms.bool( True ),
                                                                                       resolveByMatchQuality = cms.bool( True )
                                                                                       )



# matches to Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT leg for electron for electron
electronFilterMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                          #src     = cms.InputTag( "cleanPatElectrons" ),
                                                                          src     = cms.InputTag("calibratedPatElectrons"),
                                                                          matched = cms.InputTag( "patTrigger" ),
                                                                          matchedCuts = cms.string( 'filter("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter") || filter("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter")' ),
                                                                          maxDPtRel = cms.double( 0.5 ),
                                                                          maxDeltaR = cms.double( 0.1 ),
                                                                          resolveAmbiguities    = cms.bool( True ),
                                                                          resolveByMatchQuality = cms.bool( True )
                                                                          )

# matches to Ele8_Mass50 leg for electron
electronFilterMatchEle8Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                #src     = cms.InputTag( "cleanPatElectrons" ),
                                                src     = cms.InputTag("calibratedPatElectrons"),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.1 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                                )

# matches to Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT leg for electron
electronFilterMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                          #src     = cms.InputTag( "cleanPatElectrons" ),
                                                                          src     = cms.InputTag("calibratedPatElectrons"),
                                                                          matched = cms.InputTag( "patTrigger" ),
                                                                          matchedCuts = cms.string( 'filter("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsolFilter") || filter("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter")' ),
                                                                          maxDPtRel = cms.double( 0.5 ),
                                                                          maxDeltaR = cms.double( 0.1 ),
                                                                          resolveAmbiguities    = cms.bool( True ),
                                                                          resolveByMatchQuality = cms.bool( True )
                                                                          )

# matches to SC4_Mass50 leg for electron
electronFilterMatchSC4Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                               #src     = cms.InputTag( "cleanPatElectrons" ),
                                               src     = cms.InputTag("calibratedPatElectrons"),
                                               matched = cms.InputTag( "patTrigger" ),
                                               matchedCuts = cms.string( 'filter("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter")' ),
                                               maxDPtRel = cms.double( 0.5 ),
                                               maxDeltaR = cms.double( 0.1 ),
                                               resolveAmbiguities    = cms.bool( True ),
                                               resolveByMatchQuality = cms.bool( True )
                                               )

# matches to Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT leg for electron
electronFilterMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                       #src     = cms.InputTag( "cleanPatElectrons" ),
                                                                       src     = cms.InputTag("calibratedPatElectrons"),
                                                                       matched = cms.InputTag( "patTrigger" ),
                                                                       matchedCuts = cms.string( 'filter("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter") || filter("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter")' ),
                                                                       maxDPtRel = cms.double( 0.5 ),
                                                                       maxDeltaR = cms.double( 0.1 ),
                                                                       resolveAmbiguities    = cms.bool( True ),
                                                                       resolveByMatchQuality = cms.bool( True )
                                                                       )

# matches to SC17_Mass50 leg for electron
electronFilterMatchSC17Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                #src     = cms.InputTag( "cleanPatElectrons" ),
                                                src     = cms.InputTag("calibratedPatElectrons"),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.1 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                                )

# check point 2
# matches to Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT leg
eleTrgFMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                   src     = cms.InputTag( "cleanPatPhotons" ),
                                                                   matched = cms.InputTag( "patTrigger" ),
                                                                   matchedCuts = cms.string( 'filter("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter") || filter("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter")' ),
                                                                   maxDPtRel = cms.double( 0.5 ),
                                                                   maxDeltaR = cms.double( 0.3 ),
                                                                   resolveAmbiguities    = cms.bool( True ),
                                                                   resolveByMatchQuality = cms.bool( True )
                                                                   )

# matches to Ele8_Mass50 leg
eleTrgFMatchEle8Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                         src     = cms.InputTag( "cleanPatPhotons" ),
                                         matched = cms.InputTag( "patTrigger" ),
                                         matchedCuts = cms.string( 'filter("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter")' ),
                                         maxDPtRel = cms.double( 0.5 ),
                                         maxDeltaR = cms.double( 0.3 ),
                                         resolveAmbiguities    = cms.bool( True ),
                                         resolveByMatchQuality = cms.bool( True )
                                         )

# matches to Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT leg
eleTrgFMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                   src     = cms.InputTag( "cleanPatPhotons" ),
                                                                   matched = cms.InputTag( "patTrigger" ),
                                                                   matchedCuts = cms.string( 'filter("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsolFilter") || filter("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter")' ),
                                                                   maxDPtRel = cms.double( 0.5 ),
                                                                   maxDeltaR = cms.double( 0.3 ),
                                                                   resolveAmbiguities    = cms.bool( True ),
                                                                   resolveByMatchQuality = cms.bool( True )
                                                                   )

# matches to SC4_Mass50 leg
eleTrgFMatchSC4Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                        src     = cms.InputTag( "cleanPatPhotons" ),
                                        matched = cms.InputTag( "patTrigger" ),
                                        matchedCuts = cms.string( 'filter("hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter")' ),
                                        maxDPtRel = cms.double( 0.5 ),
                                        maxDeltaR = cms.double( 0.3 ),
                                        resolveAmbiguities    = cms.bool( True ),
                                        resolveByMatchQuality = cms.bool( True )
                                        )

# matches to Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT leg
eleTrgFMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                              src     = cms.InputTag( "cleanPatPhotons" ),
                                              matched = cms.InputTag( "patTrigger" ),
                                              matchedCuts = cms.string( 'filter("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter") || filter("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter")' ),
                                              maxDPtRel = cms.double( 0.5 ),
                                              maxDeltaR = cms.double( 0.3 ),
                                              resolveAmbiguities    = cms.bool( True ),
                                              resolveByMatchQuality = cms.bool( True )
                                              )

# matches to SC17_Mass50 leg
eleTrgFMatchSC17Mass50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                         src     = cms.InputTag( "cleanPatPhotons" ),
                                         matched = cms.InputTag( "patTrigger" ),
                                         matchedCuts = cms.string( 'filter("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter")' ),
                                         maxDPtRel = cms.double( 0.5 ),
                                         maxDeltaR = cms.double( 0.3 ),
                                         resolveAmbiguities    = cms.bool( True ),
                                         resolveByMatchQuality = cms.bool( True )
                                         )

# Iso-Iso 26_18
phoTrgFMatch26IdIso18IdIso18 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18CaloIdLIsoVLTrackIsolDoubleLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26IdIso18IdIso26 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                               matchedCuts = cms.string( 'filter("hltEG26CaloIdLIsoVLHcalIsolLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )
phoTrgFMatch26IdIsoXL18IdIsoXL18 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18CaloIdXLIsoXLTrackIsolDoubleLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26IdIsoXL18IdIsoXL26 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                               matchedCuts = cms.string( 'filter("hltEG26CaloIdXLIsoXLHcalIsolLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )


# R9-R9 26_18
phoTrgFMatch26R918R918 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18R9IdDoubleLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26R918R926 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG26R9IdLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )
phoTrgFMatch26R9T18R9T18 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18R9IdTDoubleLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26R9T18R9T26 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG26R9IdTLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )


# R9-Iso 26_18
phoTrgFMatch26R918IdIso18 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG26R9IdLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26R918IdIso26 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18CaloIdLIsoVLTrackIsolLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26R918IdIsoXL18 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG26R9IdTLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch26R918IdIsoXL26 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18CaloIdXLIsoXLTrackIsolLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

# Iso-R9 26_18
phoTrgFMatch26IdIso18R918 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18R9IdLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                            )

phoTrgFMatch26IdIso18R926 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG26CaloIdLIsoVLTrackIsolLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                            )

phoTrgFMatch26IdIsoXL18R918 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18R9IdTLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                            )

phoTrgFMatch26IdIsoXL18R926 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG26CaloIdXLIsoXLTrackIsolLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                            )


# Iso-Iso 36_22
phoTrgFMatch36IdIso22IdIso22 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG22CaloIdLIsoVLTrackIsolDoubleLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch36IdIso22IdIso36 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG36CaloIdLIsoVLHcalIsolLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

# R9-R9 36_22
phoTrgFMatch36R922R922 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG22R9IdDoubleLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch36R922R936 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG36R9IdLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

# R9-Iso 36_22
phoTrgFMatch36R922IdIso22 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG36R9IdLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

phoTrgFMatch36R922IdIso36 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG22CaloIdLIsoVLTrackIsolLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                               )

# Iso-R9 36_22
phoTrgFMatch36IdIso22R922 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG22R9IdLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                            )

phoTrgFMatch36IdIso22R936 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG36CaloIdLIsoVLTrackIsolLastFilter")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                            )

# EB Only photon 30 (parking data)
phoTrgFMatch30EBOnly = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                       src     = cms.InputTag( "cleanPatPhotons" ),
                                       matched = cms.InputTag( "patTrigger" ),
                                       matchedCuts = cms.string( 'filter("hltPhoton30R9Id90CaloIdHE10Iso40EBOnlyTrackIsoLastFilter")' ),
                                       maxDPtRel = cms.double( 0.5 ),
                                       maxDeltaR = cms.double( 0.3 ),
                                       resolveAmbiguities    = cms.bool( True ),
                                       resolveByMatchQuality = cms.bool( True )
                                       )

# matches to HLT_Photon26_Photon18
photonTriggerMatchHLTPhoton26Photon18= cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                        src     = cms.InputTag( "cleanPatPhotons" ),
                                                        matched = cms.InputTag( "patTrigger" ),
                                                        matchedCuts = cms.string( 'path( "HLT_Photon26_Photon18_v*",1,0)' ),
                                                        maxDPtRel = cms.double( 0.5 ),
                                                        maxDeltaR = cms.double( 0.3 ),
                                                        resolveAmbiguities    = cms.bool( True ),
                                                        resolveByMatchQuality = cms.bool( True )
                                                        )

# matches to HLT_Photon36_Photon22
photonTriggerMatchHLTPhoton36Photon22 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                        src     = cms.InputTag( "cleanPatPhotons" ),
                                                        matched = cms.InputTag( "patTrigger" ),
                                                        matchedCuts = cms.string( 'path( "HLT_Photon36_Photon22_v*",1,0)' ),
                                                        maxDPtRel = cms.double( 0.5 ),
                                                        maxDeltaR = cms.double( 0.3 ),
                                                        resolveAmbiguities    = cms.bool( True ),
                                                        resolveByMatchQuality = cms.bool( True )
                                                        )

# matches to HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60
photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass60 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                                                        src     = cms.InputTag( "cleanPatPhotons" ),
                                                                                                        matched = cms.InputTag( "patTrigger" ),
                                                                                                        matchedCuts = cms.string( 'path( "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",1,0)' ),
                                                                                                        maxDPtRel = cms.double( 0.5 ),
                                                                                                        maxDeltaR = cms.double( 0.3 ),
                                                                                                        resolveAmbiguities    = cms.bool( True ),
                                                                                                        resolveByMatchQuality = cms.bool( True )
                                                                                                        )

# matches to HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70
photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass70 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                                                        src     = cms.InputTag( "cleanPatPhotons" ),
                                                                                                        matched = cms.InputTag( "patTrigger" ),
                                                                                                        matchedCuts = cms.string( 'path( "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*",1,0)' ),
                                                                                                        maxDPtRel = cms.double( 0.5 ),
                                                                                                        maxDeltaR = cms.double( 0.3 ),
                                                                                                        resolveAmbiguities    = cms.bool( True ),
                                                                                                        resolveByMatchQuality = cms.bool( True )
                                                                                                        )

# matches to HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50
photonTriggerMatchHLTPhoton36R9Id85ORCaloId10Iso50Photon22R9Id85ORCaloId10Iso50 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                                                                  src     = cms.InputTag( "cleanPatPhotons" ),
                                                                                                  matched = cms.InputTag( "patTrigger" ),
                                                                                                  matchedCuts = cms.string( 'path( "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*",1,0)' ),
                                                                                                  maxDPtRel = cms.double( 0.5 ),
                                                                                                  maxDeltaR = cms.double( 0.3 ),
                                                                                                  resolveAmbiguities    = cms.bool( True ),
                                                                                                  resolveByMatchQuality = cms.bool( True )
                                                                                                  )

# 8TeV diphoton trigger filters
phoTrgFMatch26Id10Iso50HcalIso = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                 src     = cms.InputTag( "cleanPatPhotons" ),
                                                 matched = cms.InputTag( "patTrigger" ),
                                                 matchedCuts = cms.string( 'filter("hltEG26CaloId10Iso50HcalIsoLastFilter")' ),
                                                 maxDPtRel = cms.double( 0.5 ),
                                                 maxDeltaR = cms.double( 0.3 ),
                                                 resolveAmbiguities    = cms.bool( True ),
                                                 resolveByMatchQuality = cms.bool( True )
                                                 )

phoTrgFMatch26Id10Isol50HcalIsol = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                   src     = cms.InputTag( "cleanPatPhotons" ),
                                                   matched = cms.InputTag( "patTrigger" ),
                                                   matchedCuts = cms.string( 'filter("hltEG26CaloId10Iso50HcalIsolLastFilter")' ),
                                                   maxDPtRel = cms.double( 0.5 ),
                                                   maxDeltaR = cms.double( 0.3 ),
                                                   resolveAmbiguities    = cms.bool( True ),
                                                   resolveByMatchQuality = cms.bool( True )
                                                   )

phoTrgFMatch26R9Id85 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                       src     = cms.InputTag( "cleanPatPhotons" ),
                                       matched = cms.InputTag( "patTrigger" ),
                                       matchedCuts = cms.string( 'filter("hltEG26R9Id85LastFilter")' ),
                                       maxDPtRel = cms.double( 0.5 ),
                                       maxDeltaR = cms.double( 0.3 ),
                                       resolveAmbiguities    = cms.bool( True ),
                                       resolveByMatchQuality = cms.bool( True )
                                       )

phoTrgFMatch18R9Id85 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                       src     = cms.InputTag( "cleanPatPhotons" ),
                                       matched = cms.InputTag( "patTrigger" ),
                                       matchedCuts = cms.string( 'filter("hltEG18R9Id85LastFilterUnseeded")' ),
                                       maxDPtRel = cms.double( 0.5 ),
                                       maxDeltaR = cms.double( 0.3 ),
                                       resolveAmbiguities    = cms.bool( True ),
                                       resolveByMatchQuality = cms.bool( True )
                                       )

phoTrgFMatch18Id10Iso50TrkIso = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                src     = cms.InputTag( "cleanPatPhotons" ),
                                                matched = cms.InputTag( "patTrigger" ),
                                                matchedCuts = cms.string( 'filter("hltEG18CaloId10Iso50TrackIsoLastFilterUnseeded")' ),
                                                maxDPtRel = cms.double( 0.5 ),
                                                maxDeltaR = cms.double( 0.3 ),
                                                resolveAmbiguities    = cms.bool( True ),
                                                resolveByMatchQuality = cms.bool( True )
                                                )

phoTrgFMatch18Id10Iso50TrkIsol = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                 src     = cms.InputTag( "cleanPatPhotons" ),
                                                 matched = cms.InputTag( "patTrigger" ),
                                                 matchedCuts = cms.string( 'filter("hltEG18CaloId10Iso50TrackIsolLastFilterUnseeded")' ),
                                                 maxDPtRel = cms.double( 0.5 ),
                                                 maxDeltaR = cms.double( 0.3 ),
                                                 resolveAmbiguities    = cms.bool( True ),
                                                 resolveByMatchQuality = cms.bool( True )
                                                 )

phoTrgFMatch18Id10Iso50TrkIsoDouble = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                      src     = cms.InputTag( "cleanPatPhotons" ),
                                                      matched = cms.InputTag( "patTrigger" ),
                                                      matchedCuts = cms.string( 'filter("hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded")' ),
                                                      maxDPtRel = cms.double( 0.5 ),
                                                      maxDeltaR = cms.double( 0.3 ),
                                                      resolveAmbiguities    = cms.bool( True ),
                                                      resolveByMatchQuality = cms.bool( True )
                                                      )

phoTrgFMatch18Id10Iso50TrkIsolDouble = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                       src     = cms.InputTag( "cleanPatPhotons" ),
                                                       matched = cms.InputTag( "patTrigger" ),
                                                       matchedCuts = cms.string( 'filter("hltEG18CaloId10Iso50TrackIsolDoubleLastFilterUnseeded")' ),
                                                       maxDPtRel = cms.double( 0.5 ),
                                                       maxDeltaR = cms.double( 0.3 ),
                                                       resolveAmbiguities    = cms.bool( True ),
                                                       resolveByMatchQuality = cms.bool( True )
                                                       )

eleTrgFMatch27WP80TrkIso = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "cleanPatPhotons" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'filter("hltEle27WP80TrackIsoFilter")' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

# matches to HLT_Jet_30
jetTriggerMatchHLTJet30 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                          src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                          matched = cms.InputTag( "patTrigger" ),
                                          matchedCuts = cms.string( 'path( "HLT_Jet30_v*",1,0)' ),
                                          maxDPtRel = cms.double( 0.5 ),
                                          maxDeltaR = cms.double( 0.3 ),
                                          resolveAmbiguities    = cms.bool( True ),
                                          resolveByMatchQuality = cms.bool( True )
                                          )

# matches to HLT_Jet_60
jetTriggerMatchHLTJet60 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                          src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                          matched = cms.InputTag( "patTrigger" ),
                                          matchedCuts = cms.string( 'path( "HLT_Jet60_v*",1,0)' ),
                                          maxDPtRel = cms.double( 0.5 ),
                                          maxDeltaR = cms.double( 0.3 ),
                                          resolveAmbiguities    = cms.bool( True ),
                                          resolveByMatchQuality = cms.bool( True )
                                          )

# matches to HLT_Jet_80
jetTriggerMatchHLTJet80 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                          src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                          matched = cms.InputTag( "patTrigger" ),
                                          matchedCuts = cms.string( 'path( "HLT_Jet80_v*",1,0)' ),
                                          maxDPtRel = cms.double( 0.5 ),
                                          maxDeltaR = cms.double( 0.3 ),
                                          resolveAmbiguities    = cms.bool( True ),
                                          resolveByMatchQuality = cms.bool( True )
                                          )

# matches to HLT_Jet_110
jetTriggerMatchHLTJet110 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'path( "HLT_Jet110_v*",1,0)' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

# matches to HLT_Jet_150
jetTriggerMatchHLTJet150 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'path( "HLT_Jet150_v*",1,0)' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

# matches to HLT_Jet_190
jetTriggerMatchHLTJet190 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'path( "HLT_Jet190_v*",1,0)' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

# matches to HLT_Jet_240
jetTriggerMatchHLTJet240 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'path( "HLT_Jet240_v*",1,0)' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

# matches to HLT_Jet_300
jetTriggerMatchHLTJet300 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'path( "HLT_Jet300_v*",1,0)' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

# matches to HLT_Jet_370
jetTriggerMatchHLTJet370 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                           src     = cms.InputTag( "selectedPatJetsAK5PF" ),
                                           matched = cms.InputTag( "patTrigger" ),
                                           matchedCuts = cms.string( 'path( "HLT_Jet370_v*",1,0)' ),
                                           maxDPtRel = cms.double( 0.5 ),
                                           maxDeltaR = cms.double( 0.3 ),
                                           resolveAmbiguities    = cms.bool( True ),
                                           resolveByMatchQuality = cms.bool( True )
                                           )

ggTriggerMatcherElectron = cms.Sequence(
    electronTriggerMatchHLTEle27WP80 +
    electronTriggerMatchHLTPAPhoton15TightCaloIdVL +
    electronTriggerMatchEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL +
    electronTriggerMatchEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL +
    electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter +
    electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter + 
    electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ +
    electronTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4Mass50 +
    electronFilterMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT +
    electronFilterMatchEle8Mass50 +
    electronFilterMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT +
    electronFilterMatchSC4Mass50 +
    electronFilterMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT +
    electronFilterMatchSC17Mass50
    )

ggTriggerMatcherPhoton = cms.Sequence(
    photonTriggerMatchHLTPhoton26Photon18 +
    photonTriggerMatchHLTPhoton36Photon22 +
    photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass60 +
    photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass70 +
    photonTriggerMatchHLTPhoton36R9Id85ORCaloId10Iso50Photon22R9Id85ORCaloId10Iso50 +
    
    eleTrgFMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT +
    eleTrgFMatchEle8Mass50 +
    eleTrgFMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT +
    eleTrgFMatchSC4Mass50 +
    eleTrgFMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT +
    eleTrgFMatchSC17Mass50 +
    
    eleTrgFMatch27WP80TrkIso +

    phoTrgFMatch26IdIso18IdIso18 +
    phoTrgFMatch26IdIso18IdIso26 +
    phoTrgFMatch26R918R918 +
    phoTrgFMatch26R918R926 +
    phoTrgFMatch26R918IdIso18 +
    phoTrgFMatch26R918IdIso26 +
    phoTrgFMatch26IdIso18R918 +
    phoTrgFMatch26IdIso18R926 +

    phoTrgFMatch36IdIso22IdIso22 +
    phoTrgFMatch36IdIso22IdIso36 +
    phoTrgFMatch36R922R922 +
    phoTrgFMatch36R922R936 +
    phoTrgFMatch36R922IdIso22 +
    phoTrgFMatch36R922IdIso36 +
    phoTrgFMatch36IdIso22R922 +
    phoTrgFMatch36IdIso22R936 +

    phoTrgFMatch26IdIsoXL18IdIsoXL18 +
    phoTrgFMatch26IdIsoXL18IdIsoXL26 +
    phoTrgFMatch26R9T18R9T18 +
    phoTrgFMatch26R9T18R9T26 +
    phoTrgFMatch26R918IdIsoXL18 +
    phoTrgFMatch26R918IdIsoXL26 +
    phoTrgFMatch26IdIsoXL18R918 +
    phoTrgFMatch26IdIsoXL18R926 +

    phoTrgFMatch26Id10Iso50HcalIso +
    phoTrgFMatch26R9Id85 +
    phoTrgFMatch18R9Id85 +
    phoTrgFMatch18Id10Iso50TrkIso + 
    phoTrgFMatch18Id10Iso50TrkIsoDouble +
    phoTrgFMatch26Id10Isol50HcalIsol +
    phoTrgFMatch18Id10Iso50TrkIsol +
    phoTrgFMatch18Id10Iso50TrkIsolDouble +

    phoTrgFMatch30EBOnly
    )

ggTriggerMatcherMuon = cms.Sequence(
    muonTriggerMatchHLTIsoMu24eta2p1 +
    muonTriggerMatchHLTIsoMu24 +
    muonTriggerMatchHLTMu17Mu8 +
    muonTriggerMatchHLTMu17TkMu8 +
    muonTriggerMatchHLTMu22TkMu8 +
    muonTriggerMatchHLTMu22Mu8 +
    muonTriggerMatchHLTMu17forMu17Mu8 +
    muonTriggerMatchHLTMu8forMu17Mu8 +
    muonTriggerMatchHLTMu17forMu17TkMu8 +
    muonTriggerMatchHLTMu8forMu17TkMu8

    )

ggTriggerMatcherJet = cms.Sequence(
    jetTriggerMatchHLTJet30 +
    jetTriggerMatchHLTJet60 +
    jetTriggerMatchHLTJet80 +
    jetTriggerMatchHLTJet110 +
    jetTriggerMatchHLTJet150 +
    jetTriggerMatchHLTJet190 +
    jetTriggerMatchHLTJet240 +
    jetTriggerMatchHLTJet300 +
    jetTriggerMatchHLTJet370 
    )

ggTriggerMatcher = cms.Sequence(
    ggTriggerMatcherElectron + ggTriggerMatcherMuon + ggTriggerMatcherPhoton + ggTriggerMatcherJet
    )

ggTriggerMatcherNoJet = cms.Sequence(
    ggTriggerMatcherElectron + ggTriggerMatcherMuon + ggTriggerMatcherPhoton
    )
