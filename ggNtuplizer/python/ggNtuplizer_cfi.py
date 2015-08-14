import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

#useAOD = True

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             hggPhotonIDConfiguration = hggPhotonIDCuts,
                             doGenParticles   = cms.bool(True),
                             runOnParticleGun = cms.bool(False),
                             dumpPhotons      = cms.bool(True),
                             dumpJets         = cms.bool(False),
                             dumpSubJets      = cms.bool(False),
                             dumpTaus         = cms.bool(False),
                             isAOD            = cms.bool(False), #### actually configured through run_data_74x.py
                             runHFElectrons   = cms.bool(False), #### configured through run_data and run_mc
                             runphoIDVID      = cms.bool(True),
                             runeleIDVID      = cms.bool(True),
                             runphoMVAID      = cms.bool(False),
                             runeleMVAID      = cms.bool(False),
                             development      = cms.bool(False),
                             addFilterInfoAOD     = cms.bool(True),
                             addFilterInfoMINIAOD = cms.bool(False),
                             doNoHFMET            = cms.bool(False),

                             trgFilterDeltaPtCut = cms.double(0.5),
                             trgFilterDeltaRCut  = cms.double(0.3),

                             triggerEvent      = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                             triggerResults    = cms.InputTag("TriggerResults", "", "HLT"),
                             patTriggerResults = cms.InputTag("TriggerResults", "", "PAT"),
                             genParticleSrc    = cms.InputTag("genParticles"),
                             generatorLabel    = cms.InputTag("generator"),
                             newParticles      = cms.vint32(4000011, 4000013, 1000006, 1000022, 1000024, 1000025),
                             pileupCollection  = cms.InputTag("addPileupInfo"),
                             VtxLabel          = cms.InputTag("offlinePrimaryVertices"),
                             VtxBSLabel        = cms.InputTag("offlinePrimaryVerticesWithBS"),
                             rhoLabel          = cms.InputTag("fixedGridRhoFastjetAll"),
                             pfMETLabel        = cms.InputTag("patMETs"),
                             electronSrc       = cms.InputTag("selectedPatElectrons"),
                             hfElectronSrc     = cms.InputTag("hfRecoEcalCandidate"),
                             hfClusterShapeAssociationCollection = cms.InputTag("hfEMClusters"),
                             photonSrc         = cms.InputTag("selectedPatPhotons"),
                             muonSrc           = cms.InputTag("selectedPatMuons"),
                             gsfTrackSrc       = cms.InputTag("electronGsfTracks"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                             esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
                             recoPhotonSrc             = cms.InputTag("gedPhotons"),
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow"),
                             #ak4JetSrc                 = cms.InputTag("selectedPatJets"),
                             ak4JetSrc                 = cms.InputTag("selectedPatJetsAK4PFCHS"),
                             ak8JetSrc                 = cms.InputTag("selectedPatJetsAK8PFCHS"),
                             tauSrc                    = cms.InputTag("selectedPatTaus"),
                             pfLooseId                 = pfJetIDSelector.clone(),
                             
                             phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                             phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
                             phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
                             phoMVAValuesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrigValues"),
                             phoChargedIsolation       = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                             phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                             phoPhotonIsolation        = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                             phoWorstChargedIsolation  = cms.InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"),
                             eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                             eleLooseIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
                             eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
                             eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
                             eleHEEPIdMap    = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
                             eleMVAValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Phys14NonTrigValues")

)
