import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             hggPhotonIDConfiguration = hggPhotonIDCuts,
                             doGenParticles   = cms.bool(True),
                             runOnParticleGun = cms.bool(False),
                             dumpPhotons      = cms.bool(True), 
                             dumpJets         = cms.bool(False),
                             dumpSubJets      = cms.bool(False),
                             dumpTaus         = cms.bool(False),
                             isAOD            = cms.bool(False), #### actually configured through run_data_74x.py
                             runphoIDVID      = cms.bool(True),
                             runeleIDVID      = cms.bool(True),
                             runphoMVAID      = cms.bool(False),
                             runeleMVAID      = cms.bool(False),

                             trgFilterDeltaPtCut  = cms.double(0.05),
                             trgFilterDeltaEtaCut = cms.double(0.05),

                             triggerEvent     = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                             triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
                             genParticleSrc   = cms.InputTag("prunedGenParticles"),
                             generatorLabel   = cms.InputTag("generator"),
                             newParticles     = cms.vint32(4000011, 4000013, 1000006, 1000022, 1000024, 1000025),
                             pileupCollection = cms.InputTag("addPileupInfo"),
                             VtxLabel         = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             VtxBSLabel       = cms.InputTag("offlinePrimaryVerticesWithBS"),
                             rhoLabel         = cms.InputTag("fixedGridRhoFastjetAll"),
                             pfMETLabel       = cms.InputTag("slimmedMETs"),
                             electronSrc      = cms.InputTag("slimmedElectrons"),
                             photonSrc        = cms.InputTag("slimmedPhotons"),
                             muonSrc          = cms.InputTag("slimmedMuons"),
                             ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                             esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                             recoPhotonSrc             = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow"),
                             ak4JetSrc                 = cms.InputTag("slimmedJets"),
                             ak8JetSrc                 = cms.InputTag("slimmedJetsAK8"),
                             tauSrc                    = cms.InputTag("slimmedTaus"),
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
