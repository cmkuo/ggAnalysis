import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             hggPhotonIDConfiguration = hggPhotonIDCuts,
                             doGenParticles   = cms.bool(True),
                             runOnParticleGun = cms.bool(False),
                             dumpJets         = cms.bool(False),
                             dumpSubJets      = cms.bool(False),
                             dumpTaus         = cms.bool(False),
                             isAOD            = cms.bool(False), #### actually configured through run_data_74x.py
                             runphoIDVID   = cms.bool(False),
                             runeleIDVID   = cms.bool(False),
                             runphoMVAID   = cms.bool(False),
                             runeleMVAID   = cms.bool(False),

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

                             electrons    = cms.InputTag("gedGsfElectrons"),
                             phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-loose"),
                             phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-medium"),
                             phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-PHYS14-PU20bx25-V2-standalone-tight"),
                             eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
                             eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
                             eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
                             eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
                             eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51")
)
