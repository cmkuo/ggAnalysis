import FWCore.ParameterSet.Config as cms
from HiggsAnalysis.HiggsTo2photons.hggPhotonIDCuts_cfi import *
#from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#from CMGTools.External.pujetidproducer_cfi import stdalgos, chsalgos

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             hggPhotonIDConfiguration = hggPhotonIDCuts,
                             doGenParticles       = cms.bool(True),
                             runOnParticleGun     = cms.bool(False),
                             runOnSherpa          = cms.bool(False),
                             runL1ECALPrefire     = cms.bool(False),
                             dumpPFPhotons        = cms.bool(True), 
                             dumpJets             = cms.bool(False),
                             dumpAK8Jets          = cms.bool(False),
                             dumpTaus             = cms.bool(False),
                             dumpPDFSystWeight    = cms.bool(False),
                             dumpHFElectrons      = cms.bool(True),
                             development          = cms.bool(False),
                             addFilterInfoMINIAOD = cms.bool(False),
                             doNoHFMET            = cms.bool(False),

                             year                 = cms.int32(2017), 

                             trgFilterDeltaPtCut  = cms.double(0.5),
                             trgFilterDeltaRCut   = cms.double(0.3),

                             triggerEvent         = cms.InputTag("slimmedPatTrigger", "", ""),
                             triggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
                             #patTriggerResults    = cms.InputTag("TriggerResults", "", "PAT"),
                             patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
                             genParticleSrc       = cms.InputTag("prunedGenParticles"),
                             generatorLabel       = cms.InputTag("generator"),
                             LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                             newParticles         = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),
                             pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                             VtxLabel             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             VtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
                             rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                             rhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                             pfMETLabel           = cms.InputTag("slimmedMETs"),
                             electronSrc          = cms.InputTag("slimmedElectrons"),
                             #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
                             calibelectronSrc     = cms.InputTag("slimmedElectrons"),
                             photonSrc            = cms.InputTag("slimmedPhotons"),
                             #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
                             calibphotonSrc       = cms.InputTag("slimmedPhotons"),
                             #muonSrc              = cms.InputTag("slimmedMuons"),
                             muonSrc              = cms.InputTag("cleanedMu"),
                             basicClustersSrc     = cms.InputTag("reducedEgamma", "reducedEBEEClusters"),
                             gsfTrackSrc          = cms.InputTag("reducedEgamma", "reducedGsfTracks"),
                             conversionsSrc       = cms.InputTag("reducedEgamma", "reducedConversions"),
                             beamSpotSrc          = cms.InputTag("offlineBeamSpot"),  
                             ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                             eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                             esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                             recoPhotonSrc             = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
                             TrackLabel                = cms.InputTag("generalTracks"),
                             gsfElectronLabel          = cms.InputTag("gsfElectrons"),
                             PFAllCandidates           = cms.InputTag("particleFlow"),
                             #ak4JetSrc                 = cms.InputTag("updatedJets"),
                             ak4JetSrc                 = cms.InputTag("slimmedJets"),
                             #ak4JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
                             ak8JetSrc                 = cms.InputTag("slimmedJetsAK8"),
                             #ak8JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJECAK8"),
                             #boostedDoubleSVLabel      = cms.InputTag("pfBoostedDoubleSecondaryVertexAK8BJetTags"),
                             tauSrc                    = cms.InputTag("slimmedTaus"),
                             #pfLooseId                 = pfJetIDSelector.clone(),

                             packedPFCands             = cms.InputTag("packedPFCandidates"),
                             elePFClusEcalIsoProducer  = cms.InputTag("electronEcalPFClusterIsolationProducer"),
                             elePFClusHcalIsoProducer  = cms.InputTag("electronHcalPFClusterIsolationProducer"),
                             ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")
)
