import FWCore.ParameterSet.Config as cms

from CommonTools.PileupAlgos.Puppi_cff import *

puppi.candName = cms.InputTag('packedPFCandidates')
puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
puppi.puppiForLeptons = cms.bool(False)

egmPhotonIsolationMiniAOD = cms.EDProducer( "CITKPFIsolationSumProducer",
                                            srcToIsolate = cms.InputTag("slimmedPhotons"),
                                            srcForIsolationCone = cms.InputTag('packedPFCandidates'),
                                            isolationConeDefinitions = cms.VPSet(
        cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
                  coneSize = cms.double(0.3),
                  isolateAgainst = cms.string('h+'),
                  miniAODVertexCodes = cms.vuint32(2,3),
                  vertexIndex = cms.int32(0),
                  ),
        cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
                  coneSize = cms.double(0.3),
                  isolateAgainst = cms.string('h0'),
                  miniAODVertexCodes = cms.vuint32(2,3),
                  vertexIndex = cms.int32(0),
                  ),
        cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
                  coneSize = cms.double(0.3),
                  isolateAgainst = cms.string('gamma'),
                  miniAODVertexCodes = cms.vuint32(2,3),
                  vertexIndex = cms.int32(0),
                  )
        )
                                            )

egmPhotonIsolationMiniAODPUPPI = cms.EDProducer( "CITKPFIsolationSumProducerForPUPPI",
                                                 srcToIsolate = cms.InputTag("slimmedPhotons"),
                                                 srcForIsolationCone = cms.InputTag('packedPFCandidates'),
                                                 puppiValueMap = cms.InputTag('puppi'),
                                                 isolationConeDefinitions = cms.VPSet(
        cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
                  coneSize = cms.double(0.3),
                  isolateAgainst = cms.string('h+'),
                  miniAODVertexCodes = cms.vuint32(2,3),
                  vertexIndex = cms.int32(0),
                  ),
        cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
                  coneSize = cms.double(0.3),
                  isolateAgainst = cms.string('h0'),
                  miniAODVertexCodes = cms.vuint32(2,3),
                  vertexIndex = cms.int32(0),
                  ),
        cms.PSet( isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'), 
                  coneSize = cms.double(0.3),
                  isolateAgainst = cms.string('gamma'),
                  miniAODVertexCodes = cms.vuint32(2,3),
                  vertexIndex = cms.int32(0),
                  )
        )
                                                 )

ggPhotonIso_CITK_PUPPI_Sequence = cms.Sequence(
    puppi + egmPhotonIsolationMiniAOD + egmPhotonIsolationMiniAODPUPPI
)

ggPhotonIso_CITK_Sequence = cms.Sequence(
    egmPhotonIsolationMiniAOD
)
