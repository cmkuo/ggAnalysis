import FWCore.ParameterSet.Config as cms

from RecoMET.METFilters.CSCTightHaloFilter_cfi import *
CSCTightHaloFilter.taggingMode = cms.bool(True)

from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = cms.bool(True)

primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                   vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                   minimumNDOF = cms.uint32(4) ,
                                   maxAbsZ = cms.double(24),
                                   maxd0 = cms.double(2),
                                   filter = cms.bool(False)
                                   )

from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

ggMETFiltersSequence = cms.Sequence(
     CSCTightHaloFilter*
     HBHENoiseFilterResultProducer*
     primaryVertexFilter*
     eeBadScFilter*
     EcalDeadCellTriggerPrimitiveFilter
)
