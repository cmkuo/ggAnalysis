import FWCore.ParameterSet.Config as cms

from RecoMET.METAnalyzers.CSCHaloFilter_cfi import *
CSCTightHaloFilter.taggingMode = cms.bool(True)

#from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *
#HBHENoiseFilter.taggingMode = cms.bool(True)
from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *

from RecoMET.METFilters.hcalLaserEventFilter_cfi import *
hcalLaserEventFilter.taggingMode = cms.bool(True)
#from EventFilter.HcalRawToDigi.hcallaserhbhehffilter2012_cfi import *
#hcallaserhbhehffilter2012.forceFilterTrue = cms.untracked.bool(True)

from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = cms.bool(True)

from RecoMET.METFilters.ecalLaserCorrFilter_cfi import *
ecalLaserCorrFilter.taggingMode = cms.bool(True)

goodVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )

from RecoMET.METFilters.trackingFailureFilter_cfi import *
trackingFailureFilter.taggingMode = cms.bool(True)

from RecoMET.METFilters.trackingPOGFilters_cff import *
manystripclus53X.taggedMode = cms.untracked.bool(True)
manystripclus53X.forcedValue = cms.untracked.bool(False)
toomanystripclus53X.taggedMode = cms.untracked.bool(True)
toomanystripclus53X.forcedValue = cms.untracked.bool(False)
logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)

ggMETFiltersSequence = cms.Sequence(
    #CSCTightHaloFilter*
    #HBHENoiseFilter*
    HBHENoiseFilterResultProducer*
    hcalLaserEventFilter*
    #hcallaserhbhehffilter2012*
    EcalDeadCellTriggerPrimitiveFilter*
    eeBadScFilter*
    ecalLaserCorrFilter*
    goodVertices*
    trackingFailureFilter*
    trkPOGFilters
    )
