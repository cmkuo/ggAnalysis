import FWCore.ParameterSet.Config as cms

from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter

BadPFMuonFilterUpdateDz = BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
)

ggMETFiltersSequence = cms.Sequence(
    BadPFMuonFilterUpdateDz
)
