import FWCore.ParameterSet.Config as cms

from RecoMET.METFilters.BadPFMuonFilter_cfi import *
BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

from RecoMET.METFilters.BadChargedCandidateFilter_cfi import *
BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

ggMETFiltersSequence = cms.Sequence(
     BadPFMuonFilter *
     BadChargedCandidateFilter 
)
