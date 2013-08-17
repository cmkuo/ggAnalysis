import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from ggAnalysis.ggNtuplizer.ggTriggerMatcher_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *

patTriggerEvent.patTriggerMatches  = cms.VInputTag(
    "electronTriggerMatchHLTEle27WP80",
    "electronTriggerMatchHLTPAPhoton15TightCaloIdVL",
    "electronTriggerMatchEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL",
    "electronTriggerMatchEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL",
    "electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",
    "electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",
    "electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ",
    "electronTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4Mass50",
    "electronFilterMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "electronFilterMatchEle8Mass50",
    "electronFilterMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "electronFilterMatchSC4Mass50",
    "electronFilterMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT",
    "electronFilterMatchSC17Mass50",
    "eleTrgFMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "eleTrgFMatchEle8Mass50",
    "eleTrgFMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "eleTrgFMatchSC4Mass50",
    "eleTrgFMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT",
    "eleTrgFMatchSC17Mass50",
    "eleTrgFMatch27WP80TrkIso",
    "phoTrgFMatch26IdIso18IdIso18",
    "phoTrgFMatch26IdIso18IdIso26",
    "phoTrgFMatch26R918R918",
    "phoTrgFMatch26R918R926",
    "phoTrgFMatch26R918IdIso18",
    "phoTrgFMatch26R918IdIso26",
    "phoTrgFMatch26IdIso18R918",
    "phoTrgFMatch26IdIso18R926",
    "phoTrgFMatch36IdIso22IdIso22",
    "phoTrgFMatch36IdIso22IdIso36",
    "phoTrgFMatch36R922R922",
    "phoTrgFMatch36R922R936",
    "phoTrgFMatch36R922IdIso22",
    "phoTrgFMatch36R922IdIso36",
    "phoTrgFMatch36IdIso22R922",
    "phoTrgFMatch36IdIso22R936",
    "phoTrgFMatch26IdIsoXL18IdIsoXL18",
    "phoTrgFMatch26IdIsoXL18IdIsoXL26",
    "phoTrgFMatch26R9T18R9T18",
    "phoTrgFMatch26R9T18R9T26",
    "phoTrgFMatch26R918IdIsoXL18",
    "phoTrgFMatch26R918IdIsoXL26",
    "phoTrgFMatch26IdIsoXL18R918",
    "phoTrgFMatch26IdIsoXL18R926",
    "phoTrgFMatch30EBOnly",
    "photonTriggerMatchHLTPhoton26Photon18",
    "photonTriggerMatchHLTPhoton36Photon22",
    "photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass60",
    "photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass70",
    "photonTriggerMatchHLTPhoton36R9Id85ORCaloId10Iso50Photon22R9Id85ORCaloId10Iso50",
    "phoTrgFMatch26Id10Iso50HcalIso",
    "phoTrgFMatch26R9Id85",
    "phoTrgFMatch18R9Id85",
    "phoTrgFMatch18Id10Iso50TrkIso",
    "phoTrgFMatch18Id10Iso50TrkIsoDouble",
    "phoTrgFMatch26Id10Isol50HcalIsol",
    "phoTrgFMatch18Id10Iso50TrkIsol",
    "phoTrgFMatch18Id10Iso50TrkIsolDouble",
    "muonTriggerMatchHLTIsoMu24eta2p1",
    "muonTriggerMatchHLTIsoMu24",
    "muonTriggerMatchHLTMu17Mu8",
    "muonTriggerMatchHLTMu17TkMu8",
    "muonTriggerMatchHLTMu22TkMu8",
    "muonTriggerMatchHLTMu22Mu8" ,
    "muonTriggerMatchHLTMu17forMu17Mu8",
    "muonTriggerMatchHLTMu8forMu17Mu8",
    "muonTriggerMatchHLTMu17forMu17TkMu8" ,
    "muonTriggerMatchHLTMu8forMu17TkMu8"
    )

#patTrigger.processName = 'HISIGNAL'
#patTriggerEvent.processName = 'RECO'

ggTriggerSequence = cms.Sequence(
    patTrigger * ggTriggerMatcherNoJet * patTriggerEvent
    )

