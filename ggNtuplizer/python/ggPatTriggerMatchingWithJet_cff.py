import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import *
from ggAnalysis.ggNtuplizer.ggTriggerMatcher_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *

patTriggerEvent.patTriggerMatches  = cms.VInputTag(
    "electronTriggerMatchHLTEle27CaloIdVTCaloIsoTTrkIdTTrkIsoT",
    "electronTriggerMatchHLTEle32CaloIdVTCaloIsoTTrkIdTTrkIsoT",
    "electronTriggerMatchHLTEle27WP80PFMT50",
    "electronTriggerMatchHLTEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVL",
    "electronTriggerMatchEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL",
    "eleTrgFMatchEle17SC8TrkIso",
    "eleTrgFMatchEle17SC8Mass30",
    "eleTrgFMatchEle17Ele8TrkIso",
    "eleTrgFMatchEle17Ele8Mass30",
    "eleTrgFMatchEle32SC17TrkIso",
    "eleTrgFMatchEle32Ele17TrkIso",
    "eleTrgFMatchEle32Ele17HE",
    "eleTrgFMatch27WP80TrkIso",
    "eleTrgFMatch32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMass",
    "eleTrgFMatch32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIso",
    "phoTrgFMatch26IdIso18IdIso18",
    "phoTrgFMatch26IdIso18IdIso26",
    "phoTrgFMatch26IdIso18R926",
    "phoTrgFMatch26IdIso18R918",
    "phoTrgFMatch36IdIso22IdIso36",
    "phoTrgFMatch36IdIso22IdIso22",
    "phoTrgFMatch36IdIso22R936",
    "phoTrgFMatch36IdIso22R922",
    "phoTrgFMatch26IdIsoXL18IdIsoXL26",
    "phoTrgFMatch26IdIsoXL18IdIsoXL18",
    "phoTrgFMatch26IdIsoXL18R926",
    "phoTrgFMatch26IdIsoXL18R918",
    "phoTrgFMatch26IdIsoXL18IdIsoXLMass6026",
    "phoTrgFMatch26IdIsoXL18IdIsoXLMass6018",
    "phoTrgFMatch26IdIsoXL18R9Mass6026",
    "phoTrgFMatch26IdIsoXL18R9Mass6018",
    "phoTrgFMatch30EBOnly",
    "photonTriggerMatchHLTDoublePhoton33",
    "photonTriggerMatchHLTDoublePhoton50",
    "photonTriggerMatchHLTDoublePhoton60",
    "photonTriggerMatchHLTDoublePhoton80",
    "photonTriggerMatchHLTPhoton26Photon18",
    "photonTriggerMatchHLTPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVL",
    "photonTriggerMatchHLTPhoton26R9IdPhoton18R9Id",
    "phoTrgFMatch26Id10Iso50HcalIso",
    "phoTrgFMatch26R9Id85",
    "phoTrgFMatch18R9Id85",
    "phoTrgFMatch18Id10Iso50TrkIso",
    "phoTrgFMatch18Id10Iso50TrkIsoDouble",
    "phoTrgFMatch26Id10Isol50HcalIsol",
    "phoTrgFMatch18Id10Iso50TrkIsol",
    "phoTrgFMatch18Id10Iso50TrkIsolDouble",
    "muonTriggerMatchHLTMu30",
    "muonTriggerMatchHLTMu40",
    "muonTriggerMatchHLTMu40eta2p1",
    "muonTriggerMatchHLTMu13Mu8",
    "muonTriggerMatchHLTMu17Mu8",
    "muonTriggerMatchHLTDoubleMu7",
    "muonTriggerMatchHLTMu17forMu17Mu8" ,
    "muonTriggerMatchHLTMu17forMu17TkMu8" ,
    "jetTriggerMatchHLTJet30",
    "jetTriggerMatchHLTJet60",
    "jetTriggerMatchHLTJet80",
    "jetTriggerMatchHLTJet110",
    "jetTriggerMatchHLTJet150",
    "jetTriggerMatchHLTJet190",
    "jetTriggerMatchHLTJet240",
    "jetTriggerMatchHLTJet300",
    "jetTriggerMatchHLTJet370"
    )

#patTrigger.processName = 'RECO'
#patTriggerEvent.processName = 'RECO'

ggTriggerSequenceWithJet = cms.Sequence(
    patTrigger * ggTriggerMatcher * patTriggerEvent
    )
