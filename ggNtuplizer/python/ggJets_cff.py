import FWCore.ParameterSet.Config as cms

from QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff import *
QGTagger.srcJets = cms.InputTag("selectedPatJetsAK5PF")
QGTagger.isPatJet  = cms.untracked.bool(True)

from ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff import *
ca8Jets = cms.Sequence( PATCMGJetSequenceCA8CHS + PATCMGJetSequenceCA8CHSpruned + selectedPatJetsCA8CHSwithNsub)

