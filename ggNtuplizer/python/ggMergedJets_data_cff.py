import FWCore.ParameterSet.Config as cms

from QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff import *
QGTagger.srcJets = cms.InputTag("selectedPatJetsAK5PF")
QGTagger.isPatJet  = cms.untracked.bool(True)

from ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff import *
patJetsCA8CHSpruned.addGenJetMatch = cms.bool(False)
patJetsCA8CHSpruned.addGenPartonMatch = cms.bool(False)
patJetsCA8CHSpruned.embedGenJetMatch = cms.bool(False)
patJetsCA8CHSpruned.embedGenPartonMatch = cms.bool(False)
PATCMGJetSequenceCA8CHS = cms.Sequence(
    ca8PFJetsCHS +
#    jetMCSequenceCA8CHS +
    patJetCorrFactorsCA8CHS +
    patJetsCA8CHS +
    selectedPatJetsCA8CHS
    )
PATCMGJetSequenceCA8CHSpruned = cms.Sequence(
    ca8PFJetsCHSpruned +
#    jetMCSequenceCA8CHSpruned +
    patJetCorrFactorsCA8CHSpruned +
    patJetsCA8CHSpruned +
    selectedPatJetsCA8CHSpruned
    )

ca8Jets = cms.Sequence( PATCMGJetSequenceCA8CHS + PATCMGJetSequenceCA8CHSpruned + selectedPatJetsCA8CHSwithNsub)

