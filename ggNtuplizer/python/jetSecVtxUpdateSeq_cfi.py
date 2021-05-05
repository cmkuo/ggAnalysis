import FWCore.ParameterSet.Config as cms

from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
jetCorrFactorsNano = patJetCorrFactors.clone(src='slimmedJets',
    levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
updatedJets = updatedPatJets.clone(
	addBTagInfo=False,
	jetSource='slimmedJets',
	jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsNano") ),
)
bJetVars = cms.EDProducer("JetRegressionVarProducer",
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("updatedJets"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    gpsrc = cms.InputTag("prunedGenParticles"),
)
updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src = cms.InputTag("updatedJets"),
        userFloats = cms.PSet(
            vtxPt = cms.InputTag("bJetVars:vtxPt"),
            vtxMass = cms.InputTag("bJetVars:vtxMass"),
            vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
            vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
            ),
        userInts = cms.PSet(
            vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
            ),
        )
jetSecInfoUpdateSequence = cms.Sequence(jetCorrFactorsNano+updatedJets+bJetVars+updatedJetsWithUserData)
