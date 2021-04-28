import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v10')

#process.Tracer = cms.Service("Tracer")
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'file:data16_94X.root'
        ))

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )


### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True, 
                       era='2016-Legacy',
                       eleIDModules=['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
                                     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
                       )  

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
#runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

### update JEC
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'],
    payload = 'AK4PFchs') 

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

# MET correction and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=True,
                           postfix = "ModifiedMET"
                           )


process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.ggNtuplizer.year=cms.int32(2016)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpPFPhotons=cms.bool(True)
process.ggNtuplizer.dumpHFElectrons=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpAK8Jets=cms.bool(False)
process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.ak4JetSrc=cms.InputTag("slimmedJetsJEC")
process.ggNtuplizer.pfMETLabel=cms.InputTag("slimmedMETsModifiedMET")
process.ggNtuplizer.patTriggerResults=cms.InputTag("TriggerResults", "", "DQM")
#process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger", "", "")
process.ggNtuplizer.addFilterInfoMINIAOD=cms.bool(True)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("slimmedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))

from  PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
process.jetCorrFactorsNano = patJetCorrFactors.clone(src='slimmedJets',
    levels = cms.vstring('L1FastJet',
        'L2Relative',
        'L3Absolute',
	'L2L3Residual'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
from  PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
process.updatedJets = updatedPatJets.clone(
	addBTagInfo=False,
	jetSource='slimmedJets',
	jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactorsNano") ),
)
process.bJetVars = cms.EDProducer("JetRegressionVarProducer",
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("updatedJets"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    gpsrc = cms.InputTag("prunedGenParticles"),
    #musrc = cms.InputTag("slimmedMuons"),
    #elesrc = cms.InputTag("slimmedElectrons")
)
process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
        src = cms.InputTag("updatedJets"),
        userFloats = cms.PSet(
            #leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
            #leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
            #leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
            #leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
            #leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
            #leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
            #leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
            #leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
            #leptonPt = cms.InputTag("bJetVars:leptonPt"),
            vtxPt = cms.InputTag("bJetVars:vtxPt"),
            vtxMass = cms.InputTag("bJetVars:vtxMass"),
            vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
            vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
            #ptD = cms.InputTag("bJetVars:ptD"),
            #genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
            #qgl = cms.InputTag('qgtagger:qgLikelihood'),
            #puId94XDisc = cms.InputTag('pileupJetId94X:fullDiscriminant'),
            #puId102XDisc = cms.InputTag('pileupJetId102X:fullDiscriminant'),
            #chFPV0EF = cms.InputTag("jercVars:chargedFromPV0EnergyFraction"),
            #chFPV1EF = cms.InputTag("jercVars:chargedFromPV1EnergyFraction"),
            #chFPV2EF = cms.InputTag("jercVars:chargedFromPV2EnergyFraction"),
            #chFPV3EF = cms.InputTag("jercVars:chargedFromPV3EnergyFraction"),
            ),
        userInts = cms.PSet(
            #tightId = cms.InputTag("tightJetId"),
            #tightIdLepVeto = cms.InputTag("tightJetIdLepVeto"),
            vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
            #leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
            ),
        )
process.jetSequence = cms.Sequence(process.jetCorrFactorsNano+process.updatedJets+process.bJetVars+process.updatedJetsWithUserData)
process.ggNtuplizer.nanoUpdatedUserJetsLabel=cms.string('updatedJetsWithUserData')

process.p = cms.Path(
    process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq *
    process.cleanedMu *
    process.jetCorrFactors *
    process.slimmedJetsJEC *
    process.jetSequence*
    process.ggNtuplizer
    )

#process.out = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string('recoBPHanalysis_withFilter.root'),
#    outputCommands = cms.untracked.vstring(
#        'keep *',
#        #"drop *",
#        #"keep *_selectedMuons_MuonPreselectionEfficiencyBoolInt_myVertexingProcedure",
#        #"keep *_selectedTracks_TrackPreselectionEfficiencyBoolInt_myVertexingProcedure",
#        #"keep *_tktkVertexingProducer_tktkVertexingEfficiencyBoolInt_myVertexingProcedure",
#        #"keep *_mumuVertexingProducer_mumuVertexingEfficiencyBoolInt_myVertexingProcedure",
#        #"keep *_mumuVertexingProducer_*_myVertexingProcedure",
#        #"keep *_tktkVertexingProducer_*_myVertexingProcedure",
#        #"keep *_generalV0Candidates_*_RECO",
#        #"keep *_offlineBeamSpot_*_RECO",
#        #"keep *_offlinePrimaryVertices_*_RECO",
#        #"keep *_genParticles__HLT",
#        #"keep *_TriggerResults__HLT",
#    ),
#    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
#)

#process.e = cms.EndPath(process.out)