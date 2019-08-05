import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')

#process.Tracer = cms.Service("Tracer")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:/data4/cmkuo/testfiles/ggHZg_RunIISummer16MiniAODv3.root'
        ))

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

### L1 ECAL prefiring
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2016BtoH"), 
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)

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

### update JEC
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
    payload = 'AK4PFchs') 

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

# MET correction and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=False,
                           postfix = "ModifiedMET"
                           )

# random generator for jet smearing
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   ggNtuplizer  = cms.PSet(
        initialSeed = cms.untracked.uint32(201678),
        engineName = cms.untracked.string('TRandom3')
        )
                                                   )

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.ggNtuplizer.year=cms.int32(2016)
process.ggNtuplizer.doGenParticles=cms.bool(True)
process.ggNtuplizer.runL1ECALPrefire=cms.bool(True)
process.ggNtuplizer.dumpPFPhotons=cms.bool(True)
process.ggNtuplizer.dumpHFElectrons=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpAK8Jets=cms.bool(False)
process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.patTriggerResults=cms.InputTag("TriggerResults", "", "PAT")
process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger", "", "")
process.ggNtuplizer.ak4JetSrc=cms.InputTag("slimmedJetsJEC")
process.ggNtuplizer.pfMETLabel=cms.InputTag("slimmedMETsModifiedMET")

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("slimmedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))

process.p = cms.Path(
    process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq *
    process.cleanedMu *
    process.jetCorrFactors *
    process.slimmedJetsJEC *
    process.prefiringweight *
    process.ggNtuplizer
    )

#print process.dumpPython()
