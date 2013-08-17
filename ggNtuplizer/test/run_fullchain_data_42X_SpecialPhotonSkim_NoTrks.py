import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

process = cms.Process("ggKIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_42_V12::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/home/cmkuo/work/testfiles/A827F109-B9A8-E011-8047-001A64789E4C.root'
    ),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Trigger matching
process.load("ggAnalysis.ggNtuplizer.ggPatTriggerMatching_cff")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData(process, ['All'], outputInProcess = False)

# Jets
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.patJetCorrFactors.useRho = cms.bool(True)

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )

process.kt6PFJets44 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets44.Rho_EtaMax = cms.double(4.4)
process.kt6PFJets44.Ghost_EtaMax = cms.double(5.0)
process.fjSequence44 = cms.Sequence( process.kt6PFJets44 )

process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets44','rho')

from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                 cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = False
                 )

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

#process.patJetGenJetMatch.matched = cms.InputTag('iterativeCone5GenJets')

process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
process.ggNtuplizer.doGenParticles = cms.bool(False)
process.ggNtuplizer.jetSrc = cms.InputTag("selectedPatJetsAK5PF")
#process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::REDIGI38XPU")
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.ggNtuplizer.getBlocks=cms.bool(False)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.ggNtuplizer.dumpTrks=cms.bool(False)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#)

# Trigger requirements
process.photonHLTFilter = cms.EDFilter("HLTHighLevel",
                                       eventSetupPathsKey = cms.string(''),
                                       TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                       HLTPaths = cms.vstring("HLT_Photon75_CaloIdVL_v*", "HLT_Photon75_CaloIdVL_IsoL_v*", "HLT_Photon90_CaloIdVL_IsoL_v*", "HLT_Photon125_v*", "HLT_Photon135_v*", "HLT_DoublePhoton33_v*", "HLT_DoublePhoton50_v*", "HLT_DoublePhoton60_v*", "HLT_DoublePhoton80_v*", "HLT_Photon26_Photon18_v*"),
                                       andOr = cms.bool(True), 
                                       throw = cms.bool(False)
                                       )

process.p = cms.Path(
    process.photonHLTFilter*
    process.fjSequence25*
    process.fjSequence44*
    process.patDefaultSequence*
    process.ggTriggerSequence*
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
