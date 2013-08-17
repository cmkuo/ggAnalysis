import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

process = cms.Process("ggKIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR_R_42_V21::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/data4/cmkuo/testfiles/Photon-30Nov2011_AOD.root'
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

process.ak5PFJets.doAreaFastjet = True

process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets25','rho')

from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                 cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])),
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
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.ggNtuplizer.getBlocks=cms.bool(False)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#)

# Trigger requirements
process.diphotonHLTFilter = cms.EDFilter("HLTHighLevel",
                                       eventSetupPathsKey = cms.string(''),
                                       TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                       HLTPaths = cms.vstring("HLT_Photon.*_CaloId.*_Iso.*_Photon.*_CaloId.*_Iso.*_.*", "HLT_Photon.*_R9Id.*_Photon.*_R9Id.*_.*", "HLT_Photon.*_CaloId.*_Iso.*_Photon.*_R9Id.*_.*", "HLT_Photon.*_R9Id.*_Photon.*_CaloId.*_Iso.*_.*"),
                                       andOr = cms.bool(True),
                                       throw = cms.bool(False)
                                       )

process.p = cms.Path(
    #process.diphotonHLTFilter*
    process.fjSequence25*
    process.fjSequence44*
    process.ak5PFJets*
    process.patDefaultSequence*
    process.ggTriggerSequence* 
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
