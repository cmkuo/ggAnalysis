import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

process = cms.Process("ggKIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = autoCond['mc']
process.GlobalTag.globaltag = cms.string('GR_R_42_V12::All')
#process.GlobalTag.globaltag = cms.string('START42_V4:All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/384/80564560-24BE-E011-9747-003048F117B6.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/383/E6046463-1BBE-E011-AC3A-003048F11114.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/F60729AF-48BD-E011-A44D-003048F024FE.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/CA58D2C7-43BD-E011-8D34-003048F11C58.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/C6A3C1AE-3FBD-E011-AF4F-003048F1C420.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/8EE8561A-37BD-E011-979B-001D09F24489.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/804C1CF3-4CBD-E011-A773-003048D375AA.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/7E8BF01E-51BD-E011-AB87-003048F0258C.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/5A1FF8AE-3FBD-E011-B816-BCAEC532972E.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/286/069F03AF-3FBD-E011-92E2-E0CB4E4408D5.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/276/168694F7-D3BC-E011-9466-001D09F24F65.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/268/EAEF968D-1BBD-E011-AFD1-0019B9F70468.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/268/D68CDDC7-16BD-E011-8ABB-003048F11C58.root',
       '/store/data/Run2011A/Photon/RECO/PromptReco-v5/000/172/268/BACFC1E8-FEBC-E011-9B63-0019B9F72BAA.root'
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
process.ggNtuplizer.getBlocks=cms.bool(True)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#)

process.p = cms.Path(
    process.fjSequence25*
    process.fjSequence44*
    process.patDefaultSequence*
    process.pfBasedPhotonIsoSequence*
    process.ggTriggerSequence* 
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
