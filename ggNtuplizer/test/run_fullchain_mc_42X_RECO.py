import FWCore.ParameterSet.Config as cms

process = cms.Process("ggkIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START42_V14B::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    'file:/cms/data22/kandeen/root/H130GGgluonfusion_7TeV/H130GGgluonfusion_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_21000_22000.root'
    #'file:/data3/ncuhep/423_vgamma_42x_v4/testfiles/WToENu_TuneZ2_7TeV_Summer11_AOD.root'
    #'file:/data2/cmkuo/STEP2_RAW2DIGI_L1Reco_RECO_PU.root'
    'rfio:/castor/cern.ch/user/c/cmkuo/ggNtuple/edmFiles/3CE75CB9-317B-E011-86BE-002618943864.root'
    #'rfio:/castor/cern.ch/cms/store/relval/CMSSW_4_2_3/RelValH130GGgluonfusion/GEN-SIM-RECO/START42_V12-v2/0067/0A2A882D-087C-E011-B6E4-00248C0BE016.root',
    ), 
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Trigger matching
process.load("ggAnalysis.ggNtuplizer.ggPatTriggerMatching_cff")

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
process.ggNtuplizer.jetSrc = cms.InputTag("selectedPatJetsAK5PF")
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.ggNtuplizer.getBlocks=cms.bool(True)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#)

process.p = cms.Path(
    process.fjSequence25*
    process.fjSequence44*
    process.patDefaultSequence*
    process.ggTriggerSequence* 
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
