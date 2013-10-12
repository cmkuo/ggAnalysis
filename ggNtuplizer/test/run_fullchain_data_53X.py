import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

process = cms.Process("ggKIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.GlobalTag.globaltag = cms.string('FT_53_V21_AN3::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10027/FA331ACE-0B90-E211-9FE3-00261894393E.root'
    ), 
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Trigger matching
process.load("ggAnalysis.ggNtuplizer.ggPatTriggerMatching_cff")

from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, names=['All'], outputModules=[])
#runOnData(process, ['All'], outputInProcess = False)

# Jets
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.patJetCorrFactors.useRho = cms.bool(True)

process.ak5PFJets.doAreaFastjet = True
process.patJets.addTagInfos = cms.bool(True)

# Taus
from PhysicsTools.PatAlgos.tools.tauTools import *
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.cleanPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.cleanPatTaus.finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5 ')
process.load("ggAnalysis.ggNtuplizer.ggTau_cff")

#process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
#process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets25','rho')

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

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType0pfcp1CorrectedMet = process.pfType1CorrectedMet.clone()
process.pfType0pfcp1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')
    )

process.patMETsType0pfcp1PF = process.patMETsPF.clone()
process.patMETsType0pfcp1PF.metSource = cms.InputTag("pfType0pfcp1CorrectedMet")

process.produceType0MET = cms.Sequence(
    process.pfType0pfcp1CorrectedMet*
    process.patMETsType0pfcp1PF
    )

#process.patJetGenJetMatch.matched = cms.InputTag('iterativeCone5GenJets')

process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

# PF isolations for electrons and muons
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_crab_cfi")
process.ggNtuplizer.doGenParticles = cms.bool(False)
process.ggNtuplizer.jetSrc = cms.InputTag("selectedPatJetsAK5PF")
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.ggNtuplizer.getBlocks=cms.bool(False)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.ggNtuplizer.dumpTrks=cms.bool(True)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

# electron energy regression
process.load("EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi")
process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)

process.load("Configuration.StandardSequences.Services_cff")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
    initialSeed = cms.untracked.uint32(1),
    engineName = cms.untracked.string('TRandom3')
    ),
                                                   )

process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    )

process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
process.calibratedPatElectrons.isMC = cms.bool(False)
process.calibratedPatElectrons.inputDataset = cms.string("22Jan2013ReReco")
process.calibratedPatElectrons.correctionsType = cms.int32(2)
process.calibratedPatElectrons.combinationType = cms.int32(3)

process.load("ggAnalysis.ggNtuplizer.ggRhoFastJet_cff")
process.load("ggAnalysis.ggNtuplizer.ggMergedJets_data_cff")
process.load("ggAnalysis.ggNtuplizer.ggEleID_cff")
process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
process.load("ggAnalysis.ggNtuplizer.ggBoostedEleModIso_cff")

process.patElectrons.userIsolation.user = cms.VPSet(
    cms.PSet(src = cms.InputTag("modElectronIso","track")),
    cms.PSet(src = cms.InputTag("modElectronIso","ecal")),
    cms.PSet(src = cms.InputTag("modElectronIso","hcalDepth1"))
    )

process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    )

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#                                  )

# Trigger requirements
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.leptonHLTFilter = copy.deepcopy(process.hltHighLevel)
process.leptonHLTFilter.throw = cms.bool(False)
process.leptonHLTFilter.HLTPaths = ['HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*','HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*','HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*','HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*','HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v*','HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v*','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*','HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v*','HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*']
   
 
process.p = cms.Path(
    process.leptonHLTFilter*
    process.fjSequence*
    process.ak5PFJets*
    process.pfNoPileUpSequence* ###########
    process.pfParticleSelectionSequence*
    process.ggBoostedEleModIsoSequence*
    process.eleMVAID*
    process.type0PFMEtCorrection*
    process.patDefaultSequence*
    process.produceType0MET*
    process.eleIsoSequence*
    process.phoIsoSequence*
    process.ca8Jets* ###########
    process.QuarkGluonTagger*
    process.eleRegressionEnergy*
    process.calibratedPatElectrons*
    process.ggTriggerSequence*
    process.ggMETFiltersSequence*
    process.recoTauClassicHPSSequence*
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
