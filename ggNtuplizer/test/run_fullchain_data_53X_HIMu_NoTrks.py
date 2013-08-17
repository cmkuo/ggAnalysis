import FWCore.ParameterSet.Config as cms

from RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff import *

process = cms.Process("ggKIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.GlobalTag.globaltag = cms.string('GR_P_V41_AN1::All')
        
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'file:/data4/cmkuo/testfiles/DoublePhoton_2012C_PRv2_AOD.root'
    '/store/hidata/HIRun2013/PAHighPt/RECO/PromptReco-v1/000/211/631/00000/50B4BBE5-4B75-E211-B8BC-001D09F2960F.root'
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

# HeavyIon
from HeavyIonsAnalysis.Configuration.CommonFunctions_cff import *
overrideCentrality(process)

process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowersPlusTrunc"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("pACentrality"),
    pPbRunFlip = cms.untracked.uint32(211313)
    )

process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Appeltel.RpPbAnalysis.PAPileUpVertexFilter_cff')

# Jets
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.patJetCorrFactors.useRho = cms.bool(True)

process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets25 = process.kt4PFJets.clone()
process.kt6PFJets25.rParam = 0.6
process.kt6PFJets25.doRhoFastjet = True
process.kt6PFJets25.doAreaFastjet = True
process.kt6PFJets25.Rho_EtaMax = 2.5
process.kt6PFJets25.voronoiRfact = 0.9
process.fjSequence25 = cms.Sequence( process.kt6PFJets25 )

process.kt6PFJets25Neu = process.kt6PFJets.clone()
process.kt6PFJets25Neu.rParam = 0.6
process.kt6PFJets25Neu.doRhoFastjet = True
process.kt6PFJets25Neu.doAreaFastjet = True
process.kt6PFJets25Neu.Rho_EtaMax = 2.5
process.kt6PFJets25Neu.Ghost_EtaMax = 3.1
process.kt6PFJets25Neu.inputEtMin = 0.5
process.kt6PFJets25Neu.voronoiRfact = 0.9
process.fjSequence25Neu = cms.Sequence( process.kt6PFJets25Neu )

process.kt6PFJets44 = process.kt4PFJets.clone()
process.kt6PFJets44.rParam = 0.6
process.kt6PFJets44.doRhoFastjet = True
process.kt6PFJets44.doAreaFastjet = True
process.kt6PFJets44.Rho_EtaMax = 4.4
process.kt6PFJets44.voronoiRfact = 0.9
process.fjSequence44 = cms.Sequence( process.kt6PFJets44 )

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
process.fjSequence = cms.Sequence( process.kt6PFJetsForIsolation )

process.ak5PFJets.doAreaFastjet = True

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

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

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
process.ggNtuplizer.dumpTrks=cms.bool(False)
process.ggNtuplizer.doCentrality=cms.bool(True)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

# electron energy regression
process.load("EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi")

process.load("Configuration.StandardSequences.Services_cff")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
    initialSeed = cms.untracked.uint32(1),
    engineName = cms.untracked.string('TRandom3')
    ),
                                                   )

process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
process.calibratedPatElectrons.isMC = cms.bool(False)
process.calibratedPatElectrons.inputDataset = cms.string("2012Jul13ReReco")
process.calibratedPatElectrons.updateEnergyError = cms.bool(True)
process.calibratedPatElectrons.applyCorrections = cms.int32(10)
process.calibratedPatElectrons.debug = cms.bool(False)

process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#                                  )

# Trigger requirements
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.leptonHLTFilter = copy.deepcopy(process.hltHighLevel)
process.leptonHLTFilter.throw = cms.bool(False)
process.leptonHLTFilter.HLTPaths = ['HLT_PAL1DoubleMuOpen_v*','HLT_PAMu3_v*','HLT_PAMu7_v*','HLT_PAMu12_v*']
    
process.p = cms.Path(
    process.PAcollisionEventSelection*
    process.pileupVertexFilterCutGplus*
    process.pACentrality_step*
    process.leptonHLTFilter*
    process.fjSequence25*
    process.fjSequence25Neu*
    process.fjSequence44*
    process.fjSequence*
    process.ak5PFJets*
    process.patDefaultSequence*
    process.eleIsoSequence*
    process.phoIsoSequence*
    process.eleRegressionEnergy*
    process.calibratedPatElectrons*
    process.ggTriggerSequence*
    process.ggMETFiltersSequence*
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
