import FWCore.ParameterSet.Config as cms

process = cms.Process("ggkIT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START53_V7A::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(100)
            )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/uscms/home/makouski/nobackup/TTJets_SemiLeptMGDecays_8TeV-madgraph.root'
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

process.ak5PFJets.doAreaFastjet = True
process.patJets.addTagInfos = cms.bool(True)

# Taus
from PhysicsTools.PatAlgos.tools.tauTools import *
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.cleanPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.cleanPatTaus.finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5 ')
process.load("ggAnalysis.ggNtuplizer.ggTau_cff")

#process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
#process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets25','rho')

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

# PF isolations for electrons and muons
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_crab_cfi")
process.ggNtuplizer.jetSrc = cms.InputTag("selectedPatJetsAK5PF")
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults::HLT")
process.ggNtuplizer.getBlocks=cms.bool(False)
process.ggNtuplizer.useAllPF=cms.bool(False)
process.ggNtuplizer.dumpTrks=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

#Lvdp
# gluontagger
process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')  
process.QGTagger.srcJets = cms.InputTag("selectedPatJetsAK5PF")
process.QGTagger.isPatJet  = cms.untracked.bool(True)
#process.QGTagger.useCHS  = cms.untracked.bool(True)

#jet substructure
from RecoJets.Configuration.RecoGenJets_cff import ak7GenJetsNoNu
process.ca8GenJetsNoNu = ak7GenJetsNoNu.clone()
process.ca8GenJetsNoNu.rParam = 0.8
process.ca8GenJetsNoNu.jetAlgorithm = "CambridgeAachen"


from CommonTools.ParticleFlow.pfNoPileUp_cff import * 

# JETS  CA8 ----------------------------
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
process.ca8PFJetsCHS = ak5PFJets.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(10.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )
jetSource = 'ca8PFJetsCHS'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
process.patJetCorrFactorsCA8CHS = patJetCorrFactors.clone()
process.patJetCorrFactorsCA8CHS.src = jetSource
# will need to add L2L3 corrections in the cfg
process.patJetCorrFactorsCA8CHS.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
process.patJetCorrFactorsCA8CHS.payload = 'AK7PFchs'
process.patJetCorrFactorsCA8CHS.useRho = True
# parton and gen jet matching

from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
process.patJetPartonMatchCA8CHS = patJetPartonMatch.clone()
process.patJetPartonMatchCA8CHS.src = jetSource
process.patJetGenJetMatchCA8CHS = patJetGenJetMatch.clone()
process.patJetGenJetMatchCA8CHS.src = jetSource
process.patJetGenJetMatchCA8CHS.matched = 'ca8GenJetsNoNu'

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import *
process.patJetPartonAssociationCA8CHS = patJetPartonAssociation.clone()
process.patJetPartonAssociationCA8CHS.jets = jetSource

# pat jets

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import *
process.patJetsCA8CHS = patJets.clone()
process.patJetsCA8CHS.jetSource = jetSource
process.patJetsCA8CHS.addJetCharge = False
process.patJetsCA8CHS.embedCaloTowers = False
process.patJetsCA8CHS.embedPFCandidates = False
process.patJetsCA8CHS.addAssociatedTracks = False
process.patJetsCA8CHS.addBTagInfo = False
process.patJetsCA8CHS.addDiscriminators = False
process.patJetsCA8CHS.getJetMCFlavour = False
process.patJetsCA8CHS.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHS'))
process.patJetsCA8CHS.genPartonMatch = cms.InputTag('patJetPartonMatchCA8CHS')
process.patJetsCA8CHS.genJetMatch = cms.InputTag('patJetGenJetMatchCA8CHS')

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
process.selectedPatJetsCA8CHS = selectedPatJets.clone()
process.selectedPatJetsCA8CHS.src = 'patJetsCA8CHS'
process.selectedPatJetsCA8CHS.cut = 'pt()>20'

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.jetMCSequenceCA8CHS = cms.Sequence(
    process.patJetPartonMatchCA8CHS +
    genParticlesForJetsNoNu +
    process.ca8GenJetsNoNu +
    process.patJetGenJetMatchCA8CHS
    )

process.PATCMGJetSequenceCA8CHS = cms.Sequence(
    process.ca8PFJetsCHS +
    process.jetMCSequenceCA8CHS +
    process.patJetCorrFactorsCA8CHS +
    process.patJetsCA8CHS +
    process.selectedPatJetsCA8CHS
    )
# JETS PRUNED CA8 ----------------------------

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ca8PFJetsCHSpruned = ak5PFJetsPruned.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(10.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )

jetSource = 'ca8PFJetsCHSpruned'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
process.patJetCorrFactorsCA8CHSpruned = patJetCorrFactors.clone()
process.patJetCorrFactorsCA8CHSpruned.src = jetSource
# will need to add L2L3 corrections in the cfg
process.patJetCorrFactorsCA8CHSpruned.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
process.patJetCorrFactorsCA8CHSpruned.payload = 'AK7PFchs'
process.patJetCorrFactorsCA8CHSpruned.useRho = True
# parton and gen jet matching

from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
process.patJetPartonMatchCA8CHSpruned = patJetPartonMatch.clone()
process.patJetPartonMatchCA8CHSpruned.src = jetSource
process.patJetGenJetMatchCA8CHSpruned = patJetGenJetMatch.clone()
process.patJetGenJetMatchCA8CHSpruned.src = jetSource
process.patJetGenJetMatchCA8CHSpruned.matched = 'ca8GenJetsNoNu'

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import *
process.patJetPartonAssociationCA8CHSpruned = patJetPartonAssociation.clone()
process.patJetPartonAssociationCA8CHSpruned.jets = jetSource

# pat jets

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import *
process.patJetsCA8CHSpruned = patJets.clone()
process.patJetsCA8CHSpruned.jetSource = jetSource
process.patJetsCA8CHSpruned.addJetCharge = False
process.patJetsCA8CHSpruned.embedCaloTowers = False
process.patJetsCA8CHSpruned.embedPFCandidates = False
process.patJetsCA8CHSpruned.addAssociatedTracks = False
process.patJetsCA8CHSpruned.addBTagInfo = False
process.patJetsCA8CHSpruned.addDiscriminators = False
process.patJetsCA8CHSpruned.getJetMCFlavour = False
process.patJetsCA8CHSpruned.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHSpruned'))
process.patJetsCA8CHSpruned.genPartonMatch = cms.InputTag('patJetPartonMatchCA8CHSpruned')
process.patJetsCA8CHSpruned.genJetMatch = cms.InputTag('patJetGenJetMatchCA8CHSpruned')

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
process.selectedPatJetsCA8CHSpruned = selectedPatJets.clone()
process.selectedPatJetsCA8CHSpruned.src = 'patJetsCA8CHSpruned'
process.selectedPatJetsCA8CHSpruned.cut = 'pt()>20'

process.jetMCSequenceCA8CHSpruned = cms.Sequence(
    process.patJetPartonMatchCA8CHSpruned +
    process.patJetGenJetMatchCA8CHSpruned
    )

process.PATCMGJetSequenceCA8CHSpruned = cms.Sequence(
    process.ca8PFJetsCHSpruned +
    process.jetMCSequenceCA8CHSpruned +
    process.patJetCorrFactorsCA8CHSpruned +
    process.patJetsCA8CHSpruned +
    process.selectedPatJetsCA8CHSpruned
    )


#### Adding Nsubjetiness

process.selectedPatJetsCA8CHSwithNsub = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("selectedPatJetsCA8CHS"),
                              cone=cms.double(0.8)
                              )

ca8Jets = cms.Sequence( process.PATCMGJetSequenceCA8CHS + process.PATCMGJetSequenceCA8CHSpruned + process.selectedPatJetsCA8CHSwithNsub)

#end Lvdp


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
process.calibratedPatElectrons.isMC = cms.bool(True)
process.calibratedPatElectrons.inputDataset = cms.string("Summer12_DR53X_HCP2012")
process.calibratedPatElectrons.updateEnergyError = cms.bool(True)
process.calibratedPatElectrons.applyCorrections = cms.int32(10)
process.calibratedPatElectrons.debug = cms.bool(False)

process.load("ggAnalysis.ggNtuplizer.ggRhoFastJet_cff")
process.load("ggAnalysis.ggNtuplizer.ggEleID_cff")

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
#)

process.p = cms.Path(
    process.fjSequence*
    process.ak5PFJets*
     process.pfNoPileUpSequence* ###########
    process.pfParticleSelectionSequence*
    process.ggBoostedEleModIsoSequence*
    process.eleMVAID*
    process.patDefaultSequence*
    process.eleIsoSequence*
    process.phoIsoSequence*
    ca8Jets* ###########
    process.eleRegressionEnergy*
    process.calibratedPatElectrons*
    process.ggTriggerSequence*
    process.recoTauClassicHPSSequence*
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
