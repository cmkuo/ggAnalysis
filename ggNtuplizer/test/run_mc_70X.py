from PhysicsTools.PatAlgos.patTemplate_cfg import *
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.Tracer = cms.Service('Tracer')

#Lvdp
#process.load("PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff")
#process.load("PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi")
#process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')



process.load("PhysicsTools.PatAlgos.patSequences_cff")

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

addJetCollection(
    process,
    labelName = 'AK4PFCHS',
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    algo = 'ak',
    rParam = 0.4,
    jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = [
        'jetBProbabilityBJetTags'
      , 'jetProbabilityBJetTags'
      , 'trackCountingHighPurBJetTags'
      , 'trackCountingHighEffBJetTags'
      , 'simpleSecondaryVertexHighEffBJetTags'
      , 'simpleSecondaryVertexHighPurBJetTags'
      , 'combinedSecondaryVertexBJetTags'
      ],
)
addJetCollection(
    process,
    labelName = 'CA8PFCHS',
    jetSource = cms.InputTag('ca8PFJetsCHS'),
    algo = 'ca',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
    )

addJetCollection(
    process,
    labelName = 'AK8PFCHS',
    jetSource = cms.InputTag('ak8PFJetsCHS'),
    algo = 'ak',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
    )
process.ca8PFJetsCHSPrunedLinks = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                         src = cms.InputTag("ca8PFJetsCHS"),
                                         matched = cms.InputTag("ca8PFJetsCHSPruned"),
                                         distMax = cms.double(0.8),
                                         value = cms.string('mass')
                        )


patJetsAK4 = process.patJetsAK4PFCHS
patJetsCA8 = process.patJetsCA8PFCHS
patJetsAK8 = process.patJetsAK8PFCHS

process.out.outputCommands += ['keep *_selectedPat*_*_*',
                                'keep *_ak4PFJetsCHS_*_*',
                               'keep *_patJetsAK4PFCHS_*_*',
                               'keep *_ca8PFJetsCHS_*_*',
                               'keep *_patJetsCA8PFCHS_*_*',
                               'keep *_ak8PFJetsCHS_*_*',
                               'keep *_patJetsAK8PFCHS_*_*',
				'keep *_*ca8PFJetsCHSPruned*_*_*'
]
patJetsCA8.userData.userFloats.src += ['ca8PFJetsCHSPrunedLinks']
process.out.outputCommands += ['keep *_ca8PFJetsCHSPrunedLinks_*_*']


#### Adding Nsubjetiness

process.NjettinessCA8 = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("ca8PFJetsCHS"),
                              cone=cms.double(0.8)
                              )

patJetsCA8.userData.userFloats.src += ['NjettinessCA8:tau1','NjettinessCA8:tau2','NjettinessCA8:tau3']
process.out.outputCommands += ['keep *_NjettinessCA8_*_*']


process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")

##Taus
from PhysicsTools.PatAlgos.tools.tauTools import *
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.cleanPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.cleanPatTaus.finalCut = cms.string(' pt > 15.0 & abs(eta) < 2.5 ')
process.load("ggAnalysis.ggNtuplizer.ggTau_cff")
##Jets

process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(True)

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        'file:/uscms/home/lovedeep/eos/RSGluonToTT_M-3000_Tune4C_13TeV-pythia8_Spring14dr_PU_S14_POSTLS170_V6-v1_AODSIM.root'
#'root://cmsxrootd-site.fnal.gov//store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU_S14_POSTLS170_V6-v1/00000/0015495C-66CC-E311-8A39-00266CF33100.root'
        )
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.p = cms.Path(
process.patDefaultSequence
 *    process.NjettinessCA8
 *  process.ggNtuplizer
)

process.outpath = cms.EndPath(process.out)
process.outpath.remove(process.out)

