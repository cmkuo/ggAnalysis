## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
## switch to uncheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)
## load tau sequences up to selectedPatJets
process.load("PhysicsTools.PatAlgos.producersLayer1.jetProducer_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi")
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

addJetCollection(
    process,
    labelName = 'AK4PFCHS',
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    algo = 'ak4',
    rParam = 0.4,
    jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
    )
addJetCollection(
    process,
    labelName = 'CA8PFCHS',
    jetSource = cms.InputTag('ca8PFJetsCHS'),
    algo = 'ca8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
    )

addJetCollection(
    process,
    labelName = 'AK8PFCHS',
    jetSource = cms.InputTag('ak8PFJetsCHS'),
    algo = 'ak8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')
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
                               'keep *_patJetsAK8PFCHS_*_*']


