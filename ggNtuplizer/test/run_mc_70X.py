import FWCore.ParameterSet.Config as cms

process = cms.Process("ggkIT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START70_V7::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'file:/data4/cmkuo/testfiles/GluGluToHToZgammaTo2L_M-125_13TeV-powheg-pythia8_S14.root'
        'file:/tmp/cmkuo/GluGluToHToZgammaTo2L_M-125_13TeV-powheg-pythia8_S14.root'
        )
                            )

process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

# Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('output.root')
#                                  )

process.p = cms.Path(
    process.patDefaultSequence* 
    process.ggNtuplizer)

#process.out_step = cms.EndPath(process.output)
