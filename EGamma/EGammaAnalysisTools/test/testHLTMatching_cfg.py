import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST") 
process.load("EGamma.EGammaAnalysisTools.electronHLTMatching_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend(
    (
    'file:/home/users/matteo/Hgg_44.root',
    )
    )

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_*_TEST'),
                               fileName = cms.untracked.string('electrons.root')
                               )

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START44_V5::All"

process.myPath = cms.Sequence(process.electronHLTMatching)

process.p11 = cms.Path(process.myPath)
process.outpath = cms.EndPath(process.out)
