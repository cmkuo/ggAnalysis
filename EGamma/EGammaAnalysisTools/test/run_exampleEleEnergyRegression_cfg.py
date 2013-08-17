import FWCore.ParameterSet.Config as cms

process = cms.Process("ExREG")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V10::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('root://pcmssd12//data/mangano/MC/8TeV/hzz/reco/ggHToZZTo4L_M-120_Summer12_S7.003048678E92.root')
    )


process.load('EGamma.EGammaAnalysisTools.electronRegressionEnergyProducer_cfi')

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_*_ExREG'),
                               fileName = cms.untracked.string('electrons.root')
                                                              )
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.p = cms.Path(process.eleRegressionEnergy)
process.outpath = cms.EndPath(process.out)


