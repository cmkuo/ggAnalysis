import FWCore.ParameterSet.Config as cms

process = cms.Process("ExISO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/data3/ncuhep/44X_vgamma/testfiles/WGToENuG_START44_V9B-v1.root'),
                            )


process.load('EGamma.EGammaAnalysisTools.photonIsoProducer_cfi')
process.phoPFIso.verbose = True

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myIsohisto_522.root")
    )


process.p = cms.Path(process.phoPFIso)
#process.elePFIso.verbose = True


