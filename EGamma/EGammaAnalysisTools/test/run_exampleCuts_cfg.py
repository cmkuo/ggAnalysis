import FWCore.ParameterSet.Config as cms

process = cms.Process("EX")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START42_V12::All'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#
# input
#

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        'rfio:/castor/cern.ch/user/b/benedet/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_0.root'
    ),
    secondaryFileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

#
# Event output
#

process.load("Configuration.EventContent.EventContent_cff")
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myhisto.root")
)

#
# rho value for isolation
#

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

#
# particle flow isolation
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

#
# example analyzer
#

process.eleIdAnalyzer = cms.EDAnalyzer("EGammaCutBasedEleIdAnalyzer",
    electronsInputTag       = cms.InputTag("gsfElectrons"),
    conversionsInputTag     = cms.InputTag("allConversions"),
    beamSpotInputTag        = cms.InputTag("offlineBeamSpot"),
    rhoIsoInputTag          = cms.InputTag("kt6PFJetsForIsolation", "rho"),
    primaryVertexInputTag   = cms.InputTag("offlinePrimaryVertices"),
    isoValInputTags         = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
    printDebug              = cms.bool(True)
)

#
# analyzer
#

#
# processing
#

process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso * process.eleIdAnalyzer)
process.schedule = cms.Schedule(process.p)

