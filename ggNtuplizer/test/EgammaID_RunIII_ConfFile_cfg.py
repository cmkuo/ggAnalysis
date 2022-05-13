import FWCore.ParameterSet.Config as cms

process = cms.Process("ggNtuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.load("Geometry.HcalEventSetup.CaloTowerTopology_cfi")
#process.load("Configuration.Geometry.GeometryExtended2017_cff")
#process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Geometry.HcalEventSetup.hcalTopologyIdeal_cfi")
process.load("Configuration.Geometry.GeometryExtended2021Reco_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '122X_mcRun3_2021_realistic_v9', '')     # 2022 Winter22MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2024_realistic_v4', '')     # 2024 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2023_realistic_v3', '')     # 2023 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun3_2021_realistic_v3', '')     # 2021 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')     # 2018 MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '')     # 2018 data
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v6', '')     # 2017 data 



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#            'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root'
                                    '/store/mc/Run3Winter22MiniAOD/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/MINIAODSIM/FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/2430000/019c9ef2-86db-4258-a076-bdb5169dc3d0.root'
                                    #'/store/mc/Run3Summer19MiniAOD/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/MINIAODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/130000/00DF0005-F507-2C4B-BF8B-C46342D7194E.root' #Run3 MINIAOD 
#                                    '/store/mc/Run3Summer19DRPremix/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/AODSIM/2021Scenario_106X_mcRun3_2021_realistic_v3-v2/130000/9CAB93A4-1C7A-9F48-A707-6ED2FF2AAC6C.root' #Run3 AOD
#                                    '/store/mc/RunIISummer19UL17RECO/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/AODSIM/106X_mc2017_realistic_v6-v2/70000/FFF299F2-182B-B742-978A-A1E27778EFB0.root'  #Run2 UL17 AOD
#                                    '/store/mc/RunIISummer19UL18RECO/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v1/70001/2A3DF112-2092-AC4D-82E6-1CFBEC62B830.root'  #Run2 UL18 AOD DYJetsToLL
#                                    '/store/mc/RunIISummer19UL16RECOAPV/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_mcRun2_asymptotic_preVFP_v8-v1/270002/FEBBB375-C5E8-7446-AAFB-7BC9225DECB5.root' ##Run2 UL16 AOD DYJetsToLL
                                )
                            )

process.TFileService = cms.Service("TFileService", fileName = cms.string('egmTree.root'))

process.ggNtuplizer = cms.EDAnalyzer('ggNtuplizer',
                                      doGenParticles       = cms.bool(True),
                                      runOnParticleGun     = cms.bool(False),
                                      runOnSherpa          = cms.bool(False),
                                      dumpPDFSystWeight    = cms.bool(False),
                                      dumpGenScaleSystWeights = cms.bool(False),
                                      dumpCrystalinfo      = cms.bool(False), 

                                      VtxLabel             = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                                      VtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                      rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                                      rhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                                      triggerEvent         = cms.InputTag("slimmedPatTrigger", "", ""),
                                      triggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
                                      #patTriggerResults    = cms.InputTag("TriggerResults", "", "PAT"),
                                      patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),

                                      genParticleSrc       = cms.InputTag("prunedGenParticles"),
                                      generatorLabel       = cms.InputTag("generator"),
                                      LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                                      pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                                      newParticles         = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),

#                                     #electronSrc = cms.InputTag('gedGsfElectrons'), ## for AOD
                                      electronSrc            = cms.InputTag("slimmedElectrons"), ## for MINIAOD
                                      photonSrc            = cms.InputTag("slimmedPhotons"), ## for MINIAOD
                                      barrelClusterCollection = cms.InputTag("reducedEgamma","reducedEBEEClusters"),
                                      endcapClusterCollection = cms.InputTag("reducedEgamma","reducedEBEEClusters"),
                                      ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                                      eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                                      esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                                  )

process.p = cms.Path(process.ggNtuplizer)
