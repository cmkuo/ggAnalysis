##
#ggTau_cff :specify tauID discriminators to be retained
#           in patTau formation.
#
#Lovedeep Kaur, KSU, 09/2014
#
#Note: 
# 1. Adapted from Ming's old ggTau_cff for CMSSW_5X
#################################################################

import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import *

patTaus.tauIDSources.byMVArawMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByMVArawMuonRejection:category")
patTaus.tauIDSources.byMVA5rawElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByMVArawMuonRejection:category")
patTaus.tauIDSources.byIsolationMVA3newDMwoLTraw   = cms.InputTag("hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw:category")
patTaus.tauIDSources.byIsolationMVA3oldDMwLTraw   = cms.InputTag("hpsPFTauDiscriminationByIsolationMVA3oldDMwLTraw:category")
patTaus.tauIDSources.byIsolationMVA3newDMwoLTraw   = cms.InputTag("hpsPFTauDiscriminationByIsolationMVA3newDMwoLTraw:category")
patTaus.tauIDSources.byIsolationMVA3newDMwLTraw   = cms.InputTag("hpsPFTauDiscriminationByIsolationMVA3newDMwLTraw:category")
patTaus.tauIDSources.byLooseElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection")
patTaus.tauIDSources.byMediumElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByMediumElectronRejection")
patTaus.tauIDSources.byTightElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByTightElectronRejection")
patTaus.tauIDSources.byMVA5LooseElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByMVA5LooseElectronRejection")
patTaus.tauIDSources.byMVA5MediumElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByMVA5MediumElectronRejection")
patTaus.tauIDSources.byMVA5TightElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByMVA5TightElectronRejection")
patTaus.tauIDSources.byMVA5VTightElectronRejection   = cms.InputTag("hpsPFTauDiscriminationByMVA5VTightElectronRejection")
patTaus.tauIDSources.byLooseMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection")
patTaus.tauIDSources.byMediumMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection")
patTaus.tauIDSources.byTightMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection")
patTaus.tauIDSources.byLooseMuonRejection3   = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3")
patTaus.tauIDSources.byTightMuonRejection3   = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3")
patTaus.tauIDSources.byMVALooseMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByMVALooseMuonRejection")
patTaus.tauIDSources.byMVAMediumMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByMVAMediumMuonRejection")
patTaus.tauIDSources.byMVATightMuonRejection   = cms.InputTag("hpsPFTauDiscriminationByMVATightMuonRejection")
###patTaus.tauIDSources.pfTausDiscriminationByDecayModeFinding   = cms.InputTag("pfTausDiscriminationByDecayModeFinding")
patTaus.tauIDSources.byVLooseIsolation   = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation")
patTaus.tauIDSources.byVLooseCombinedIsolationDBSumPtCorr   = cms.InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr")
patTaus.tauIDSources.byLooseCombinedIsolationDBSumPtCorr   = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr")
patTaus.tauIDSources.byMediumCombinedIsolationDBSumPtCorr   = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr")
patTaus.tauIDSources.byTightCombinedIsolationDBSumPtCorr   = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr")
patTaus.tauIDSources.byLooseCombinedIsolationDBSumPtCorr3Hits   = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits")
patTaus.tauIDSources.byMediumCombinedIsolationDBSumPtCorr3Hits   = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits")
patTaus.tauIDSources.byTightCombinedIsolationDBSumPtCorr3Hits   = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits")
patTaus.tauIDSources.byVLooseIsolationMVA3newDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwoLT")
patTaus.tauIDSources.byLooseIsolationMVA3newDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3newDMwoLT")
patTaus.tauIDSources.byMediumIsolationMVA3newDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3newDMwoLT")
patTaus.tauIDSources.byTightIsolationMVA3newDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3newDMwoLT")
patTaus.tauIDSources.byVTightIsolationMVA3newDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVA3newDMwoLT")
patTaus.tauIDSources.byVVTightIsolationMVA3newDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwoLT")
patTaus.tauIDSources.byVLooseIsolationMVA3oldDMwLT   = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwLT")
patTaus.tauIDSources.byLooseIsolationMVA3oldDMwLT   = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwLT")
patTaus.tauIDSources.byMediumIsolationMVA3oldDMwLT   = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3oldDMwLT")
patTaus.tauIDSources.byTightIsolationMVA3oldDMwLT   = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwLT")

patTaus.tauIDSources.byVTightIsolationMVA3oldDMwLT   = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwLT")
patTaus.tauIDSources.byVVTightIsolationMVA3oldDMwLT   = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwLT")
patTaus.tauIDSources.byVLooseIsolationMVA3oldDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVA3oldDMwoLT")
patTaus.tauIDSources.byLooseIsolationMVA3oldDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3oldDMwoLT")
patTaus.tauIDSources.byTightIsolationMVA3oldDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3oldDMwoLT")
patTaus.tauIDSources.byVTightIsolationMVA3oldDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVA3oldDMwoLT")
patTaus.tauIDSources.byVVTightIsolationMVA3oldDMwoLT   = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVA3oldDMwoLT")
patTaus.tauIDSources.byLooseIsolationMVA3newDMwLT   = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA3newDMwLT")
patTaus.tauIDSources.byVLooseIsolationMVA3newDMwLT   = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVA3newDMwLT")
patTaus.tauIDSources.byMediumIsolationMVA3newDMwLT   = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA3newDMwLT")
patTaus.tauIDSources.byTightIsolationMVA3newDMwLT   = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA3newDMwLT")
patTaus.tauIDSources.byVTightIsolationMVA3newDMwLT   = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVA3newDMwLT")
patTaus.tauIDSources.byVVTightIsolationMVA3newDMwLT   = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVA3newDMwLT")


##DONT KNOW WHAT TO DO WITH THESE BRANCHES: do not begin with hps...
#reco::PFTauDiscriminator              "pfTausDiscriminationByDecayModeFinding"   ""                "RECO"    
#reco::PFTauDiscriminator              "pfTausDiscriminationByIsolation"   ""                "RECO"    
