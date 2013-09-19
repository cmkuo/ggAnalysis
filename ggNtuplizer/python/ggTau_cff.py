import FWCore.ParameterSet.Config as cms

# ------------------------------------------------------------------------------------
# This python cfi file contains the HPS Tau ID sources for 2012 analysis.
# ------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------
# Out-of-the-box available IDs are:
# 'againstElectronDeadECAL' 'againstElectronLoose' 'againstElectronLooseMVA3' 'againstElectronMVA3category' 'againstElectronMVA3raw'
# 'againstElectronMedium' 'againstElectronMediumMVA3' 'againstElectronTight' 'againstElectronTightMVA3' 'againstElectronVTightMVA3'
# 'againstMuonLoose' 'againstMuonLoose2' 'againstMuonMedium' 'againstMuonMedium2' 'againstMuonTight' 'againstMuonTight2'
# 'byCombinedIsolationDeltaBetaCorrRaw' 'byCombinedIsolationDeltaBetaCorrRaw3Hits' 'byIsolationMVA2raw' 'byIsolationMVAraw'
# 'byLooseCombinedIsolationDeltaBetaCorr' 'byLooseCombinedIsolationDeltaBetaCorr3Hits' 'byLooseIsolationMVA' 'byLooseIsolationMVA2'
# 'byMediumCombinedIsolationDeltaBetaCorr' 'byMediumCombinedIsolationDeltaBetaCorr3Hits' 'byMediumIsolationMVA' 'byMediumIsolationMVA2'
# 'byTightCombinedIsolationDeltaBetaCorr' 'byTightCombinedIsolationDeltaBetaCorr3Hits' 'byTightIsolationMVA' 'byTightIsolationMVA2'
# 'byVLooseCombinedIsolationDeltaBetaCorr' 'decayModeFinding' .
# ------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import *

patTaus.tauIDSources.byVLooseIsolation              = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation")
patTaus.tauIDSources.byLooseIsolation               = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")
patTaus.tauIDSources.byMediumIsolation              = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation")
patTaus.tauIDSources.byTightIsolation               = cms.InputTag("hpsPFTauDiscriminationByTightIsolation")
patTaus.tauIDSources.byVLooseIsolationDeltaBetaCorr = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr")
patTaus.tauIDSources.byLooseIsolationDeltaBetaCorr  = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr")
patTaus.tauIDSources.byMediumIsolationDeltaBetaCorr = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr")
patTaus.tauIDSources.byTightIsolationDeltaBetaCorr  = cms.InputTag("hpsPFTauDiscriminationByTightIsolationDBSumPtCorr")
#patTaus.tauIDSources.byIsolationMVAraw             = cms.InputTag("hpsPFTauDiscriminationByIsolationMVAraw")     #included
#patTaus.tauIDSources.byLooseIsolationMVA           = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA")   #included
#patTaus.tauIDSources.byMediumIsolationMVA          = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA")  #included
#patTaus.tauIDSources.byTightIsolationMVA           = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA")   #included
patTaus.tauIDSources.againstElectronMVA2raw         = cms.InputTag("hpsPFTauDiscriminationByMVA2rawElectronRejection")
patTaus.tauIDSources.againstElectronMVA2category    = cms.InputTag("hpsPFTauDiscriminationByMVA2rawElectronRejection:category")
patTaus.tauIDSources.againstElectronVLooseMVA2      = cms.InputTag("hpsPFTauDiscriminationByMVA2VLooseElectronRejection")
patTaus.tauIDSources.againstElectronLooseMVA2       = cms.InputTag("hpsPFTauDiscriminationByMVA2LooseElectronRejection")
patTaus.tauIDSources.againstElectronMediumMVA2      = cms.InputTag("hpsPFTauDiscriminationByMVA2MediumElectronRejection")
patTaus.tauIDSources.againstElectronTightMVA2       = cms.InputTag("hpsPFTauDiscriminationByMVA2TightElectronRejection")
patTaus.tauIDSources.againstElectronMVA             = cms.InputTag("hpsPFTauDiscriminationByMVAElectronRejection")
patTaus.tauIDSources.againstMuonLoose3              = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3")
#patTaus.tauIDSources.againstMuonMedium3            = cms.InputTag("hpsPFTauDiscriminationByMediumMuonRejection3") #does not exist
patTaus.tauIDSources.againstMuonTight3              = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3")
