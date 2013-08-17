import FWCore.ParameterSet.Config as cms

# 2012 MVA eID
from EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi import *
#eleMVAID = cms.Sequence(mvaTrigV0 + mvaTrigNoIPV0 + mvaNonTrigV0)
eleMVAID = cms.Sequence(mvaTrigV0 + mvaNonTrigV0)

# modified electron isolation
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
                                     eleLabel = cms.InputTag("gsfElectrons"),
                                     barrelCuts = cms.PSet(heepBarrelCuts),
                                     endcapCuts = cms.PSet(heepEndcapCuts),
                                     eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
                                     eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
                                     verticesLabel = cms.InputTag("offlinePrimaryVertices"),
                                     applyRhoCorrToEleIsol = cms.bool(True),
                                     writeIdAsInt = cms.bool(True)
                                     )
heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")

heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer", cutValueMap = cms.InputTag("heepIdNoIso"),
                                         inputGsfEles = cms.InputTag("gsfElectrons")  )

# Boosted Z ModEleIso: 1b) Calculating the modified iso. values using BstdZeeTools EDProducer

from TSWilliams.BstdZeeTools.bstdzeemodisolproducer_cff import *

modElectronIso = cms.EDProducer("BstdZeeModIsolProducer",
                                        bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )

ggBoostedEleModIsoSequence = cms.Sequence(
    heepIdNoIso*
    heepIdNoIsoEles*
    modElectronIso
    )
