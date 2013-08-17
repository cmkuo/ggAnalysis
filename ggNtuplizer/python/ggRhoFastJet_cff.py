import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.kt4PFJets_cfi import *

kt6PFJets25 = kt4PFJets.clone()
kt6PFJets25.rParam = 0.6
kt6PFJets25.doRhoFastjet = True
kt6PFJets25.doAreaFastjet = True
kt6PFJets25.Rho_EtaMax = 2.5
kt6PFJets25.voronoiRfact = 0.9

kt6PFJets25Neu = kt4PFJets.clone()
kt6PFJets25Neu.rParam = 0.6
kt6PFJets25Neu.doRhoFastjet = True
kt6PFJets25Neu.doAreaFastjet = True
kt6PFJets25Neu.Rho_EtaMax = 2.5
kt6PFJets25Neu.Ghost_EtaMax = 3.1
kt6PFJets25Neu.inputEtMin = 0.5
kt6PFJets25Neu.voronoiRfact = 0.9

kt6PFJets44 = kt4PFJets.clone()
kt6PFJets44.rParam = 0.6
kt6PFJets44.doRhoFastjet = True
kt6PFJets44.doAreaFastjet = True
kt6PFJets44.Rho_EtaMax = 4.4
kt6PFJets44.voronoiRfact = 0.9

kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

fjSequence = cms.Sequence( kt6PFJets25 *
                           kt6PFJets25Neu *
                           kt6PFJets44 *
                           kt6PFJetsForIsolation
                           )

