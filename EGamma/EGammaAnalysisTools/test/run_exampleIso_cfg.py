import FWCore.ParameterSet.Config as cms

process = cms.Process("ExISO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C61C71DC-8061-E111-AAEB-0025B32036D2.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C49AC4E5-7F61-E111-8FB2-003048F1C836.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/C27CBF93-7E61-E111-AB30-003048F024C2.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/B46B6A67-8D61-E111-82C0-003048F118AA.root',
    '/store/relval/CMSSW_5_1_2/RelValZEE/GEN-SIM-RECO/PU_START50_V15A-v1/0003/36B93F53-7F61-E111-AA27-002481E0D524.root'
    
 
    #'/store/relval/CMSSW_5_2_2/RelValZEE/GEN-SIM-RECO/START52_V4-v1/0003/5A681CA6-9273-E111-B557-003048F117B6.root',
    #'/store/relval/CMSSW_5_2_2/RelValZEE/GEN-SIM-RECO/START52_V4-v1/0003/8AE43C38-9373-E111-A208-003048F117B6.root',
    #'/store/relval/CMSSW_5_2_2/RelValZEE/GEN-SIM-RECO/START52_V4-v1/0004/D4D2BD38-1F74-E111-A5CB-BCAEC53296F8.root'
    
    )
    )




process.load('EGamma.EGammaAnalysisTools.electronIsoProducer_cfi')


process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("myIsohisto_522.root")
    )
process.eleIsoAnalyzer = cms.EDAnalyzer("ElectronIsoAnalyzer",
                                        verbose = cms.untracked.bool(True),
                                        Electrons = cms.InputTag('gsfElectrons'),
                                        IsoValElectrons = cms.VInputTag(cms.InputTag('elePFIso:chIsoForGsfEle'),
                                                                        cms.InputTag('elePFIso:phIsoForGsfEle'),
                                                                        cms.InputTag('elePFIso:nhIsoForGsfEle')),
                                        deltaR = cms.string('03'),
                                        effectiveAreaTarget = cms.string('Data2011'),
                                        rhoIsoInputTag          = cms.InputTag("kt6PFJetsForIsolation","rho"),   ##NOTE: to be changed in the HZZ-4l
                                        )


# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
# all the analyses
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)


#NOTE: special case only for HZZ-4l and synchronization
#from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
#kt6PFJets = kt4PFJets.clone( rParam = cms.double(0.6),
#                             doAreaFastjet = cms.bool(True),
#                             doRhoFastjet = cms.bool(True)
#                             )
#process.kt6PFJetsForIso = kt6PFJets.clone( rParam = cms.double(0.6),
#                                           Rho_EtaMax = cms.double(2.5),
#                                           Ghost_EtaMax = cms.double(2.5),
#                                           )



process.p = cms.Path(process.kt6PFJetsForIsolation * process.elePFIso * process.eleIsoAnalyzer )
#process.elePFIso.verbose = True


