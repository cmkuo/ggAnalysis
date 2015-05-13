import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = cms.string('GR_P_V49B::All')
process.GlobalTag.globaltag = cms.string('GR_R_74_V8::All')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'file:/tmp/cmkuo/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_Phys14.root'
        #'file:/data4/cmkuo/testfiles/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_Phys14.root'
        #'root://cmsxrootd.fnal.gov///store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
        "root://eoscms//eos/cms/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleElectron/RECO/GR_R_74_V8_1Apr_RelVal_zEl2012D-v1/00000/C04717C4-48D9-E411-9E88-002618943901.root"
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/1EB73E6B-0AF4-E411-8AAB-02163E012276.root'
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/229CC9FA-0BF4-E411-852D-02163E0137FA.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/22ECB7D5-06F4-E411-831D-02163E0139C6.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/328F8702-0CF4-E411-A80D-02163E0135BE.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/36B7140D-0CF4-E411-BA17-02163E01354B.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/36C8B0E6-07F4-E411-9B0C-02163E012332.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/3AEBD3A8-0AF4-E411-A403-02163E0129ED.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/4AA07411-05F4-E411-AD92-02163E01370E.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/4C6A1AEA-05F4-E411-8EEA-02163E012324.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/5685DD77-0DF4-E411-B0B9-02163E01354B.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/98B718CC-09F4-E411-A471-02163E013780.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/BE2329D7-06F4-E411-8590-02163E012847.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/C67F2E22-0CF4-E411-A6EB-02163E012AE4.root',
#'/store/express/Commissioning2015/ExpressPhysics/FEVT/Express-v1/000/243/639/00000/E2FCCA72-0DF4-E411-B071-02163E0137FA.root'
)
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])



process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(False)
process.ggNtuplizer.dumpJets=cms.bool(False)
process.ggNtuplizer.dumpTaus=cms.bool(False)
process.ggNtuplizer.isAOD=cms.bool(useAOD)



process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

###9th April, 2015 SHILPI JAIN
###ADDING ELECTRON ID MVA IN THE PATH: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2

#process.load("EgammaAnalysis.ElectronTools.electronIdMVAProducer_CSA14_cfi")

#setattr(process.patElectrons.electronIDSources, "trigMVAid",cms.InputTag("mvaTrigV050nsCSA14"))

useAOD = True
#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD


switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)



process.p = cms.Path(
    ###process.egmGsfElectronIDSequence
    # process.mvaTrigV050nsCSA14
    # + process.mvaTrigV025nsCSA14
    # + process.mvaNonTrigV050nsCSA14
    # + process.mvaNonTrigV025nsCSA14
    # + process.mvaNonTrigV025nsPHYS14 
#    process.patDefaultSequence *
    process.egmGsfElectronIDSequence 
    * process.egmPhotonIDSequence 
    * process.ggNtuplizer
    )


#print process.dumpPython()
