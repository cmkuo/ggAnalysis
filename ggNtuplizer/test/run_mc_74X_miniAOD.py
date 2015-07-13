import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'file:/tmp/cmkuo/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_Phys14.root'
        #'file:/data4/cmkuo/testfiles/GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_Phys14.root'
        #'root://cmsxrootd.fnal.gov///store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00CC714A-F86B-E411-B99A-0025904B5FB8.root'
'/store/relval/CMSSW_7_4_0_pre9_ROOT6/RelValWpToENu_M-2000_13TeV/MINIAODSIM/MCRUN2_74_V7-v1/00000/4A75C5D1-DCD1-E411-BE48-002618943951.root'
#        "/store/relval/CMSSW_7_4_0_pre9_ROOT6/DoubleElectron/RECO/GR_R_74_V8_1Apr_RelVal_zEl2012D-v1/00000/C04717C4-48D9-E411-9E88-002618943901.root"
)
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])



process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
process.ggNtuplizer.doGenParticles=cms.bool(True)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

###9th April, 2015 SHILPI JAIN
###ADDING ELECTRON ID MVA IN THE PATH: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2

#process.load("EgammaAnalysis.ElectronTools.electronIdMVAProducer_CSA14_cfi")

#setattr(process.patElectrons.electronIDSources, "trigMVAid",cms.InputTag("mvaTrigV050nsCSA14"))

useAOD = False
#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

process.ggNtuplizer.isAOD=cms.bool(useAOD)

switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V0_cff']

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
