import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v4')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.Tracer = cms.Service("Tracer")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/user/lovedeep/RunIISpring15MiniAODv2DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia802CDE51A-726D-E511-B2C4-0025905C96EA.root'
#'/store/mc/RunIISpring15MiniAODv2/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/30000/20C7145D-056F-E511-977B-02163E012EFA.root'
#'/store/mc/RunIISpring15MiniAODv2/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/00DA4606-CF6D-E511-AE2C-0025905964A6.root'
        ))

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                  calibratedPatElectrons = cms.PSet(
    initialSeed = cms.untracked.uint32(12345),
    engineName = cms.untracked.string('TRandom3')
    ),
                                                  calibratedPatPhotons = cms.PSet(
    initialSeed = cms.untracked.uint32(12345),
    engineName = cms.untracked.string('TRandom3')
    )
                                                   )

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
process.calibratedPatElectrons.isMC = cms.bool(True)
process.calibratedPatPhotons.isMC = cms.bool(True)
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

# this loads all available b-taggers
process.load("RecoBTag.Configuration.RecoBTag_cff") 
process.load("RecoBTag.SecondaryVertex.pfBoostedDoubleSecondaryVertexAK8BJetTags_cfi")
process.pfImpactParameterTagInfosAK8.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfImpactParameterTagInfosAK8.candidates = cms.InputTag("packedPFCandidates")
process.pfImpactParameterTagInfosAK8.jets = cms.InputTag("slimmedJetsAK8")
process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfosAK8_cfi")
process.pfInclusiveSecondaryVertexFinderTagInfosAK8.extSVCollection = cms.InputTag("slimmedSecondaryVertices")

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#from PhysicsTools.PatAlgos.tools.coreTools import *
#runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
jecLevels = [
  'Summer15_25nsV6_DATA_L2Relative_AK8PFchs.txt',
  'Summer15_25nsV6_DATA_L3Absolute_AK8PFchs.txt'
]

useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
    from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute'] )
    jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, bTagDiscriminators=['pfBoostedDoubleSecondaryVertexAK8BJetTags'],JETCorrPayload='AK8PFchs',JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute']  )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(False)

else :
    dataFormat = DataFormat.MiniAOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
    from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
#    process.load("JetMETCorrections.Type1MET.pfMETmultShiftCorrections_cfi");
#    process.pfMEtMultShiftCorr.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')

process.ggNtuplizer.jecAK8PayloadNames=cms.vstring(jecLevels)
process.ggNtuplizer.runHFElectrons=cms.bool(True)
process.ggNtuplizer.isAOD=cms.bool(useAOD)
process.ggNtuplizer.doGenParticles=cms.bool(True)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(False)

switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.p = cms.Path(
    process.pfImpactParameterTagInfosAK8 *
    process.pfInclusiveSecondaryVertexFinderTagInfosAK8 *
    process.pfBoostedDoubleSecondaryVertexAK8BJetTags 
    * process.calibratedPatElectrons 
    * process.calibratedPatPhotons 
    * process.egmGsfElectronIDSequence 
    * process.egmPhotonIDSequence 
    * process.ggNtuplizer
    )


#print process.dumpPython()
