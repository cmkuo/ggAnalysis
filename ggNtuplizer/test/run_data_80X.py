import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff" )
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v5')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/data/Run2016B/HLTPhysics/MINIAOD/PromptReco-v1/000/272/022/00000/98263424-E20F-E611-864C-02163E011D0A.root'
        )
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015'),
    electrons = cms.InputTag("slimmedElectrons"),
    gbrForestName = cms.string('gedelectron_p4combination_25ns'),
    isMC = cms.bool(False),
    isSynchronization = cms.bool(False)
)

process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015'),
    isMC = cms.bool(False),
    isSynchronization = cms.bool(False),
    photons = cms.InputTag("slimmedPhotons")
)

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])

# this loads all available b-taggers
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoBTag.SecondaryVertex.pfBoostedDoubleSecondaryVertexAK8BJetTags_cfi")
process.pfImpactParameterTagInfosAK8.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfImpactParameterTagInfosAK8.candidates = cms.InputTag("packedPFCandidates")
process.pfImpactParameterTagInfosAK8.jets = cms.InputTag("slimmedJetsAK8")
process.load("RecoBTag.SecondaryVertex.pfInclusiveSecondaryVertexFinderTagInfosAK8_cfi")
process.pfInclusiveSecondaryVertexFinderTagInfosAK8.extSVCollection = cms.InputTag("slimmedSecondaryVertices")

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

jecLevels = [
  'Fall15_25nsV2_DATA_L2Relative_AK8PFchs.txt',
  'Fall15_25nsV2_DATA_L3Absolute_AK8PFchs.txt',
  'Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt'
]

useAOD = False
#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
    process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
    process.ggNtuplizer.addFilterInfo=cms.bool(True)
    process.ggNtuplizer.development=cms.bool(False)
    from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', runOnMC = False, miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'] )
    jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', runOnMC = False, miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, bTagDiscriminators=['pfBoostedDoubleSecondaryVertexAK8BJetTags'] )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(False)

else :
    dataFormat = DataFormat.MiniAOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
    process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
    process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
    ###process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    ###from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    ###process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        ###src = cms.InputTag("slimmedJets"),
        ###levels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
        ###payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

    ###from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    ###process.patJetsReapplyJEC = process.patJetsUpdated.clone(
        ###jetSource = cms.InputTag("slimmedJets"),
        ###jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        ###)

    ###process.patJetAK8CorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        ###src = cms.InputTag("slimmedJetsAK8"),
        ###levels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
        ###payload = 'AK8PFchs' ) # Make sure to choose the appropriate levels and payload here!

    ###process.patJetsAK8ReapplyJEC = process.patJetsUpdated.clone(
        ###jetSource = cms.InputTag("slimmedJetsAK8"),
        ###jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetAK8CorrFactorsReapplyJEC"))
        ###)

    ###process.reapplyJEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC +process.patJetAK8CorrFactorsReapplyJEC + process. patJetsAK8ReapplyJEC )

process.ggNtuplizer.jecAK8PayloadNames=cms.vstring(jecLevels)
process.ggNtuplizer.runHFElectrons=cms.bool(True)
process.ggNtuplizer.isAOD=cms.bool(useAOD)
process.ggNtuplizer.doGenParticles=cms.bool(False)
process.ggNtuplizer.dumpSubJets=cms.bool(True)
process.ggNtuplizer.dumpJets=cms.bool(True)
process.ggNtuplizer.dumpTaus=cms.bool(True)

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

if useAOD == True:
    process.p = cms.Path(
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.ggMETFiltersSequence* 
        process.ggNtuplizer
        )
else:
    process.p = cms.Path(
        process.HBHENoiseFilterResultProducer* # produces HBHE bools
#        process.ApplyBaselineHBHENoiseFilter*  # reject events 
        ###process.reapplyJEC*
        ###process.pfImpactParameterTagInfosAK8 *
        ###process.pfInclusiveSecondaryVertexFinderTagInfosAK8 *
        ###process.pfBoostedDoubleSecondaryVertexAK8BJetTags *
        process.calibratedPatElectrons*
        process.calibratedPatPhotons*
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.ggNtuplizer
        )
    
#print process.dumpPython()
