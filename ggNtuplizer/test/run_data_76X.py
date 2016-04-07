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
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_dataRun2_16Dec2015_v0')

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#jec from sqlite
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                           connect = cms.string('sqlite:Fall15_25nsV2_DATA.db'),
                           toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_DATA_AK4PFchs'),
            label = cms.untracked.string('AK4PFchs')
            ),
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_DATA_AK8PFchs'),
            label = cms.untracked.string('AK8PFchs')
            )))
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/3E38C5C7-87A6-E511-8620-002618943910.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/7EB36900-77A6-E511-ACCB-0CC47A4D764A.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/C67A894C-83A6-E511-A90D-0025905A6068.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FED54E49-83A6-E511-AB65-0CC47A4C8ED8.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/3E65CB18-77A6-E511-A3DE-0025905A612C.root',  
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/7EF48202-89A6-E511-9089-0025905A6110.root',  
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/C68195C6-73A6-E511-A893-0CC47A4C8E34.root',  
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/3E9B66CD-87A6-E511-943C-0025905938A4.root',  
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/802ED605-79A6-E511-A555-0CC47A4D762A.root',  
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/C8E870D0-61A6-E511-8D4C-0CC47A74524E.root', 
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/ECA73F4E-82A6-E511-8843-0CC47A4D76A2.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/EE3246BE-6AA6-E511-A292-002590A2CCFE.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/EEC83EE3-75A6-E511-AE15-0025905A6136.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/EEE889E4-77A6-E511-A7D7-0025905B85D8.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F066634D-83A6-E511-8874-0CC47A4C8E1E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F22451C8-87A6-E511-BEE7-0CC47A4D769A.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F2515AE1-74A6-E511-9524-0CC47A74527A.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F419DFC8-7AA6-E511-BCF7-0CC47A4D769E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F4D15347-66A6-E511-A73F-E41D2D08DE00.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F4F4C3E9-76A6-E511-9572-0025905A6110.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F60AEA20-66A6-E511-BDC4-0CC47A78A3EC.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F67E523E-80A6-E511-836A-0025905964B6.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F687B99A-81A6-E511-9BE4-0CC47A78A440.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F6901A39-6AA6-E511-A634-0025901D4B20.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F6E918C9-87A6-E511-B3D3-0CC47A4D76B2.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F80C3E00-89A6-E511-8B04-0CC47A78A33E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F886D4C6-7AA6-E511-A3DF-00261894398B.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F890A515-5BA6-E511-A1FD-0CC47A4C8E86.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FA5697C8-7AA6-E511-A010-0CC47A78A4B8.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FA78CBE0-75A6-E511-80EB-0025905AA9CC.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FABF4750-7AA6-E511-BBBF-0CC47A4D762A.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FAC03F87-6DA6-E511-90F5-0025901D4A0E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FC122FD8-74A6-E511-ADBA-0025905A6110.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FC6E00D4-78A6-E511-BB82-003048FFD76E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FCAB2A25-57A6-E511-95F9-0CC47A4C8E16.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FCC55B5B-5DA6-E511-9505-0CC47A78A33E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FE092324-64A6-E511-B90A-00266CFFA0B0.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FE2483E5-77A6-E511-A182-0025905A609E.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FEAD8C6B-64A6-E511-B1FC-002618943985.root',
        '/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/FED54E49-83A6-E511-AB65-0CC47A4C8ED8.root'
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
    process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
        payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.patJetsReapplyJEC = process.patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    process.patJetAK8CorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJetsAK8"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'],
        payload = 'AK8PFchs' ) # Make sure to choose the appropriate levels and payload here!

    process.patJetsAK8ReapplyJEC = process.patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJetsAK8"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetAK8CorrFactorsReapplyJEC"))
        )

    process.reapplyJEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC +process.patJetAK8CorrFactorsReapplyJEC + process. patJetsAK8ReapplyJEC )

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
        process.reapplyJEC*
        process.pfImpactParameterTagInfosAK8 *
        process.pfInclusiveSecondaryVertexFinderTagInfosAK8 *
        process.pfBoostedDoubleSecondaryVertexAK8BJetTags *
        process.calibratedPatElectrons*
        process.calibratedPatPhotons*
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.ggNtuplizer
        )
    
#print process.dumpPython()
