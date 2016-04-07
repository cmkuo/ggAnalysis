import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_RunIIFall15DR76_v1')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.Tracer = cms.Service("Tracer")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#jec from sqlite
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
  connect = cms.string('sqlite:Fall15_25nsV2_MC.db'),
  toGet = cms.VPSet(
  cms.PSet(
   record = cms.string('JetCorrectionsRecord'),
   tag = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'),
   label = cms.untracked.string('AK4PFchs')
   ),
  cms.PSet(
   record = cms.string('JetCorrectionsRecord'),
   tag = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK8PFchs'),
   label = cms.untracked.string('AK8PFchs')
   )))
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/ACA927E2-04C8-E511-971A-E41D2D08DDF0.root'
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/E010AD0D-27B9-E511-B6BC-3417EBE705EB.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/E02A181E-28B9-E511-BC2B-003048D4794E.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/E467BD62-20B9-E511-916D-3417EBE6FFFD.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/E876A49F-20B9-E511-A547-002590DE39F0.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/ECBC1E37-F7B8-E511-A35E-002590DE6E30.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/EEB51ECC-25B9-E511-9AD2-002590DE6DE4.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/EEE5367B-24B9-E511-87BB-3417EBE70729.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/F0AAB51B-25B9-E511-94EE-002590494C88.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/F600345C-05B9-E511-B7C2-3417EBE70729.root',
        '/store/mc/RunIIFall15MiniAODv2/GJet_Pt-15to6000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOff_76X_mcRun2_asymptotic_v12-v1/70000/F61F065C-05B9-E511-95D1-3417EBE743C0.root'
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

#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.load('EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi')
#process.calibratedPatElectrons.isMC = cms.bool(True)
#process.calibratedPatPhotons.isMC = cms.bool(True)
process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015'),
    electrons = cms.InputTag("slimmedElectrons"),
    gbrForestName = cms.string('gedelectron_p4combination_25ns'),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(False)
)

process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015'),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(False),
    photons = cms.InputTag("slimmedPhotons")
)

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

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
  'Fall15_25nsV2_MC_L2Relative_AK8PFchs.txt',
  'Fall15_25nsV2_MC_L3Absolute_AK8PFchs.txt'
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
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
    process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
    process.patJetCorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
        payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
    process.patJetsReapplyJEC = process.patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    process.patJetAK8CorrFactorsReapplyJEC = process.patJetCorrFactorsUpdated.clone(
        src = cms.InputTag("slimmedJetsAK8"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
        payload = 'AK8PFchs' ) # Make sure to choose the appropriate levels and payload here!

    process.patJetsAK8ReapplyJEC = process.patJetsUpdated.clone(
        jetSource = cms.InputTag("slimmedJetsAK8"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetAK8CorrFactorsReapplyJEC"))
        )

    process.reapplyJEC = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC +process.patJetAK8CorrFactorsReapplyJEC + process. patJetsAK8ReapplyJEC )

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
    process.reapplyJEC*
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
