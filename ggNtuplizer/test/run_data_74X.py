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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

#jec from sqlite
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
        connect = cms.string('sqlite:DBFiles/Summer15_50nsV4_MC.db'),
        toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV4_MC_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            )))
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/data/Run2015B/SinglePhoton/AOD/PromptReco-v1/000/251/562/00000/000B31EB-262A-E511-8D03-02163E011976.root'
        #'/store/data/Run2015B/SinglePhoton/MINIAOD/17Jul2015-v1/30000/9AA235B1-C12E-E511-91BB-002618943902.root'
        '/store/data/Run2015B/DoubleEG/AOD/PromptReco-v1/000/251/244/00000/FA8EC6E3-D528-E511-B1D5-02163E012BD2.root'
#        '/store/express/Run2015C/ExpressPhysics/FEVT/Express-v1/000/254/879/00000/FA465069-4D49-E511-AD11-02163E011E1E.root'
        #'file:DoubleEG_Run2015B_251562_miniAOD.root'
        #'/store/data/Run2015B/SinglePhoton/AOD/PromptReco-v1/000/251/562/00000/000B31EB-262A-E511-8D03-02163E011976.root'
        #'/store/data/Run2015B/DoubleEG/AOD/PromptReco-v1/000/251/244/00000/FA8EC6E3-D528-E511-B1D5-02163E012BD2.root'
        #'/store/data/Run2015B/DoubleEG/MINIAOD/PromptReco-v1/000/251/562/00000/CCEFEFB1-402A-E511-BF33-02163E0136E2.root'
)
                            )

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))

useAOD = True
#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
    process.load("ggAnalysis.ggNtuplizer.ggMETFilters_cff")
    process.ggNtuplizer.addFilterInfo=cms.bool(True)
    process.ggNtuplizer.runHFElectrons=cms.bool(True)
    doNoHFMet = False
    from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'] )
    jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(False)

else :
    dataFormat = DataFormat.MiniAOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
    doNoHFMet = False
    process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")
    process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
    process.ApplyBaselineHBHENoiseFilter = cms.EDFilter("BooleanFlagFilter",
                                                        inputLabel = cms.InputTag("HBHENoiseFilterResultProducer",
                                                                                  "HBHENoiseFilterResult"),
                                                        reverseDecision = cms.bool(False)
                                                       )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)

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
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.ggNtuplizer.doNoHFMET=cms.bool(doNoHFMet)
tellMETData = True

if doNoHFMet == True:
    process.noHFCands = cms.EDFilter("CandPtrSelector",
                                     src=cms.InputTag("packedPFCandidates"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )
    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
    #default configuration for miniAOD reprocessing, change the isData flag to run on data
    #for a full met computation, remove the pfCandColl input
    runMetCorAndUncFromMiniAOD(process,
                               isData=tellMETData,
                               pfCandColl=cms.InputTag("noHFCands"),
                               postfix="NoHF"
                               )
    process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")


if useAOD == True:
    process.p = cms.Path(
        ###process.egmGsfElectronIDSequence
        # process.mvaTrigV050nsCSA14
        # + process.mvaTrigV025nsCSA14
        # + process.mvaNonTrigV050nsCSA14
        # + process.mvaNonTrigV025nsCSA14
        # + process.mvaNonTrigV025nsPHYS14 
        #    process.patDefaultSequence *
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.ggMETFiltersSequence* 
        process.ggNtuplizer
        )
else:
    process.p = cms.Path(
        ###process.egmGsfElectronIDSequence
        # process.mvaTrigV050nsCSA14
        # + process.mvaTrigV025nsCSA14
        # + process.mvaNonTrigV050nsCSA14
        # + process.mvaNonTrigV025nsCSA14
        # + process.mvaNonTrigV025nsPHYS14 
        #    process.patDefaultSequence *
        process.HBHENoiseFilterResultProducer* # produces HBHE bools
#        process.ApplyBaselineHBHENoiseFilter*  # reject events 
        process.egmGsfElectronIDSequence*
        process.egmPhotonIDSequence*
        process.ggNtuplizer
        )
    

#print process.dumpPython()
