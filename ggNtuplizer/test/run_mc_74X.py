import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

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
        #'file:DYJetsToLL_M50_13TeV_Spring15_Asympt50ns_MCRUN2_74_V9A-v3.root'
'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/00000/04051FA7-CB05-E511-91BE-00266CFAE264.root'
        #'file:/data4/cmkuo/testfiles/GJet_pt_15to6000_13TeV_Spring15_Asympt50ns_MCRUN2_74_V9A-v3.root'
        #'file:/data4/cmkuo/testfiles/DYJetsToLL_M50_13TeV_Spring15_Asympt50ns_MCRUN2_74_V9A-v3.root'
        #'file:/data4/cmkuo/testfiles/WJetsToLNu_13TeV_Spring15_Asympt25ns_MCRUN2_74_V9-v1.root'
        #'file:/data4/cmkuo/testfiles/TTJets_amcatnloFXFX-pythia813TeV_Asympt25ns_MCRUN2_74_V9-v1.root'
        #'file:/data4/cmkuo/testfiles/WZ_13TeV_Spring15_Asympt25ns_MCRUN2_74_V9-v1.root'
        #'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/00000/04051FA7-CB05-E511-91BE-00266CFAE264.root'
       #'/store/mc/RunIISpring15DR74/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v2/70000/04AFB8D8-900C-E511-8FA1-3417EBE6471A.root'
        ))

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

#from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
#from PhysicsTools.PatAlgos.tools.coreTools import *
#runOnData( process, outputModules = [] )
#removeMCMatching(process, names=['All'], outputModules=[])

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_mc.root'))

#####VID framework####################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 

useAOD = True

if useAOD == True :
    dataFormat = DataFormat.AOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
    process.ggNtuplizer.runHFElectrons=cms.bool(True)
    doNoHFMet = False
    from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute'] )
    jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(False)

else :
    dataFormat = DataFormat.MiniAOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
    doNoHFMet = True
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)

process.ggNtuplizer.isAOD=cms.bool(useAOD)
process.ggNtuplizer.doGenParticles=cms.bool(True)
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
tellMETData = False

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
