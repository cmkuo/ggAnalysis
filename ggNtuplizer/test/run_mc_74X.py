import FWCore.ParameterSet.Config as cms

process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v3')
process.load("Configuration.StandardSequences.MagneticField_cff")

#process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/02CDE51A-726D-E511-B2C4-0025905C96EA.root'
        #'/store/mc/RunIISpring15DR74/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/0E885FA4-FF01-E511-B2DB-002590A3C95E.root'
        #'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/00000/04051FA7-CB05-E511-91BE-00266CFAE264.root'
        #'file:/data4/cmkuo/testfiles/GJet_pt_15to6000_13TeV_Spring15_Asympt50ns_MCRUN2_74_V9A-v3.root'
        #'file:/data4/cmkuo/testfiles/DYJetsToLL_M50_13TeV_Spring15_Asympt50ns_MCRUN2_74_V9A-v3.root'
        #'file:/data4/cmkuo/testfiles/WJetsToLNu_13TeV_Spring15_Asympt25ns_MCRUN2_74_V9-v1.root'
        #'file:/data4/cmkuo/testfiles/TTJets_amcatnloFXFX-pythia813TeV_Asympt25ns_MCRUN2_74_V9-v1.root'
        #'file:/data4/cmkuo/testfiles/WZ_13TeV_Spring15_Asympt25ns_MCRUN2_74_V9-v1.root'
        #'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/00000/04051FA7-CB05-E511-91BE-00266CFAE264.root'
#        '/store/mc/RunIISpring15DR74/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v2/70000/04AFB8D8-900C-E511-8FA1-3417EBE6471A.root'
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

useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_cfi")
    from JMEAnalysis.JetToolbox.jetToolbox_cff import *
    jetToolbox( process, 'ak4', 'ak4PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, addPUJetID=True, JETCorrPayload='AK4PFchs', JETCorrLevels=['L1FastJet','L2Relative', 'L3Absolute'] )
    jetToolbox( process, 'ak8', 'ak8PFJetsCHS', 'out', miniAOD= False, addSoftDrop=True, addSoftDropSubjets=True, addNsub=True, bTagDiscriminators=['pfBoostedDoubleSecondaryVertexAK8BJetTags'] )
    process.ggNtuplizer.dumpSoftDrop= cms.bool(False)

else :
    dataFormat = DataFormat.MiniAOD
    process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
    process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
#    process.load("JetMETCorrections.Type1MET.pfMETmultShiftCorrections_cfi");
#    process.pfMEtMultShiftCorr.vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')

process.ggNtuplizer.runHFElectrons=cms.bool(True)
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
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff',
                    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.p = cms.Path(
    process.egmGsfElectronIDSequence 
    * process.egmPhotonIDSequence 
    * process.ggNtuplizer
    )


#print process.dumpPython()
