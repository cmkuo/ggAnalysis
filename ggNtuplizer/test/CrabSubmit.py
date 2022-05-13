from CRABClient.UserUtilities import config
import sys
config = config()

submitVersion = "EgammaID_ntuples_RunIII"

mainOutputDir = '/store/group/phys_egamma/prrout/Egamma_Photon_Identification/RunIII/ggNtuples/%s' % submitVersion

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName  = '/afs/cern.ch/work/p/prrout/public/Egamma_RunIII_Ntuplizer_CMKuo/CMSSW_12_2_1/src/ggAnalysis/ggNtuplizer/test/EgammaID_RunIII_ConfFile_cfg.py'

config.JobType.sendExternalFolder = True
config.JobType.allowUndistributedCMSSW = True

config.Data.allowNonValidInputDataset = True

config.Data.inputDBS = 'global'
config.Data.publication = False

#config.Data.publishDataName = 
config.Site.storageSite = 'T2_CH_CERN'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_%s' % submitVersion

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    ##### now submit DATA
    config.Data.outLFNDirBase = '%s/%s/' % (mainOutputDir,'data')
    config.Data.splitting     = 'LumiBased'
    config.Data.totalUnits      = -1
#    config.Data.lumiMask      = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt' #UL2018  
    config.Data.unitsPerJob   = 40
#    config.Data.unitsPerJob   = 200

    ### For different DataSets
#    config.General.requestName  = 'EGamma_UL2018A'
#    config.Data.inputDataset    = '/EGamma/Run2018A-12Nov2019_UL2018-v2/AOD'
#    submit(config)

    config.General.requestName = 'job_GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_Run3Winter22MiniAoD_2022'
    config.Data.inputDataset   = '/GJet_Pt-10to40_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/Run3Winter22MiniAOD-FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/MINIAODSIM'
    submit(config)    

    config.General.requestName = 'job_GJet_Pt-40toInf_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8_Run3Winter22MiniAoD_2022'
    config.Data.inputDataset   = '/GJet_Pt-40toInf_DoubleEMEnriched_TuneCP5_13p6TeV_pythia8/Run3Winter22MiniAOD-FlatPU0to70_122X_mcRun3_2021_realistic_v9-v2/MINIAODSIM'
    submit(config)    

    #config.General.requestName = 'job_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2021'
    #config.Data.inputDataset   = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2021Scenario_106X_mcRun3_2021_realistic_v3-v1/MINIAODSIM'
    #submit(config)    

    #config.General.requestName = 'job_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2023'
    #config.Data.inputDataset   = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v1/MINIAODSIM'
    #submit(config)    

    #config.General.requestName = 'job_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_Run3Summer19_2023'
    #config.Data.inputDataset   = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v2/MINIAODSIM'
    #submit(config)    

    #config.General.requestName = 'job_GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_Run3Summer19_2023'
    #config.Data.inputDataset   = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_14TeV_Pythia8/Run3Summer19MiniAOD-2023Scenario_106X_mcRun3_2023_realistic_v3-v2/MINIAODSIM'
    #submit(config)    

