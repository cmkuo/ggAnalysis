if __name__ == '__main__':

# Usage : python crabConfig.py (to create jobs)
#         ./multicrab -c status -d <work area> (to check job status)

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    from CRABClient.UserUtilities import config
    config = config()
    
    from multiprocessing import Process

    # Common configuration

    config.General.workArea     = 'crab_projects_ntuples_25ns'
    config.General.transferLogs = False
    config.JobType.pluginName   = 'Analysis' # PrivateMC
    config.JobType.psetName     = 'run_mc_74X.py'
    config.Data.inputDBS        = 'global'    
    config.Data.splitting       = 'LumiBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
    config.Data.totalUnits      = -1
    config.Data.publication     = False
    config.Site.storageSite     = 'T2_CH_CERN'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers

    # dataset dependent configuration

    config.General.requestName = 'job_spring15_gjet_pt20to40_MGG_80toInf_25ns'
    config.Data.unitsPerJob    = 40
    config.Data.inputDataset   = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_gjet_pt20to40_MGG_80toInf_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_gjet_pt40_MGG_80toInf_25ns'
    config.Data.unitsPerJob    = 60
    config.Data.inputDataset   = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_gjet_pt40_MGG_80toInf_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_gjet_pt20to40_MGG_40to80_25ns'
    config.Data.unitsPerJob    = 40
    config.Data.inputDataset   = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_gjet_pt20to40_MGG_40to80_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_gjet_pt15to6000_25ns'
    config.Data.unitsPerJob    = 50
    config.Data.inputDataset   = '/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_gjet_pt15to6000_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_qcd_pt30to40_MGG_80toInf_25ns'
    config.Data.unitsPerJob    = 60
    config.Data.inputDataset   = '/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_qcd_pt30to40_MGG_80toInf_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_qcd_pt40_MGG_80toInf_25ns'
    config.Data.unitsPerJob    = 40
    config.Data.inputDataset   = '/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_qcd_pt40_MGG_80toInf_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_qcd_pt30_MGG_40to80_25ns'
    config.Data.unitsPerJob    = 60
    config.Data.inputDataset   = '/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_qcd_pt30_MGG_40to80_25ns'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'job_spring15_DYJetsToLL_m50'
    config.Data.unitsPerJob    = 150
    config.Data.inputDataset   = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM'
    config.Data.outLFNDirBase  = '/store/user/cmkuo/job_spring15_DYJetsToLL_m50_25ns'
    p = Process(target=submit, args=(config,))
    p.start()        
    p.join()    


