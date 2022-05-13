# Ntuple for Egamma electron and photon identification in CMS
# This ggtuplizer is mainly developed by cmkuo and contributed by A. Roy (https://github.com/ashimroy/EGMObjectDumper), A. Purohit (https://github.com/ArnabPurohit/PhotonDemoAnalyzer) and others from time to time. 
# Key points in CMSSW_12_X: Some changes from accessing "ECalClusterLaZyTools" was made from CMSSW_11_X. These changes are included in this branch 12_X.  

cmsrel CMSSW_12_2_1

cd CMSSW_12_2_1/src

cmsenv

git clone -b 12_X https://github.com/cmkuo/ggAnalysis

scram b -j8

Interactive run:
cd ggAnalysis/ggNtuplizer/test

# set up grid environment.

voms-proxy-init --voms cms --valid 192:00

Check the config file

cmsRun EgammaID_RunIII_ConfFile_cfg.py

For CRAB jobs:

cd ggAnalysis/ggNtuplizer/test

source set_env.sh ## It will ask for grid password.
python CrabSubmit.py ## You can edit this file on your requirement
