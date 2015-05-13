To work with CMSSW_7_2_0 or 7_2_3, you do:

git cms-merge-topic ikrav:egm_id_phys14 # for photon ID recipe <br>
git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720 # for electron ID recipe <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone https://github.com/cmkuo/ggAnalysis.git <br>

####FOR VID Framework this is needed
cd CMSSW_7_4_X/src
cmsenv
git cms-merge-topic ikrav:egm_id_74X_v0
scram b -j 10
