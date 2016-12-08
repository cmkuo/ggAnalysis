#### Current production tag : V08_00_20_00
#### Newest tag for testing : master
#### Note that the current head version can be run with CMSSW_8_0_11

##### To work with CMSSW_8_0_20, you do :
cd CMSSW_8_0_20/src <br>
cmsenv <br>
git cms-init <br>
git cms-merge-topic emanueledimarco:ecal_smear_fix_80X <br>
cd EgammaAnalysis/ElectronTools/data <br>
git clone -b ICHEP2016_v2 https://github.com/ECALELFS/ScalesSmearings.git <br>
cd ../../../ <br>
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter <br>
git cms-merge-topic cms-met:METRecipe_8020 <br>
git cms-merge-topic ikrav:egm_id_80X_v2 <br>
cd $CMSSW_BASE/external <br>
cd slc6_amd64_gcc530/ <br>
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data <br>
cd data/RecoEgamma/ElectronIdentification/data <br>
git checkout egm_id_80X_v1 <br>
cd $CMSSW_BASE/src <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V08_00_20_00 https://github.com/cmkuo/ggAnalysis.git <br>
scram b -j 10 <br>

##### To work with CMSSW_8_0_11 and V08_00_11_01, you do :
cd CMSSW_8_0_11/src <br>
cmsenv <br>
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily <br>
git cms-init <br>
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git <br>
git fetch --tags btv-cmssw <br>
git cms-merge-topic cms-btv-pog:BoostedDoubleSVTaggerV3-WithWeightFiles-v1_from-CMSSW_8_0_8_patch1 <br>
git remote add -f -t ecal_smear_fix_80X emanueledimarco https://github.com/emanueledimarco/cmssw.git <br>
git cms-addpkg EgammaAnalysis/ElectronTools <br>
git checkout -b from-52f192a 52f192a <br>
cd EgammaAnalysis/ElectronTools/data <br>
git clone -b ICHEP2016_v2 https://github.com/ECALELFS/ScalesSmearings.git <br>
cd ../../../ <br>
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V08_00_11_01 https://github.com/cmkuo/ggAnalysis.git <br>
scram b -j 10 <br>

##### To work with CMSSW_8_0_11 and V08_00_11_00, you do :
cd CMSSW_8_0_11/src <br>
cmsenv <br>
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily <br>
git cms-init <br>
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git <br>
git fetch --tags btv-cmssw <br>
git cms-merge-topic cms-btv-pog:BoostedDoubleSVTaggerV3-WithWeightFiles-v1_from-CMSSW_8_0_8_patch1 <br>
git remote add -f -t ecal_smear_fix_80X emanueledimarco https://github.com/emanueledimarco/cmssw.git <br>
git cms-addpkg EgammaAnalysis/ElectronTools <br>
git checkout -b from-277de3c 277de3c <br>
cd EgammaAnalysis/ElectronTools/data <br>
git clone -b ICHEP2016_approval_7p65fb https://github.com/emanueledimarco/ScalesSmearings.git <br>
cd ../../../ <br>
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V08_00_11_00 https://github.com/cmkuo/ggAnalysis.git <br>
scram b -j 10 <br>

##### To work with CMSSW_8_0_10, you do:
cd CMSSW_8_0_10/src <br>
cmsenv <br>
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily <br>
git cms-init <br>
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git <br>
git fetch --tags btv-cmssw <br>
git cms-merge-topic cms-btv-pog:BoostedDoubleSVTaggerV3-WithWeightFiles-v1_from-CMSSW_8_0_8_patch1 <br>
git cms-merge-topic -u matteosan1:smearer_76X <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V08_00_10_00 https://github.com/cmkuo/ggAnalysis.git <br>
scram b -j 10 <br>

##### To work with CMSSW_7_6_3_patch2, you do:
cd CMSSW_7_6_3_patch2/src <br>
cmsenv <br>
git cms-merge-topic -u matteosan1:smearer_76X <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V07_06_03_01 https://github.com/cmkuo/ggAnalysis.git <br>
scram b -j 10 <br>

##### To work with CMSSW_7_4_16, you do:
cd CMSSW_7_4_16/src <br>
cmsenv <br>
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily <br>
git cms-init <br>
git cms-merge-topic matteosan1:smearer <br>
git clone https://github.com/cms-jet/JetToolbox JMEAnalysis/JetToolbox <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V07-04-16-02 https://github.com/cmkuo/ggAnalysis.git <br>
scram b -j 10 <br>

The above code stores the decision in 64 integer. Each bit represents a decision<br>
for ELECRON ID: 5 IDs (Veto, Loose, Medium, Tight and HEEP) so only 5 bits are imp for us (59 bits of this integer  we are not using so may be we can change that to 16 bit integer later)<br>
Representing that integer in 5 bits: b4 b3 b2 b1 b0<br>
b0: Veto; b1: Loose; b2: Medium; b3: Tight and b4: HEEP<br>
To access the decision for <br>
(a) veto: eleIDbit[]>>0&1 ---> gives 0 or 1. if 0--> this eID is failed. if 1--> this eID is passed<br>
(b) Loose: eleIDbit[]>>1&1<br>
(c) Medium: eleIDbit[]>>2&1<br>
(d) Tight: eleIDbit[]>>3&1<br>
(e) HEEP: eleIDbit[]>>4&1<br>

for photons it is done the same way: it has 3 IDs<br>
so 3 bits represent the decision<br>
Representing that integer in 3 bits:  b2 b1 b0<br>
b0: Loose; b1: Medium; b2: Tight<br>
To access the decision for <br>
(a) Loose: phoIDbit[]>>0&1 ---> gives 0 or 1. if 0--> this phoID is failed. if 1--> this phoID is passed<br>
(b) Medium: phoIDbit[]>>1&1<br>
(c) Tight: phoIDbit[]>>2&1<br>

to access the MC status flag with GEN particles <br>
(a) fromHardProcessFinalState : mcStatusFlag[]>>0&1 ---> gives 0 (no) or 1 (yes). <br>
(b) isPromptFinalState        : mcStatusFlag[]>>1&1 ---> gives 0 (no) or 1 (yes). <br>
(c) fromHardProcessBeforeFSR  : mcStatusFlag[]>>2&1 ---> gives 0 (no) or 1 (yes). <br>

