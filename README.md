#### Current production tag : V08_00_26_03
#### Newest tag for testing : 
#### Note that the current head version can be run with CMSSW_8_0_26_patch1

##### To work with CMSSW_8_0_26_patch1 and head version, you do :
cd CMSSW_8_0_26_patch1/src <br>
cmsenv <br>
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily <br>
git cms-init <br>
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git <br>
git fetch --tags btv-cmssw <br>
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21 <br>
git cms-merge-topic rafaellopesdesa:RegressionCheckNegEnergy <br>
git cms-merge-topic cms-egamma:EGM_gain_v1 <br>
cd EgammaAnalysis/ElectronTools/data <br>
git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git <br>
cd $CMSSW_BASE/src <br>
git cms-merge-topic cms-met:METRecipe_8020 -u <br>
git cms-merge-topic cms-met:METRecipe_80X_part2 -u <br>
git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout <br>
git cms-merge-topic Sam-Harper:PackedCandNoPuppi <br>
git cms-merge-topic ikrav:egm_id_80X_v2 <br>
git cms-merge-topic ikrav:egm_id_80X_v3_photons <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V08_00_26_03 https://github.com/cmkuo/ggAnalysis.git <br>

scram b -j 10 <br>

##### To work with CMSSW_8_0_26_patch1 and V08_00_26_01, you do :
cd CMSSW_8_0_26_patch1/src <br>
cmsenv <br>
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily <br>
git cms-init <br>
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git <br>
git fetch --tags btv-cmssw <br>
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21 <br>
git cms-merge-topic rafaellopesdesa:EgammaAnalysis80_EGMSmearer_Moriond17_23Jan <br>
cd EgammaAnalysis/ElectronTools/data <br>
git clone git@github.com:ECALELFS/ScalesSmearings.git <br>
cd ../../../ <br>
git cms-merge-topic cms-met:METRecipe_8020 -u <br>
git cms-merge-topic cms-met:METRecipe_80X_part2 -u <br>
git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout <br>
git cms-merge-topic Sam-Harper:PackedCandNoPuppi <br>
git cms-merge-topic ikrav:egm_id_80X_v2 <br>
git cms-merge-topic ikrav:egm_id_80X_v3_photons <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V08_00_26_01 https://github.com/cmkuo/ggAnalysis.git <br>

scram b -j 10 <br>

##### To work with CMSSW_7_6_3_patch2, you do:
cd CMSSW_7_6_3_patch2/src <br>
cmsenv <br>
git cms-merge-topic -u matteosan1:smearer_76X <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V07_06_03_01 https://github.com/cmkuo/ggAnalysis.git <br>
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

