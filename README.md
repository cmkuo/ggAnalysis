##### To work with CMSSW_7_4_1, you do:

cd CMSSW_7_4_1/src <br>
cmsenv <br>
git cms-merge-topic ikrav:egm_id_74X_v0 <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone https://github.com/cmkuo/ggAnalysis.git <br>
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

##### To work with CMSSW_7_2_0 or 7_2_3, you do:

git cms-merge-topic ikrav:egm_id_phys14 # for photon ID recipe <br>
git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720 # for electron ID recipe <br>
git clone https://github.com/cmkuo/HiggsAnalysis.git <br>
git clone -b V07-02-03-00 https://github.com/cmkuo/ggAnalysis.git <br>


