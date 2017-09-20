#include "map"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// local variables: per-filter per-electron/muon/photon/jet arrays of matched trigger objects
// NOTE: number of elements in the arrays equals sizeof(ULong64_t)
vector<float> trgSingleElePt[64], trgSingleEleEta[64], trgSingleElePhi[64];
vector<float> trgDoubleElePt[64], trgDoubleEleEta[64], trgDoubleElePhi[64];
vector<float> trgSinglePhoPt[64], trgSinglePhoEta[64], trgSinglePhoPhi[64];
vector<float> trgDoublePhoPt[64], trgDoublePhoEta[64], trgDoublePhoPhi[64];
vector<float> trgTriplePhoPt[64], trgTriplePhoEta[64], trgTriplePhoPhi[64];
vector<float> trgMuPt[64],  trgMuEta[64],  trgMuPhi[64];
vector<float> trgJetPt[64], trgJetEta[64], trgJetPhi[64];
vector<float> trgL1Eta[64],  trgL1Phi[64];

void ggNtuplizer::initTriggerFilters(const edm::Event &e) {
  // Fills the arrays above.

  // cleanup from previous execution
  for (size_t i = 0; i < 64; ++i) {
    trgSingleElePt [i].clear();
    trgSingleEleEta[i].clear();
    trgSingleElePhi[i].clear();
    trgDoubleElePt [i].clear();
    trgDoubleEleEta[i].clear();
    trgDoubleElePhi[i].clear();
    trgSinglePhoPt [i].clear();
    trgSinglePhoEta[i].clear();
    trgSinglePhoPhi[i].clear();
    trgDoublePhoPt [i].clear();
    trgDoublePhoEta[i].clear();
    trgDoublePhoPhi[i].clear();
    trgTriplePhoPt [i].clear();
    trgTriplePhoEta[i].clear();
    trgTriplePhoPhi[i].clear();
    trgMuPt  [i].clear();
    trgMuEta [i].clear();
    trgMuPhi [i].clear();
    trgJetPt [i].clear();
    trgJetEta[i].clear();
    trgJetPhi[i].clear();
    trgL1Eta[i].clear();
    trgL1Phi[i].clear();
  }

  // filter => index (in trg*[] arrays) mappings
  static std::map<string,size_t> eleSingleFilters;
  static std::map<string,size_t> eleDoubleFilters;
  static std::map<string,size_t> phoSingleFilters;
  static std::map<string,size_t> phoDoubleFilters;
  static std::map<string,size_t> phoTripleFilters;
  static std::map<string,size_t> muFilters;
  static std::map<string,size_t> jetFilters;
  static std::map<string,size_t> l1Filters;

  // one-time initialization
  if (eleSingleFilters.size() == 0) {
    eleSingleFilters["hltEle27WPLooseGsfTrackIsoFilter"] = 0;
    eleSingleFilters["hltEle25WP60SC4HcalIsoFilter"] = 1;
    // for HLT_Ele23_WPLoose_Gsf_v and EventTree
    eleSingleFilters["hltEGL1SingleEG40ORSingleIsoEG22erOrSingleIsoEG24erORSingleIsoEG24OrSingleIsoEG26Filter"] = 2;
    //HLT_Ele22_eta2p1_WPLoose_Gsf_v3
    eleSingleFilters["hltSingleEle22WPLooseGsfTrackIsoFilter"] = 3;
    //HLT_Ele24_eta2p1_WPLoose_Gsf_v1
    eleSingleFilters["hltSingleEle24WPLooseGsfTrackIsoFilter"] = 4;
    //HLT_Ele25_WPTight_Gsf_v1
    eleSingleFilters["hltEle25WPTightGsfTrackIsoFilter"] = 5;
    //HLT_Ele25_eta2p1_WPLoose_Gsf_v1
    eleSingleFilters["hltEle25erWPLooseGsfTrackIsoFilter"] = 6;
    //HLT_Ele25_eta2p1_WPTight_Gsf_v1
    eleSingleFilters["hltEle25erWPTightGsfTrackIsoFilter"] = 7;
    //HLT_Ele27_WPLoose_Gsf_v1
    eleSingleFilters["hltEle27noerWPLooseGsfTrackIsoFilter"] = 8;
    //HLT_Ele27_WPLoose_Gsf_v
    eleSingleFilters["hltEG27L1IsoEG22erORIsoEG24erORIsoEG24ORIsoEG26OREG40EtFilter"] = 9;
    //HLT_Ele27_eta2p1_WPLoose_Gsf_v2, HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau*                             
    eleSingleFilters["hltEle27erWPLooseGsfTrackIsoFilter"] = 10;
    //HLT_Ele27_eta2p1_WPTight_Gsf                                                                             
    eleSingleFilters["hltEle27erWPTightGsfTrackIsoFilter"] = 11;
    //HLT_Ele27_WPTight_Gsf                                                                                    
    eleSingleFilters["hltEle27WPTightGsfTrackIsoFilter"] = 12;
    //HLT_Ele32_eta2p1_WPTight_Gsf                                                                             
    eleSingleFilters["hltEle32WPTightGsfTrackIsoFilter"] = 13;
    //HLT_Ele35_WPLoose_Gsf                                                                                    
    eleSingleFilters["hltEle35WPLooseGsfTrackIsoFilter"] = 14;
    //HLT_Ele45_WPLoose_Gsf                                                                                    
    eleSingleFilters["hltEle45WPLooseGsfTrackIsoFilter"] = 15;
    //HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1 
    eleSingleFilters["hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter"] = 16;
    //HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1                                                    
    eleSingleFilters["hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20"] = 17;
    //HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1                                                    
    eleSingleFilters["hltEle24WPLooseL1SingleIsoEG22erGsfTrackIsoFilter"] = 18;
    //HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20                                                             
    eleSingleFilters["hltEle24WPLooseL1IsoEG22erTau20erGsfTrackIsoFilter"] = 19;
    //HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1                                                          
    eleSingleFilters["hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20"] = 20;
    //HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90                                                   
    eleSingleFilters["hltEG18Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter"] = 21;
    eleSingleFilters["hltEG18R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter"] = 22; 
    //HLT_DiMu9_Ele9_CaloIdL_TrackIdL                                                                     
    eleSingleFilters["hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter"] = 23; 
    //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL                                                    
    eleSingleFilters["hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"] = 24;
    //HLT_Mu8_DiEle12_CaloIdL_TrackIdL 
    eleSingleFilters["hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter"] = 25;
    //HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL 
    eleSingleFilters["hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"] = 26;
    //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL 
    eleSingleFilters["hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"] = 27;
    //HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_v
    eleSingleFilters["hltEle27WPTightGsfTrackIsoL1TauJetSeededFilter"] = 28;
    //HLT_Ele30_WPTight_Gsf_L1JetTauSeeded_v
    eleSingleFilters["hltEle30WPTightGsfTrackIsoFilter"] = 29;
    //HLT_Ele32_WPTight_Gsf_L1JetTauSeeded_v
    eleSingleFilters["hltEle32noerWPTightGsfTrackIsoFilter"] = 30;
    //HLT_Ele115_CaloIdVT_GsfTrkIdT_v
    eleSingleFilters["hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter"] = 31;
    //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v
    eleSingleFilters["hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter"] = 32;
    //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v
    eleSingleFilters["hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter"] = 33;
    eleSingleFilters["hltEle8CaloIdLTrackIdLIsoVLTrackIsoFilter"] = 34;
    eleSingleFilters["hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter"] = 35;
    eleSingleFilters["hltMu9Ele9DZFilter"] = 36;
    eleSingleFilters["hltMu8Ele12DZFilter"] = 37;
    //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v
    eleSingleFilters["hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"] = 38;

    //Double electron triggers
    //HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v
    eleDoubleFilters["hltEle24Ele22WPLooseGsfleg1TrackIsoFilter"]                         = 0;
    eleDoubleFilters["hltEle24Ele22WPLooseGsfleg2TrackIsoFilter"]                         = 1;
    //HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v and HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
    eleDoubleFilters["hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter"]               = 2;
    eleDoubleFilters["hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter"]               = 3;
    eleDoubleFilters["hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter"]                         = 4;
    // HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v
    eleDoubleFilters["hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1L1TauJetSeededFilter"] = 5;
    eleDoubleFilters["hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2L1TauJetSeededFilter"] = 6;
    eleDoubleFilters["hltEle23Ele12CaloIdLTrackIdLIsoVLDZL1TauJetSeededFilter"]           = 7;
    eleDoubleFilters["hltEle12Ele12DZFilter"]                                             = 8;

    //triple electron triggers
    eleDoubleFilters["hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter"]                    = 9;
    eleDoubleFilters["hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter"]                    = 10;
    eleDoubleFilters["hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter"]                    = 11;

    muFilters["hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"] = 0; //HLT_IsoMu22_v2
    muFilters["hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"] = 1; //HLT_IsoMu24_v1
    muFilters["hltL3fL1sL1Mu5IsoEG18L1f5L2f7L3Filtered17"] = 2; //HLT_Mu17_Photon*  muon
    muFilters["hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9"] = 3; //HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v3 muon
    muFilters["hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23"] = 4; //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
    muFilters["hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8"] = 5; //LT_Mu8_DiEle12_CaloIdL_TrackIdL_v3 muon
    muFilters["hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"] = 6; //HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL
    muFilters["hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8"] = 7; //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
    muFilters["hltDiMuonGlb17Trk8DzFiltered0p2"] = 8; //HLT_Mu17_TkMu8_DZ
    muFilters["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4"] = 9; //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
    muFilters["hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2"] = 10; //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ
    muFilters["hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4"] = 11; //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
    muFilters["hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2"] = 12; //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ
    muFilters["hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4"] = 13; //HLT_Mu17_TrkIsoVVL
    muFilters["hltDiMuonGlb27Trk8DzFiltered0p2"] = 14; //HLT_Mu27_TkMu8
    muFilters["hltDiMuonGlb30Trk11DzFiltered0p2"] = 15; //HLT_Mu30_TkMu11
    muFilters["hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3Dz0p2L3crIsoRhoFiltered0p15IterTrk02"] = 16; //HLT_DoubleIsoMu17_eta2p1
    muFilters["hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"] = 17; //HLT_IsoMu27
    muFilters["hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09"] = 18; //HLT_IsoTkMu22
    muFilters["hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09"] = 19; //HLT_IsoTkMu24
    muFilters["hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09"] = 20; //HLT_IsoTkMu27
    muFilters["hltL3fL1sL1Mu5IsoEG18ORL1Mu5IsoEG20L1f5L2f7L3Filtered17"] = 21; //HLT_Mu17_Photon*  muon
    muFilters["hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter"] = 23; //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
    muFilters["hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter"] = 24; //HLT_Mu23_TrkIsoVVL_Ele12*
    muFilters["hltMu12TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter"] = 25; //DZ of HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
    muFilters["hltL3fL1DoubleMu7EG7f0Filtered9"] = 26;   // HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v
    muFilters["hltDoubleMu9DZFilter"] = 27;
    muFilters["hltL3fL1Mu6DoubleEG10f0Filtered8"] = 28; // HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_
    muFilters["hltDoubleMu207Mass0to30Photon23L3Filtered"] = 29;  //filter of DoubleMu20_Pho23 (2017)
    muFilters["hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"] = 30; //HLT_IsoMu24 (2017)
    muFilters["hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"] = 31; //HLT_IsoMu27 (2017)
    muFilters["hltL3fL1sMu12Diphoton20L1f0L2f8QL3Filtered12"] = 32; //HLTMu12_DIPho20 (2017)

    phoSingleFilters["hltEG33L1EG26HEFilter"]                = 0;
    phoSingleFilters["hltEG50HEFilter"]                      = 1;
    phoSingleFilters["hltEG75HEFilter"]                      = 2;
    phoSingleFilters["hltEG90HEFilter"]                      = 3;
    phoSingleFilters["hltEG120HEFilter"]                     = 4;
    phoSingleFilters["hltEG150HEFilter"]                     = 5;
    phoSingleFilters["hltEG175HEFilter"]                     = 6;
    phoSingleFilters["hltEG200HEFilter"]                     = 7; 
    phoSingleFilters["hltEG300erEtFilter"]                   = 8;
    phoSingleFilters["hltHtEcal800"]                         = 9;
    phoSingleFilters["hltEG50R9Id90HE10IsoMTrackIsoFilter"]  = 10;
    phoSingleFilters["hltEG75R9Id90HE10IsoMTrackIsoFilter"]  = 11;
    phoSingleFilters["hltEG90R9Id90HE10IsoMTrackIsoFilter"]  = 12;
    phoSingleFilters["hltEG120R9Id90HE10IsoMTrackIsoFilter"] = 13;
    phoSingleFilters["hltEG165R9Id90HE10IsoMTrackIsoFilter"] = 14;

    //L1 seed for diphoton triggers  
    phoDoubleFilters["hltSingleEGL1SingleEG40ORL1SingleEG25ORL1DoubleEG2210ORL1DoubleEG1510Filter"]         = 0;
    //For path HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v    
    phoDoubleFilters["hltEG18Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededLastFilter"]                    = 1;
    phoDoubleFilters["hltEG18R9Id85b90eHE12R9Id50b80eR9UnseededLastFilter"]                                 = 2;
    //For HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v
    phoDoubleFilters["hltEG18Iso60CaloId15b35eHE12R9Id50b80eTrackIsoUnseededpixSeedLastFilter"]             = 3;
    phoDoubleFilters["hltEG18R9Id85b90eHE12R9Id50b80eR9pixSeedUnseededLastFilter"]                          = 4;
    //For both HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v and 
    //HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v
    phoDoubleFilters["hltEG30LR9Id85b90eHE12R9Id50b80eR9IdLastFilter"]                                      = 5;
    phoDoubleFilters["hltEG30LIso60CaloId15b35eHE12R9Id50b80eEcalIsoLastFilter"]                            = 6;
    //For path HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55_v           
    phoDoubleFilters["hltEG18Iso60CaloId15b35eHE10R9Id50b80eTrackIsoSolidUnseededLastFilter"]               = 7;
    phoDoubleFilters["hltEG18R9Id85b90eHE10R9Id50b80eR9UnseededLastFilter"]                                 = 8;
    phoDoubleFilters["hltEG30R9Id85b90eHE10R9Id50b80eR9IdLastFilter"]                                       = 9;
    phoDoubleFilters["hltEG30Iso60CaloId15b35eHE10R9Id50b80eEcalIsoFilter"]                                 = 10;
    //For path HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v 
    phoDoubleFilters["hltEG18PVIso60CaloId15b35eHE10R9Id50b80eTrackIsoUnseededDoublePixelVetoLastFilter"]   = 11;
    phoDoubleFilters["hltEG18PVR9Idb85e90HE10R9Id50b80eR9DoublePixelVetoUnseededLastFilter"]                = 12;
    phoDoubleFilters["hltEG30PVR9Idb85e90HE10R9Id50b80eR9IdLastFilter"]                                     = 13;
    phoDoubleFilters["hltEG30PVIso60CaloId15b35eHE10R9Id50b80eEcalIsoFilter"]                               = 14;
    //For path HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v
    phoDoubleFilters["hltEG18EBR9Idb85e90HE10R9Id50b80eR9DoublePixelVetoUnseededLastFilter"]                = 15;
    phoDoubleFilters["hltEG18EBIso60CaloId15b35eHE10R9Id50b80eTrackIsoUnseededDoublePixelVetoLastFilter"]   = 16;
    phoDoubleFilters["hltEG30EBR9Idb85e90HE10R9Id50b80eR9IdLastFilter"]                                     = 17;
    phoDoubleFilters["hltEG30EBIso60CaloId15b35eHE10R9Id50b80eEcalIsoFilter"]                               = 18;
    //For path HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60_v      
    phoDoubleFilters["hltEG16R9Id85HE10R9Id65R9IdUnseededLastFilter"]                                       = 19;
    phoDoubleFilters["hltEG16Iso50T80LCaloId24b40eHE10R9Id65TrackIsoUnseededLastFilter"]                    = 20;
    phoDoubleFilters["hltEG26Iso50T80LCaloId24b40eHE10R9Id65Eta2HcalIsoLastFilter"]                         = 21;
    phoDoubleFilters["hltEG26R9Id85HE10R9Id65R9IdEta2LastFilter"]                                           = 22;
    //For path HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_v    
    phoDoubleFilters["hltEG22R9Id85HE10R9Id65R9IdUnseededLastFilter"]                                       = 23;
    phoDoubleFilters["hltEG22Iso50T80LCaloId24b40eHE10R9Id65TrackIsoUnseededLastFilter"]                    = 24;
    phoDoubleFilters["hltEG36R9Id85HE10R9Id65R9IdEta2LastFilter"]                                           = 25;
    phoDoubleFilters["hltEG36Iso50T80LCaloId24b40eHE10R9Id65Eta2HcalIsoLastFilter"]                         = 26;
    phoDoubleFilters["hltMu17Photon22CaloIdLL1ISOHEFilter"]                                                 = 27;
    phoDoubleFilters["hltMu17Photon30CaloIdLL1ISOHEFilter"]                                                 = 28;
    phoDoubleFilters["hltMu17Photon30CaloIdLL1ISOORHEFilter"]                                               = 29;
    phoDoubleFilters["hltEG70HEFilter"]                                                                     = 30; //seeded leg of DoublePho70
    phoDoubleFilters["hltDiEG70HEUnseededFilter"]                                                           = 31; //unseeded leg of DoublePho70
    phoDoubleFilters["hltMu12DiEG20HEUnseededFilter"]                                                       = 32; //filter of Mu12_DiPho20

    //triple photon filter                                                           
    phoTripleFilters["hltEG20CaloIdLV2ClusterShapeL1TripleEGFilter"]         = 0;  //1st leg of TriplePho20_Calo
    phoTripleFilters["hltTriEG20CaloIdLV2ClusterShapeUnseededFilter"]        = 1;  //3rd leg of TriplePho20_Calo

    phoTripleFilters["hltEG20CaloIdLV2R9IdVLR9IdL1TripleEGFilter"]           = 2;  //1st leg of TriplePho20_Calo_R9
    phoTripleFilters["hltTriEG20CaloIdLV2R9IdVLR9IdUnseededFilter"]          = 3;  //3rd leg of TriplePho20_Calo_R9

    phoTripleFilters["hltEG30CaloIdLV2ClusterShapeL1TripleEGFilter"]         = 4;  //1st leg Seeded_30 of TriplePho30_30_10_Calo
    phoTripleFilters["hltEG10CaloIdLV2ClusterShapeUnseededFilter"]           = 5;  //1st leg UnSeeded_10 of TriplePho30_30_10_Calo
    phoTripleFilters["hltDiEG30CaloIdLV2EtUnseededFilter"]                   = 6;  //2nd leg UnSeeded_30 of TriplePho30_30_10_Calo
    phoTripleFilters["hltDiEG30TriEG10CaloIdLV2EtUnseededFilter"]            = 7;  //3 legs of TriplePho30_30_10_Calo (last filter)

    phoTripleFilters["hltEG30CaloIdLV2R9IdVLR9IdL1TripleEGFilter"]           = 8;  //1st leg Seeded_30 of TriplePho30_30_10_Calo_R9
    phoTripleFilters["hltEG10CaloIdLV2R9IdVLR9IdUnseededFilter"]             = 9;  //1st leg UnSeeded_10 of TriplePho30_30_10_Calo_R9
    phoTripleFilters["hltDiEG30CaloIdLV2R9IdVLEtUnseededFilter"]             = 10; //2nd leg UnSeeded_30 of TriplePho30_30_10_Calo_R9
    phoTripleFilters["hltDiEG30TriEG10HECaloIdV2R9IdVLEtUnseededFilter"]     = 11; //3 legs of TriplePho30_30_10_Calo_R9 (last filter)

    phoTripleFilters["hltEG35CaloIdLV2R9IdVLR9IdL1TripleEGFilter"]           = 12; //1st leg Seeded_35 of TriplePho35_35_5_Calo_R9
    phoTripleFilters["hltEG5CaloIdLV2R9IdVLR9IdUnseededFilter"]              = 13; //1st leg UnSeeded_5 of TriplePho35_35_5_Calo_R9
    phoTripleFilters["hltDiEG35CaloIdLV2R9IdVLEtUnseededFilter"]             = 14; //2nd leg UnSeeded_35 of TriplePho35_35_5_Calo_R9
    phoTripleFilters["hltDiEG35TriEG5CaloIdLV2R9IdVLEtUnseededFilter"]       = 15; //3 legs of TriplePho35_35_5_Calo_R9 (last filter)
    
    jetFilters["hltSinglePFJet40"]  =  0;
    jetFilters["hltSinglePFJet60"]  =  1;
    jetFilters["hltSinglePFJet80"]  =  2;
    jetFilters["hltSinglePFJet140"] =  3;
    jetFilters["hltSinglePFJet200"] =  4;
    jetFilters["hltSinglePFJet260"] =  5;
    jetFilters["hltSinglePFJet320"] =  6;
    jetFilters["hltSinglePFJet400"] =  7;
    jetFilters["hltSinglePFJet450"] =  8;
    jetFilters["hltSinglePFJet500"] =  9;

    l1Filters["hltL1sSingleEG26"] = 0;
    l1Filters["hltL1sSingleEG40"] = 1;
    l1Filters["hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24er"] = 2;
    l1Filters["hltL1sSingleIsoEG20erIorSingleIsoEG22erIorSingleEG40"] = 3;
    l1Filters["hltL1sSingleIsoEG22erIorSingleIsoEG24erIorSingleEG40"] = 4;
    l1Filters["hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26"] = 5;
    l1Filters["hltL1sHTT200IorHTT220IorHTT240lorHTT255IorHTT270lorHTT280lorHTT300lorHTT320"] = 6;
    l1Filters["hltL1sEG27erHTT200IorHTT280IorHTT300"] = 7;
    l1Filters["hltL1sIsoEG22erTau20erdEtaMin0p2"] = 8;
    l1Filters["hltL1sSingleMu18"] = 9;
    l1Filters["hltL1sSingleMu22"] = 10;
    l1Filters["hltL1sSingleMu22Or25"] = 11;
    l1Filters["hltL1sMu5IsoEG18"] = 12; //HLT_Mu17_Photon*
    l1Filters["hltL1sMu20EG10"] = 13; //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
    l1Filters["hltL1sDoubleMu7EG7"] = 14; //HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v3
    l1Filters["hltL1sMu6DoubleEG10"] = 15; //LT_Mu8_DiEle12_CaloIdL_TrackIdL_v3
    l1Filters["hltL1sMu5EG15"] = 16; //HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3
    l1Filters["hltL1sMu5EG20IorMu5IsoEG18"] = 17; //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3
    l1Filters["hltL1sDoubleMu114IorDoubleMu125"] = 18; //DiMuon
    l1Filters["hltL1sSingleMu10LowQ"] = 19; //HLT_Mu17_TrkIsoVVL_v2
    l1Filters["hltL1sL1SingleEG40ORL1SingleEG35ORL1DoubleEG2210ORL1DoubleEG1510"] = 20;
    l1Filters["hltL1sDiPhotonSeeds"] = 21;
    l1Filters["hltL1sSingleAndDoubleEGor"] = 22; //HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v and HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
    l1Filters["hltL1sSingleEG40IorSingleIsoEG26IorSingleIsoEG24IorDoubleEG1510IorDoubleEG2212IorDoubleEG1817"] = 23;
    l1Filters["hltL1sSingleEGor"] = 24; //HLT_Ele27_WPTight_Gsf_v
    l1Filters["hltL1sSingleJetAndTauHighPtOr"] = 25; //HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_v
    l1Filters["hltL1sSingleEGNonIsoOrWithJetAndTau"] = 26; //HLT_Ele115_CaloIdVT_GsfTrkIdT_v
    l1Filters["hltL1sDoubleEG2210IorDoubleEG2512"] = 27; //HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v
    l1Filters["hltL1sSingleEG34IorSingleEG40"] = 28; //HLT_Photon90_CaloIdL_PFHT600_v
    l1Filters["hltL1sSingleEG18IorSingleEG24"] = 29;
    l1Filters["hltL1sSingleEG10IorSingleEG15"] = 30;
    l1Filters["hltL1sSingleAndDoubleEGNonIsoOrWithEG26WithJetAndTau"] = 31;
    
  }
 
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerHandleMiniAOD;
  e.getByToken(triggerObjectsLabel_, triggerHandleMiniAOD);

  edm::Handle<edm::TriggerResults> trgResultsHandle;
  e.getByToken(trgResultsLabel_, trgResultsHandle);

  edm::Handle<std::string> filterLabels_;
  e.getByLabel("slimmedPatTrigger:filterLabels", filterLabels_);

  //const edm::TriggerNames &names = e.triggerNames(*trgResultsHandle);

  for (pat::TriggerObjectStandAlone obj : *triggerHandleMiniAOD) {
    //obj.unpackPathNames(names);
    //obj.unpackPathNames(e);
    obj.unpackFilterLabels(e, *trgResultsHandle);

    // loop over filters    
    for (size_t iF = 0; iF < obj.filterLabels().size(); ++iF) {
      string label = obj.filterLabels()[iF];

      //cout<<"label : "<<iF<<" "<<label<<endl;

      std::map<string,size_t>::iterator idxEleSingle = eleSingleFilters.find(label);
      std::map<string,size_t>::iterator idxEleDouble = eleDoubleFilters.find(label);
      std::map<string,size_t>::iterator idxPhoSingle = phoSingleFilters.find(label);
      std::map<string,size_t>::iterator idxPhoDouble = phoDoubleFilters.find(label);
      std::map<string,size_t>::iterator idxPhoTriple = phoTripleFilters.find(label);
      std::map<string,size_t>::iterator idxMu  = muFilters.find(label);
      std::map<string,size_t>::iterator idxJet = jetFilters.find(label);
      std::map<string,size_t>::iterator idxL1 = l1Filters.find(label);

      // single electron filters
      if (idxEleSingle != eleSingleFilters.end()) {
        size_t idx = idxEleSingle->second;
        trgSingleElePt [idx].push_back(obj.pt());
        trgSingleEleEta[idx].push_back(obj.eta());
        trgSingleElePhi[idx].push_back(obj.phi());
	//cout<<idx<<" "<<obj.pt()<<" "<<obj.eta()<<" "<<obj.phi()<<endl;
      }

      // double electron filters
      if (idxEleDouble != eleDoubleFilters.end()) {
        size_t idx = idxEleDouble->second;
        trgDoubleElePt [idx].push_back(obj.pt());
        trgDoubleEleEta[idx].push_back(obj.eta());
        trgDoubleElePhi[idx].push_back(obj.phi());
      }

      // single photon filters
      if (idxPhoSingle != phoSingleFilters.end()) {
        size_t idx = idxPhoSingle->second;
        trgSinglePhoPt [idx].push_back(obj.pt());
        trgSinglePhoEta[idx].push_back(obj.eta());
        trgSinglePhoPhi[idx].push_back(obj.phi());
      }

      // double photon filters
      if (idxPhoDouble != phoDoubleFilters.end()) {
        size_t idx = idxPhoDouble->second;
        trgDoublePhoPt [idx].push_back(obj.pt());
        trgDoublePhoEta[idx].push_back(obj.eta());
        trgDoublePhoPhi[idx].push_back(obj.phi());
      }

      //triple photon filters
      if (idxPhoTriple != phoTripleFilters.end()) {
        size_t idx = idxPhoTriple->second;
        trgTriplePhoPt [idx].push_back(obj.pt());
        trgTriplePhoEta[idx].push_back(obj.eta());
	trgTriplePhoPhi[idx].push_back(obj.phi());
      }

      // muon filters
      if (idxMu != muFilters.end()) {
        size_t idx = idxMu->second;
        trgMuPt [idx].push_back(obj.pt());
        trgMuEta[idx].push_back(obj.eta());
        trgMuPhi[idx].push_back(obj.phi());
      }

      // jet filters
      if (idxJet != jetFilters.end()) {
        size_t idx = idxJet->second;
        trgJetPt [idx].push_back(obj.pt());
        trgJetEta[idx].push_back(obj.eta());
        trgJetPhi[idx].push_back(obj.phi());
      }

      // L1 filters
      if (idxL1 != l1Filters.end()) {
        size_t idx = idxL1->second;
        trgL1Eta[idx].push_back(obj.eta());
        trgL1Phi[idx].push_back(obj.phi());
      }
    }
  }

}

ULong64_t ggNtuplizer::matchSingleElectronTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgSingleElePt[f].size(); ++v)
      if (fabs(pt - trgSingleElePt[f][v])/trgSingleElePt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgSingleEleEta[f][v], trgSingleElePhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

ULong64_t ggNtuplizer::matchDoubleElectronTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgDoubleElePt[f].size(); ++v)
      if (fabs(pt - trgDoubleElePt[f][v])/trgDoubleElePt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgDoubleEleEta[f][v], trgDoubleElePhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

ULong64_t ggNtuplizer::matchSinglePhotonTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgSinglePhoPt[f].size(); ++v)
      if (fabs(pt - trgSinglePhoPt[f][v])/trgSinglePhoPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgSinglePhoEta[f][v], trgSinglePhoPhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

ULong64_t ggNtuplizer::matchDoublePhotonTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgDoublePhoPt[f].size(); ++v)
      if (fabs(pt - trgDoublePhoPt[f][v])/trgDoublePhoPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgDoublePhoEta[f][v], trgDoublePhoPhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

ULong64_t ggNtuplizer::matchTriplePhotonTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above 
  UInt_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgTriplePhoPt[f].size(); ++v)
      if (fabs(pt - trgTriplePhoPt[f][v])/trgTriplePhoPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgTriplePhoEta[f][v], trgTriplePhoPhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}



ULong64_t ggNtuplizer::matchMuonTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgMuPt[f].size(); ++v)
      if (fabs(pt - trgMuPt[f][v])/trgMuPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgMuEta[f][v], trgMuPhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

ULong64_t ggNtuplizer::matchJetTriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgJetPt[f].size(); ++v)
      if (fabs(pt - trgJetPt[f][v])/trgJetPt[f][v] < trgFilterDeltaPtCut_ &&
          deltaR(eta, phi, trgJetEta[f][v], trgJetPhi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

ULong64_t ggNtuplizer::matchL1TriggerFilters(double pt, double eta, double phi) {

  // bits in the return value correspond to decisions from filters defined above
  ULong64_t result = 0;

  for (size_t f = 0; f < 64; ++f)
    for (size_t v = 0; v < trgL1Eta[f].size(); ++v)
      if (deltaR(eta, phi, trgL1Eta[f][v], trgL1Phi[f][v]) < trgFilterDeltaRCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

Double_t ggNtuplizer::deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}

Double_t ggNtuplizer::deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}
