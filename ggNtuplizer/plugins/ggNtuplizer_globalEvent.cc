#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t       run_;
Long64_t    event_;
Int_t       lumis_;
Bool_t      isData_;
Int_t       nVtx_;
Int_t       nGoodVtx_;
Bool_t      isPVGood_;
float       vtx_;
float       vty_;
float       vtz_;
float       rho_;
float       rhoCentral_;
ULong64_t   HLTEleMuX_;
ULong64_t   HLTPho_;
ULong64_t   HLTPhoRejectedByPS_;
ULong64_t   HLTJet_;
ULong64_t   HLTEleMuXIsPrescaled_;
ULong64_t   HLTPhoIsPrescaled_;
ULong64_t   HLTJetIsPrescaled_;
vector<int> phoPrescale_;

void ggNtuplizer::branchesGlobalEvent(TTree* tree) {

  tree->Branch("run",     &run_);
  tree->Branch("event",   &event_);
  tree->Branch("lumis",   &lumis_);
  tree->Branch("isData",  &isData_);
  tree->Branch("nVtx",                 &nVtx_);
  tree->Branch("nGoodVtx",             &nGoodVtx_);
  tree->Branch("isPVGood",             &isPVGood_);
  tree->Branch("vtx",                  &vtx_); 
  tree->Branch("vty",                  &vty_); 
  tree->Branch("vtz",                  &vtz_); 
  tree->Branch("rho",                  &rho_);
  tree->Branch("rhoCentral",           &rhoCentral_);
  tree->Branch("HLTEleMuX",            &HLTEleMuX_);
  tree->Branch("HLTPho",               &HLTPho_);
  tree->Branch("HLTPhoRejectedByPS",   &HLTPhoRejectedByPS_);
  tree->Branch("HLTJet",               &HLTJet_);
  tree->Branch("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled_);
  tree->Branch("HLTPhoIsPrescaled",    &HLTPhoIsPrescaled_);
  tree->Branch("HLTJetIsPrescaled",    &HLTJetIsPrescaled_);

  if (!doGenParticles_)
    tree->Branch("phoPrescale",          &phoPrescale_); 
}

void ggNtuplizer::fillGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {

  if (!doGenParticles_) {
    phoPrescale_.clear();
    phoPrescale_.reserve(9);
  }

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);

  edm::Handle<double> rhoCentralHandle;
  e.getByToken(rhoCentralLabel_, rhoCentralHandle);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();
  rho_    = *(rhoHandle.product());
  if (rhoCentralHandle.isValid()) rhoCentral_ = *(rhoCentralHandle.product());
  else rhoCentral_ = -99.;

  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);

  nVtx_     = -1;
  nGoodVtx_ = -1;
  if (vtxHandle.isValid()) {
    nVtx_     = 0;
    nGoodVtx_ = 0;    

    for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {

      if (nVtx_ == 0) {
	vtx_     = v->x();
	vty_     = v->y();
	vtz_     = v->z();

	isPVGood_ = false;
	if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) isPVGood_ = true;
      }

      if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) nGoodVtx_++;
      nVtx_++;

    }
  } else
    edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";

  // HLT treatment
  HLTEleMuX_            = 0;
  HLTPho_               = 0;
  HLTPhoRejectedByPS_   = 0;
  HLTJet_               = 0;
  HLTEleMuXIsPrescaled_ = 0;
  HLTPhoIsPrescaled_    = 0;
  HLTJetIsPrescaled_    = 0;

  edm::Handle<edm::TriggerResults> trgResultsHandle;
  e.getByToken(trgResultsLabel_, trgResultsHandle);

  bool cfg_changed = true;
  hltPrescaleProvider_.init(e.getRun(), es, trgResultsProcess_, cfg_changed);
  HLTConfigProvider const&  hltCfg = hltPrescaleProvider_.hltConfigProvider();
  const int prescaleSet = hltPrescaleProvider_.prescaleSet(e,es);

  const edm::TriggerNames &trgNames = e.triggerNames(*trgResultsHandle);

  for (size_t i = 0; i < trgNames.size(); ++i) {
    const string &name = trgNames.triggerName(i);

    // HLT name => bit correspondence
    int bitEleMuX = -1;
    int bitPho    = -1;
    int bitJet    = -1;

    if (year_ == 2016) {
      if      (name.find("HLT_Ele25_eta2p1_WPTight_Gsf_v")                      != string::npos) bitEleMuX =  0;
      else if (name.find("HLT_Ele27_eta2p1_WPTight_Gsf_v")                      != string::npos) bitEleMuX =  1; 
      else if (name.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v")                      != string::npos) bitEleMuX =  2;
      else if (name.find("HLT_Ele32_eta2p1_WPTight_Gsf_v")                      != string::npos) bitEleMuX =  3; 
      else if (name.find("HLT_Ele27_WPTight_Gsf_v")                             != string::npos) bitEleMuX =  4; 
      else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")         != string::npos) bitEleMuX =  5; 
      else if (name.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")             != string::npos) bitEleMuX =  6; 
      else if (name.find("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")           != string::npos) bitEleMuX =  7; 
      else if (name.find("HLT_Mu17_Photon30_CaloIdL_L1ISO_v")                   != string::npos) bitEleMuX =  8; 
      else if (name.find("HLT_Mu17_Photon35_CaloIdL_L1ISO_v")                   != string::npos) bitEleMuX =  9; 
      else if (name.find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v")             != string::npos) bitEleMuX = 10;
      else if (name.find("HLT_DoubleEle33_CaloIdL_MW_v")                        != string::npos) bitEleMuX = 11;
      else if (name.find("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v")          != string::npos) bitEleMuX = 12;
      else if (name.find("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v")             != string::npos) bitEleMuX = 13;
      else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")               != string::npos) bitEleMuX = 14;
      else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")             != string::npos) bitEleMuX = 15;
      else if (name.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")           != string::npos) bitEleMuX = 16;
      else if (name.find("HLT_Mu30_TkMu11_v")                                   != string::npos) bitEleMuX = 17;
      else if (name.find("HLT_DoubleIsoMu17_eta2p1_noDzCut_v")                  != string::npos) bitEleMuX = 18;
      else if (name.find("HLT_IsoMu24_v")                                       != string::npos) bitEleMuX = 19;
      else if (name.find("HLT_IsoTkMu24_v")                                     != string::npos) bitEleMuX = 20;
      else if (name.find("HLT_Mu50_v")                                          != string::npos) bitEleMuX = 21;
      else if (name.find("HLT_TripleMu_12_10_5_v")                              != string::npos) bitEleMuX = 22;
      else if (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")    != string::npos) bitEleMuX = 23;
      else if (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != string::npos) bitEleMuX = 24;
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v")    != string::npos) bitEleMuX = 25;
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v") != string::npos) bitEleMuX = 26;
      else if (name.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")                  != string::npos) bitEleMuX = 27;
      else if (name.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")                   != string::npos) bitEleMuX = 28;
      else if (name.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v")      != string::npos) bitEleMuX = 29;
      else if (name.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v")      != string::npos) bitEleMuX = 30;
      else if (name.find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v")                != string::npos) bitEleMuX = 31;
      else if (name.find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v")                != string::npos) bitEleMuX = 32;
      else if (name.find("HLT_Ele17_Ele12_CaloId_TrackId_Iso_DZ_v")             != string::npos) bitEleMuX = 33;
      else if (name.find("HLT_DoubleEle33_CaloId_GsfTrackIdVL_v")               != string::npos) bitEleMuX = 34;
      else if (name.find("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_v")              != string::npos) bitEleMuX = 35; 
      else if (name.find("HLT_Ele30_WPTight_Gsf_L1JetTauSeeded_v")              != string::npos) bitEleMuX = 36; 
      else if (name.find("HLT_Ele32_WPTight_Gsf_L1JetTauSeeded_v")              != string::npos) bitEleMuX = 37; 
      else if (name.find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v")                     != string::npos) bitEleMuX = 38;
      else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v") != string::npos) bitEleMuX = 39;
      else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")            != string::npos) bitEleMuX = 40;
      else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")                  != string::npos) bitEleMuX = 41;
      else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v")                != string::npos) bitEleMuX = 42;
      else if (name.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")                  != string::npos) bitEleMuX = 43;
      else if (name.find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")                  != string::npos) bitEleMuX = 44;
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")   != string::npos) bitEleMuX = 45;
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")!= string::npos) bitEleMuX = 46;
      else if (name.find("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")                  != string::npos) bitEleMuX = 47;
      else if (name.find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")          != string::npos) bitEleMuX = 48;
      else if (name.find("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v")          != string::npos) bitEleMuX = 49;
      else if (name.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")          != string::npos) bitEleMuX = 50;

      // Photon triggers
      if      (name.find("HLT_Photon22_v")                    != string::npos) bitPho =  0; //bit0(lowest)
      else if (name.find("HLT_Photon30_v")                    != string::npos) bitPho =  1; 
      else if (name.find("HLT_Photon36_v")                    != string::npos) bitPho =  2; 
      else if (name.find("HLT_Photon50_v")                    != string::npos) bitPho =  3; 
      else if (name.find("HLT_Photon75_v")                    != string::npos) bitPho =  4; 
      else if (name.find("HLT_Photon90_v")                    != string::npos) bitPho =  5; 
      else if (name.find("HLT_Photon120_v")                   != string::npos) bitPho =  6; 
      else if (name.find("HLT_Photon175_v")                   != string::npos) bitPho =  7; 
      else if (name.find("HLT_Photon250_NoHE_v")              != string::npos) bitPho =  8; 
      else if (name.find("HLT_Photon300_NoHE_v")              != string::npos) bitPho =  9; 
      else if (name.find("HLT_Photon500_v")                   != string::npos) bitPho = 10; 
      else if (name.find("HLT_Photon600_v")                   != string::npos) bitPho = 11; 
      else if (name.find("HLT_Photon165_HE10_v")              != string::npos) bitPho = 12; 
      else if (name.find("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v") != string::npos) bitPho = 13; // exist
      else if (name.find("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v")                             != string::npos) bitPho = 14; // used
      else if (name.find("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v")        != string::npos) bitPho = 15; // exist
      else if (name.find("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v")        != string::npos) bitPho = 16; // used
      else if (name.find("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v")         != string::npos) bitPho = 17; // used
      else if (name.find("HLT_Photon135_PFMET100_v")                          != string::npos) bitPho = 18; 
      else if (name.find("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v")  != string::npos) bitPho = 19;
      else if (name.find("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v")       != string::npos) bitPho = 20;
      else if (name.find("HLT_Photon90_CaloIdL_PFHT600_v")                    != string::npos) bitPho = 21;
      else if (name.find("HLT_DoublePhoton60_v")                              != string::npos) bitPho = 22;
      else if (name.find("HLT_DoublePhoton85_v")                              != string::npos) bitPho = 23;
      else if (name.find("HLT_Photon22_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 24;
      else if (name.find("HLT_Photon50_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 25;
      else if (name.find("HLT_Photon75_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 26;
      else if (name.find("HLT_Photon90_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 27;
      else if (name.find("HLT_Photon120_R9Id90_HE10_IsoM_v")                  != string::npos) bitPho = 28;
      else if (name.find("HLT_Photon165_R9Id90_HE10_IsoM_v")                  != string::npos) bitPho = 29;
      else if (name.find("HLT_ECALHT800_v")                                   != string::npos) bitPho = 30;

      // Jet triggers
      if      (name.find("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_v")                    != string::npos) bitJet =  0;
      else if (name.find("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_v")                    != string::npos) bitJet =  1; 
      else if (name.find("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_v")                != string::npos) bitJet =  2; 
      else if (name.find("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v")                != string::npos) bitJet =  3; 
      else if (name.find("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v")                      != string::npos) bitJet =  4;
      else if (name.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v")                    != string::npos) bitJet =  5;
      else if (name.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v")                    != string::npos) bitJet =  6;
      else if (name.find("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != string::npos) bitJet =  7;
      else if (name.find("HLT_PFMET170_HBHECleaned_v")                                 != string::npos) bitJet =  8;
      else if (name.find("HLT_CaloJet500_NoJetID_v")                                   != string::npos) bitJet =  9;
      else if (name.find("HLT_PFJet40_v")                                              != string::npos) bitJet = 10; 
      else if (name.find("HLT_PFJet60_v")                                              != string::npos) bitJet = 11; 
      else if (name.find("HLT_PFJet80_v")                                              != string::npos) bitJet = 12; 
      else if (name.find("HLT_PFJet140_v")                                             != string::npos) bitJet = 13; 
      else if (name.find("HLT_PFJet200_v")                                             != string::npos) bitJet = 14; 
      else if (name.find("HLT_PFJet260_v")                                             != string::npos) bitJet = 15; 
      else if (name.find("HLT_PFJet320_v")                                             != string::npos) bitJet = 16; 
      else if (name.find("HLT_PFJet400_v")                                             != string::npos) bitJet = 17; 
      else if (name.find("HLT_PFJet450_v")                                             != string::npos) bitJet = 18;     
      else if (name.find("HLT_PFJet500_v")                                             != string::npos) bitJet = 19; 
      else if (name.find("HLT_AK8PFHT700_TrimR0p1PT0p3Mass50_v")                       != string::npos) bitJet = 20; 
      else if (name.find("HLT_AK8PFJet360_TrimMass30_v")                               != string::npos) bitJet = 21;
      else if (name.find("HLT_PFHT300_PFMET110_v")                                     != string::npos) bitJet = 22;
      else if (name.find("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v")   != string::npos) bitJet = 23;
      else if (name.find("HLT_PFMET170_HBHE_BeamHaloCleaned_v")                        != string::npos) bitJet = 24;
      else if (name.find("HLT_PFMET300_v")                                             != string::npos) bitJet = 25;
      else if (name.find("HLT_PFMET110_PFMHT110_IDTight_v")                            != string::npos) bitJet = 26;
      else if (name.find("HLT_PFMET120_PFMHT120_IDTight_v")                            != string::npos) bitJet = 27;
      else if (name.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v")                    != string::npos) bitJet = 28;
      else if (name.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v")                    != string::npos) bitJet = 29;
      else if (name.find("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != string::npos) bitJet = 30;
      else if (name.find("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != string::npos) bitJet = 31;
      else if (name.find("HLT_PFHT800_v")                                              != string::npos) bitJet = 32;
      else if (name.find("HLT_PFHT900_v")                                              != string::npos) bitJet = 33;
      else if (name.find("HLT_PFHT750_4JetPt50_v")                                     != string::npos) bitJet = 34;
      else if (name.find("HLT_PFHT750_4JetPt70_v")                                     != string::npos) bitJet = 35;
      else if (name.find("HLT_PFHT800_4JetPt50_v")                                     != string::npos) bitJet = 36;
    }

    if (year_ == 2017) {
      // Electron or Muon or Cross triggers
      if      (name.find("HLT_DoubleMu20_7_Mass0to30_L1_DM4_v")                  != string::npos) bitEleMuX =  0; // 2017
      else if (name.find("HLT_DoubleMu20_7_Mass0to30_Photon23_v")                != string::npos) bitEleMuX =  1; // 2017
      else if (name.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v")                       != string::npos) bitEleMuX =  2;
      else if (name.find("HLT_Ele35_WPTight_Gsf_v")                              != string::npos) bitEleMuX =  3; // 2017
      else if (name.find("HLT_Ele27_WPTight_Gsf_v")                              != string::npos) bitEleMuX =  4; // 2017
      else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")          != string::npos) bitEleMuX =  5; // 2017
      else if (name.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")              != string::npos) bitEleMuX =  6; // 2017
      else if (name.find("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v")            != string::npos) bitEleMuX =  7; 
      else if (name.find("HLT_Mu17_Photon30_CaloIdL_L1ISO_v")                    != string::npos) bitEleMuX =  8; 
      else if (name.find("HLT_Mu17_Photon35_CaloIdL_L1ISO_v")                    != string::npos) bitEleMuX =  9; 
      else if (name.find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v")              != string::npos) bitEleMuX = 10;
      else if (name.find("HLT_DoubleEle33_CaloIdL_MW_v")                         != string::npos) bitEleMuX = 11; // 2017
      else if (name.find("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v")           != string::npos) bitEleMuX = 12;
      else if (name.find("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v")              != string::npos) bitEleMuX = 13;
      else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v")          != string::npos) bitEleMuX = 14; // 2017
      else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v")        != string::npos) bitEleMuX = 15; // 2017
      else if (name.find("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")            != string::npos) bitEleMuX = 16;
      else if (name.find("HLT_TripleMu_10_5_5_DZ_v")                             != string::npos) bitEleMuX = 17; // 2017
      else if (name.find("HLT_DoubleIsoMu20_eta2p1_v")                           != string::npos) bitEleMuX = 18; // 2017
      else if (name.find("HLT_IsoMu27_v")                                        != string::npos) bitEleMuX = 19; // 2017
      else if (name.find("HLT_IsoTkMu24_v")                                      != string::npos) bitEleMuX = 20; // 2017
      else if (name.find("HLT_Mu50_v")                                           != string::npos) bitEleMuX = 21; // 2017
      else if (name.find("HLT_TripleMu_12_10_5_v")                               != string::npos) bitEleMuX = 22; // 2017
      else if (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")    != string::npos) bitEleMuX = 23; // 2017*
      else if (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != string::npos) bitEleMuX = 24; // 2017
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v")     != string::npos) bitEleMuX = 25; 
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v")  != string::npos) bitEleMuX = 26;
      else if (name.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")                   != string::npos) bitEleMuX = 27; // 2017
      else if (name.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")                    != string::npos) bitEleMuX = 28; // 2017
      else if (name.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v")       != string::npos) bitEleMuX = 29;
      else if (name.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v")       != string::npos) bitEleMuX = 30;
      else if (name.find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v")                 != string::npos) bitEleMuX = 31;
      else if (name.find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v")                 != string::npos) bitEleMuX = 32;
      else if (name.find("HLT_Ele17_Ele12_CaloId_TrackId_Iso_DZ_v")              != string::npos) bitEleMuX = 33;
      else if (name.find("HLT_DoubleEle33_CaloId_GsfTrackIdVL_v")                != string::npos) bitEleMuX = 34;
      else if (name.find("HLT_Ele27_WPTight_Gsf_L1JetTauSeeded_v")               != string::npos) bitEleMuX = 35; 
      else if (name.find("HLT_Ele30_WPTight_Gsf_L1JetTauSeeded_v")               != string::npos) bitEleMuX = 36; 
      else if (name.find("HLT_Ele32_WPTight_Gsf_L1JetTauSeeded_v")               != string::npos) bitEleMuX = 37; 
      else if (name.find("HLT_Ele115_CaloIdVT_GsfTrkIdT_v")                      != string::npos) bitEleMuX = 38; // 2017
      else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_L1JetTauSeeded_v") != string::npos) bitEleMuX = 39;
      else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")             != string::npos) bitEleMuX = 40; // 2017
      else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")                   != string::npos) bitEleMuX = 41; // 2017*
      else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")                != string::npos) bitEleMuX = 42; // 2017*
      else if (name.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")                   != string::npos) bitEleMuX = 43;
      else if (name.find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")                   != string::npos) bitEleMuX = 44;
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")    != string::npos) bitEleMuX = 45; // 2017
      else if (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != string::npos) bitEleMuX = 46; // 2017
      else if (name.find("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")                   != string::npos) bitEleMuX = 47;
      else if (name.find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")           != string::npos) bitEleMuX = 48;
      else if (name.find("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v")                 != string::npos) bitEleMuX = 49; // 2017
      else if (name.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")           != string::npos) bitEleMuX = 50; // 2017
      else if (name.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v")                 != string::npos) bitEleMuX = 51; // 2017
      else if (name.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v")                != string::npos) bitEleMuX = 52; // 2017
      else if (name.find("HLT_Mu12_DoublePhoton20_v")                            != string::npos) bitEleMuX = 53; // 2017
      
      // Photon triggers
      if      (name.find("HLT_Photon22_v")                    != string::npos) bitPho =  0; //bit0(lowest)
      else if (name.find("HLT_Photon30_v")                    != string::npos) bitPho =  1; 
      else if (name.find("HLT_Photon33_v")                    != string::npos) bitPho =  2; // 2017
      else if (name.find("HLT_Photon50_v")                    != string::npos) bitPho =  3; // 2017
      else if (name.find("HLT_Photon75_v")                    != string::npos) bitPho =  4; // 2017
      else if (name.find("HLT_Photon90_v")                    != string::npos) bitPho =  5; // 2017
      else if (name.find("HLT_Photon120_v")                   != string::npos) bitPho =  6; // 2017
      else if (name.find("HLT_Photon175_v")                   != string::npos) bitPho =  7; // 2017
      else if (name.find("HLT_Photon250_v")                   != string::npos) bitPho =  8; // 2017
      else if (name.find("HLT_Photon300_NoHE_v")              != string::npos) bitPho =  9; // 2017
      else if (name.find("HLT_Photon200_v")                   != string::npos) bitPho = 10; // 2017
      else if (name.find("HLT_Photon600_v")                   != string::npos) bitPho = 11; 
      else if (name.find("HLT_Photon165_HE10_v")              != string::npos) bitPho = 12; 
      else if (name.find("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15_v") != string::npos) bitPho = 13; // exist
      else if (name.find("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v")                             != string::npos) bitPho = 14; // used
      else if (name.find("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v")        != string::npos) bitPho = 15; // exist
      else if (name.find("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v")        != string::npos) bitPho = 16; // used
      else if (name.find("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v")         != string::npos) bitPho = 17; // used
      else if (name.find("HLT_Photon135_PFMET100_v")                          != string::npos) bitPho = 18; 
      else if (name.find("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v")  != string::npos) bitPho = 19;
      else if (name.find("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v")       != string::npos) bitPho = 20;
      else if (name.find("HLT_Photon90_CaloIdL_PFHT600_v")                    != string::npos) bitPho = 21;
      else if (name.find("HLT_DoublePhoton70_v")                              != string::npos) bitPho = 22; // 2017
      else if (name.find("HLT_DoublePhoton85_v")                              != string::npos) bitPho = 23; // 2017
      else if (name.find("HLT_Photon22_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 24;
      else if (name.find("HLT_TriplePhoton_20_20_20_CaloIdLV2_v")             != string::npos) bitPho = 25; //2017
      else if (name.find("HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL_v")      != string::npos) bitPho = 26; //2017
      else if (name.find("HLT_TriplePhoton_30_30_10_CaloIdLV2_v")             != string::npos) bitPho = 27; //2017
      else if (name.find("HLT_TriplePhoton_30_30_10_CaloIdLV2_R9IdVL_v")      != string::npos) bitPho = 28; //2017
      else if (name.find("HLT_TriplePhoton_35_35_5_CaloIdLV2_R9IdVL_v")       != string::npos) bitPho = 29; //2017
      else if (name.find("HLT_Photon50_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 30;
      else if (name.find("HLT_Photon75_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 31;
      else if (name.find("HLT_Photon90_R9Id90_HE10_IsoM_v")                   != string::npos) bitPho = 32;
      else if (name.find("HLT_Photon120_R9Id90_HE10_IsoM_v")                  != string::npos) bitPho = 33;
      else if (name.find("HLT_Photon165_R9Id90_HE10_IsoM_v")                  != string::npos) bitPho = 34;
      else if (name.find("HLT_ECALHT800_v")                                   != string::npos) bitPho = 35;
      
      // Jet triggers
      if      (name.find("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460_v")                    != string::npos) bitJet =  0;
      else if (name.find("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500_v")                    != string::npos) bitJet =  1; 
      else if (name.find("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200_v")                != string::npos) bitJet =  2; 
      else if (name.find("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240_v")                != string::npos) bitJet =  3; 
      else if (name.find("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v")                      != string::npos) bitJet =  4;
      else if (name.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v")                    != string::npos) bitJet =  5;
      else if (name.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v")                    != string::npos) bitJet =  6;
      else if (name.find("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != string::npos) bitJet =  7;
      else if (name.find("HLT_PFMET170_HBHECleaned_v")                                 != string::npos) bitJet =  8;
      else if (name.find("HLT_CaloJet500_NoJetID_v")                                   != string::npos) bitJet =  9;
      else if (name.find("HLT_PFJet40_v")                                              != string::npos) bitJet = 10; 
      else if (name.find("HLT_PFJet60_v")                                              != string::npos) bitJet = 11; 
      else if (name.find("HLT_PFJet80_v")                                              != string::npos) bitJet = 12; 
      else if (name.find("HLT_PFJet140_v")                                             != string::npos) bitJet = 13; 
      else if (name.find("HLT_PFJet200_v")                                             != string::npos) bitJet = 14; 
      else if (name.find("HLT_PFJet260_v")                                             != string::npos) bitJet = 15; 
      else if (name.find("HLT_PFJet320_v")                                             != string::npos) bitJet = 16; 
      else if (name.find("HLT_PFJet400_v")                                             != string::npos) bitJet = 17; 
      else if (name.find("HLT_PFJet450_v")                                             != string::npos) bitJet = 18;     
      else if (name.find("HLT_PFJet500_v")                                             != string::npos) bitJet = 19; 
      else if (name.find("HLT_AK8PFHT700_TrimR0p1PT0p3Mass50_v")                       != string::npos) bitJet = 20; 
      else if (name.find("HLT_AK8PFJet360_TrimMass30_v")                               != string::npos) bitJet = 21;
      else if (name.find("HLT_PFHT300_PFMET110_v")                                     != string::npos) bitJet = 22;
      else if (name.find("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067_v")   != string::npos) bitJet = 23;
      else if (name.find("HLT_PFMET170_HBHE_BeamHaloCleaned_v")                        != string::npos) bitJet = 24;
      else if (name.find("HLT_PFMET300_v")                                             != string::npos) bitJet = 25;
      else if (name.find("HLT_PFMET110_PFMHT110_IDTight_v")                            != string::npos) bitJet = 26;
      else if (name.find("HLT_PFMET120_PFMHT120_IDTight_v")                            != string::npos) bitJet = 27;
      else if (name.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v")                    != string::npos) bitJet = 28;
      else if (name.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v")                    != string::npos) bitJet = 29;
      else if (name.find("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != string::npos) bitJet = 30;
      else if (name.find("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != string::npos) bitJet = 31;
      else if (name.find("HLT_PFHT800_v")                                              != string::npos) bitJet = 32;
      else if (name.find("HLT_PFHT900_v")                                              != string::npos) bitJet = 33;
      else if (name.find("HLT_PFHT750_4JetPt50_v")                                     != string::npos) bitJet = 34;
      else if (name.find("HLT_PFHT750_4JetPt70_v")                                     != string::npos) bitJet = 35;
      else if (name.find("HLT_PFHT800_4JetPt50_v")                                     != string::npos) bitJet = 36;
    }

    // indicates prescaling and whether trigger was fired or not
    ULong64_t isPrescaled = (hltCfg.prescaleValue(prescaleSet, name)!=1) ? 1 : 0;
    ULong64_t isFired     = (trgResultsHandle->accept(i)) ? 1 : 0;
    ULong64_t isrejectedByHLTPS = (hltCfg.moduleType(hltCfg.moduleLabel(i,trgResultsHandle->index(i)))=="HLTPrescaler") ? 1: 0;

    if (bitEleMuX >= 0) {
      HLTEleMuX_            |= (isFired << bitEleMuX);
      HLTEleMuXIsPrescaled_ |= (isPrescaled << bitEleMuX);
    }

    if (bitPho >= 0) {
      HLTPho_            |= (isFired << bitPho);
      HLTPhoIsPrescaled_ |= (isPrescaled << bitPho);
      HLTPhoRejectedByPS_|= (isrejectedByHLTPS << bitPho);
    }

    if (bitJet >= 0) {
      HLTJet_            |= (isFired << bitJet);
      HLTJetIsPrescaled_ |= (isPrescaled << bitJet);
    }

    //if (name.find("HLT_PFJet450_v") == string::npos) 
    //cout<<"HLT : "<<i<<" "<<name<<" "<<isPrescaled<<" "<<isFired<<endl;

    if (!doGenParticles_) {
      /* to be fixed as this only works for 2017
      if      (name.find("HLT_Photon33_v")       != string::npos) phoPrescale_.insert(phoPrescale_.begin(),hltCfg.prescaleValue(prescaleSet, name));
      else if (name.find("HLT_Photon50_v")       != string::npos) phoPrescale_.insert(phoPrescale_.begin()+1,hltCfg.prescaleValue(prescaleSet, name));
      else if (name.find("HLT_Photon75_v")       != string::npos) phoPrescale_.insert(phoPrescale_.begin()+2,hltCfg.prescaleValue(prescaleSet, name));
      else if (name.find("HLT_Photon90_v")       != string::npos) phoPrescale_.insert(phoPrescale_.begin()+3,hltCfg.prescaleValue(prescaleSet, name));
      else if (name.find("HLT_Photon120_v")       != string::npos) phoPrescale_.insert(phoPrescale_.begin()+4,hltCfg.prescaleValue(prescaleSet, name));
      else if (name.find("HLT_Photon150_v")       != string::npos) phoPrescale_.insert(phoPrescale_.begin()+5,hltCfg.prescaleValue(prescaleSet, name));      
      */
      //else if (name.find("HLT_Photon120_v")      != string::npos) phoPrescale_.insert(phoPrescale_.begin()+6,hltCfg.prescaleValue(prescaleSet, name));
      //else if (name.find("HLT_Photon175_v")      != string::npos) phoPrescale_.insert(phoPrescale_.begin()+7,hltCfg.prescaleValue(prescaleSet, name));
      //else if (name.find("HLT_Photon250_NoHE_v") != string::npos) phoPrescale_.insert(phoPrescale_.end(),hltCfg.prescaleValue(prescaleSet, name));
    }

  }

}
