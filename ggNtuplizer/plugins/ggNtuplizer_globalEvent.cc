#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t     run_;
Long64_t  event_;
Int_t     lumis_;
Bool_t    isData_;
Int_t     nVtx_;
Int_t     nTrksPV_;
float     vtx_;
float     vty_;
float     vtz_;
float     rho_;
ULong64_t HLT_;
ULong64_t HLTIsPrescaled_;
ULong64_t HLT50ns_;
ULong64_t HLTIsPrescaled50ns_;
float     genMET_;
float     genMETPhi_;
float     pfMET_;
float     pfMETPhi_;
float     pfMETsumEt_;
float     pfMETmEtSig_;
float     pfMETSig_;

void ggNtuplizer::branchesGlobalEvent(TTree* tree) {

  tree->Branch("run",     &run_);
  tree->Branch("event",   &event_);
  tree->Branch("lumis",   &lumis_);
  tree->Branch("isData",  &isData_);
  tree->Branch("nVtx",    &nVtx_);
  tree->Branch("nTrksPV", &nTrksPV_);
  tree->Branch("vtx",     &vtx_); 
  tree->Branch("vty",     &vty_); 
  tree->Branch("vtz",     &vtz_); 
  tree->Branch("rho",     &rho_);
  tree->Branch("HLT",     &HLT_);
  tree->Branch("HLTIsPrescaled", &HLTIsPrescaled_);
  tree->Branch("HLT50ns",     &HLT50ns_);
  tree->Branch("HLTIsPrescaled50ns", &HLTIsPrescaled50ns_);
  if (doGenParticles_) {
    tree->Branch("genMET",      &genMET_);
    tree->Branch("genMETPhi",   &genMETPhi_);
  }
  tree->Branch("pfMET",       &pfMET_);
  tree->Branch("pfMETPhi",    &pfMETPhi_);
  tree->Branch("pfMETsumEt",  &pfMETsumEt_);
  tree->Branch("pfMETmEtSig", &pfMETmEtSig_);
  tree->Branch("pfMETSig",    &pfMETSig_);
}

void ggNtuplizer::fillGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();
  rho_    = *(rhoHandle.product());

  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  
  nVtx_ = -1;
  if (vtxHandle.isValid()) {
    nVtx_ = 0;
    
    for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
      bool isFake = isAOD_ ? v->isFake() : (v->chi2() == 0 && v->ndof() == 0);
      if (!isFake) {
        if (nVtx_ == 0) {
	  nTrksPV_ = v->nTracks();
	  vtx_     = v->x();
	  vty_     = v->y();
	  vtz_     = v->z();
	}
        nVtx_++;
      }
    }
  } else
    edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";

  // HLT treatment
  HLT_ = 0;

  edm::Handle<edm::TriggerResults> trgResultsHandle;
  e.getByToken(trgResultsLabel_, trgResultsHandle);

  bool cfg_changed = true;
  HLTConfigProvider hltCfg;
  hltCfg.init(e.getRun(), es, trgResultsProcess_, cfg_changed);

  const edm::TriggerNames &trgNames = e.triggerNames(*trgResultsHandle);

  for (size_t i = 0; i < trgNames.size(); ++i) {
    const string &name = trgNames.triggerName(i);

    // HLT name => bit correspondence
    int bit = -1;
    if      (name.find("HLT_Physics_v")                    != string::npos) bit = 0;  // bit 0 (lowest)
    else if (name.find("HLT_L1SingleMu3p5_WideWindow_v" )  != string::npos) bit = 1;  // bit 1
    else if (name.find("HLT_L1SingleMu3p5_v" )             != string::npos) bit = 2;
    else if (name.find("HLT_L1SingleMuOpen_WideWindow_v" ) != string::npos) bit = 3;
    else if (name.find("HLT_L1SingleMuOpen_v" )            != string::npos) bit = 4;
    else if (name.find("HLT_L1SingleEG5_WideWindow_v")     != string::npos) bit = 5;
    else if (name.find("HLT_L1SingleEG5_v")                != string::npos) bit = 6;
    else if (name.find("HLT_L1SingleEG20_WideWindow_v")    != string::npos) bit = 7;
    else if (name.find("HLT_L1SingleEG20_v")               != string::npos) bit = 8;
    else if (name.find("HLT_L1SingleJet36_WideWindow_v")   != string::npos) bit = 9;
    else if (name.find("HLT_L1SingleJet36_v")              != string::npos) bit = 10;
    else if (name.find("HLT_L1SingleJet68_WideWindow_v")   != string::npos) bit = 11;
    else if (name.find("HLT_L1SingleJet68_v")              != string::npos) bit = 12;
    else if (name.find("HLT_ZeroBias_")                    != string::npos) bit = 13;

    //////////////////////////////triggers for 50ns////////////////////////////////////////
    
    ///muon triggers
    Long64_t bit50ns = -1;
    if      (name.find("HLT_Mu50_v")                    != string::npos) bit50ns = 0;  // bit 0 (lowest)
    else if (name.find("HLT_Mu45_eta2p1_v" )  != string::npos) bit50ns = 1; 
    else if (name.find("HLT_IsoMu20_v" )  != string::npos) bit50ns = 2; 
    else if (name.find("HLT_IsoMu17_eta2p1_v" )  != string::npos) bit50ns = 3; 
    else if (name.find("HLT_IsoTkMu20_v" )  != string::npos) bit50ns = 4; 
    else if (name.find("HLT_IsoTkMu20_eta2p1_v" )  != string::npos) bit50ns = 5; 
    else if (name.find("HLT_Mu27_TkMu8_v" )  != string::npos) bit50ns = 6; 
    else if (name.find("HLT_Mu17_TrkIsoVVL_v" )  != string::npos) bit50ns = 7; 
    else if (name.find("HLT_Mu17_Photon30_CaloIdL_L1ISO_v" )  != string::npos) bit50ns = 8; 
    else if (name.find("HLT_DoubleIsoMu17_eta2p1_v" )  != string::npos) bit50ns = 9; 
    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v" )  != string::npos) bit50ns = 10; 
    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v" )  != string::npos) bit50ns = 11; 
    
    ///electron triggers
    else if (name.find("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v" )  != string::npos) bit50ns = 12; 
    else if (name.find("HLT_Ele105_CaloIdVT_GsfTrkIdT_v" )  != string::npos) bit50ns = 13; 
    else if (name.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v" )  != string::npos) bit50ns = 14; 
    else if (name.find("HLT_Ele27_eta2p1_WPTight_Gsf_v" )  != string::npos) bit50ns = 15; 
    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v" )  != string::npos) bit50ns = 16; 
    else if (name.find("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v" )  != string::npos) bit50ns = 17; 
    else if (name.find("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v" )  != string::npos) bit50ns = 18; 
    else if (name.find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v" )  != string::npos) bit50ns = 19; 
    else if (name.find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v" )  != string::npos) bit50ns = 20; 

    ///photon triggers
    else if (name.find("HLT_Photon175_v")                    != string::npos) bit50ns = 21;  
    else if (name.find("HLT_Photon250_NoHE_v" )  != string::npos) bit50ns = 22; 
    else if (name.find("HLT_Photon165_HE10_v" )  != string::npos) bit50ns = 23; 
    else if (name.find("HLT_Photon165_R9Id90_HE10_IsoM_v" )  != string::npos) bit50ns = 24; 
    else if (name.find("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v" )  != string::npos) bit50ns = 25; 
    else if (name.find("HLT_DoublePhoton85_v" )  != string::npos) bit50ns = 26; 
    else if (name.find("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15_v" )  != string::npos) bit50ns = 27; 
    else if (name.find("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass95_v" )  != string::npos) bit50ns = 28; 
    else if (name.find("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70_v" )  != string::npos) bit50ns = 29; 
    else if (name.find("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v" )  != string::npos) bit50ns = 30; 
    else if (name.find("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v" )  != string::npos) bit50ns = 31; 
    else if (name.find("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55_v" )  != string::npos) bit50ns = 32;
    // tau triggers
    else if (name.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v" )  != string::npos) bit50ns = 33; 
    else if (name.find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v")  != string::npos) bit50ns = 34;
    // Tri-muon triggers
    else if (name.find("HLT_TripleMu_12_10_5_v" )  != string::npos) bit50ns = 35; 
    // cross triggers
    else if (name.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v" )  != string::npos) bit50ns = 36; 
    else if (name.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v" )  != string::npos) bit50ns = 37; 
    else if (name.find("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v" )  != string::npos) bit50ns = 38;
    else if (name.find("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v" )  != string::npos) bit50ns = 39;
    else if (name.find("HLT_Mu17_Photon30_CaloIdL_L1ISO_v" )  != string::npos) bit50ns = 40;
    else if (name.find("HLT_Mu17_Photon35_CaloIdL_L1ISO_v" )  != string::npos) bit50ns = 41;

    // indicates prescaling and whether trigger was fired or not
    ULong64_t isPrescaled = (hltCfg.prescaleValue(e, es, name) > 0 ? 1 : 0);
    ULong64_t isFired = (trgResultsHandle->accept(i) && isPrescaled) ? 1 : 0;

    if (bit >= 0) {
      HLT_            |= (isFired << bit);
      HLTIsPrescaled_ |= (isPrescaled << bit);
    }

    if (bit50ns >= 0) {
      HLT50ns_            |= (isFired << bit50ns);
      HLTIsPrescaled50ns_ |= (isPrescaled << bit50ns);
    }

    //////////////////////////end of HLT 50ns menu////////////////////////////
    
    

  }//for (size_t i = 0; i < trgNames.size(); ++i)

  edm::Handle<edm::View<pat::MET> > pfMETHandle;
  e.getByToken(pfMETlabel_, pfMETHandle);

  genMET_      = -99;
  genMETPhi_   = -99;
  pfMET_       = -99;
  pfMETPhi_    = -99;
  pfMETsumEt_  = -99;
  pfMETmEtSig_ = -99;
  pfMETSig_    = -99;

  if (pfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET = &(pfMETHandle->front());

    pfMET_       = pfMET->et();
    pfMETPhi_    = pfMET->phi();
    pfMETsumEt_  = pfMET->sumEt();
    pfMETmEtSig_ = (pfMET->mEtSig() < 1.e10) ? pfMET->mEtSig() : 0;
    pfMETSig_    = (pfMET->significance() < 1.e10) ? pfMET->significance() : 0;;

    if (!isData_) {
      genMET_    = pfMET->genMET()->et();
      genMETPhi_ = pfMET->genMET()->phi();
    }

  } else
    edm::LogWarning("ggNtuplizer") << "MET info unavailable";

}
