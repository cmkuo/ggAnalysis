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
float     rho_;
ULong64_t HLT_;
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
  tree->Branch("rho",     &rho_);
  tree->Branch("HLT",     &HLT_);

  if (!isData_) {
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
        if (nVtx_ == 0) nTrksPV_ = v->nTracks();
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

    // indicates whether trigger was fired or not
    ULong64_t bit = (trgResultsHandle->accept(i) && hltCfg.prescaleValue(e, es, name) != 0) ? 1 : 0;

    // HLT name => bit correspondence
    if      (!name.compare("HLT_Physics_v")                              ) HLT_ |= (bit<<  0);  // bit 0 (lowest)
    else if (!name.compare("HLT_L1SingleMu3p5_WideWindow_v" )            ) HLT_ |= (bit<<  1);  // bit 1
    else if (!name.compare("HLT_L1SingleMu3p5_v" )                       ) HLT_ |= (bit<<  2);  
    else if (!name.compare("HLT_L1SingleMuOpen_WideWindow_v" )           ) HLT_ |= (bit<<  3);  
    else if (!name.compare("HLT_L1SingleMuOpen_v" )                      ) HLT_ |= (bit<<  4);  
    else if (!name.compare("HLT_L1SingleEG5_WideWindow_v")               ) HLT_ |= (bit<<  5); 
    else if (!name.compare("HLT_L1SingleEG5_v")                          ) HLT_ |= (bit<<  6);  
    else if (!name.compare("HLT_L1SingleEG20_WideWindow_v")              ) HLT_ |= (bit<<  7);  
    else if (!name.compare("HLT_L1SingleEG20_v")                         ) HLT_ |= (bit<<  8);  
    else if (!name.compare("HLT_L1SingleJet36_WideWindow_v")             ) HLT_ |= (bit<<  9);
    else if (!name.compare("HLT_L1SingleJet36_v")                        ) HLT_ |= (bit<< 10);
    else if (!name.compare("HLT_L1SingleJet68_WideWindow_v")             ) HLT_ |= (bit<< 11);
    else if (!name.compare("HLT_L1SingleJet68_v")                        ) HLT_ |= (bit<< 12);
    else if (!name.compare("HLT_ZeroBias_")                              ) HLT_ |= (bit<< 13);
  }

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
