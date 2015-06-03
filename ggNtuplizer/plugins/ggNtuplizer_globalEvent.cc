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

float    genMET_;
float    genMETPhi_;
float    pfMET_;
float    pfMETPhi_;
float    pfMETsumEt_;
float    pfMETmEtSig_;
float    pfMETSig_;

void ggNtuplizer::branchesGlobalEvent(TTree* tree) {

  tree->Branch("run",     &run_);
  tree->Branch("event",   &event_);
  tree->Branch("lumis",   &lumis_);
  tree->Branch("isData",  &isData_);
  tree->Branch("nVtx",    &nVtx_);
  tree->Branch("nTrksPV", &nTrksPV_);
  tree->Branch("rho",     &rho_);
  tree->Branch("HLT",     &HLT_);

  tree->Branch("genMET",      &genMET_);
  tree->Branch("genMETPhi",   &genMETPhi_);
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
    if      (!name.compare("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v")      ) HLT_ |= (bit<< 0);  // bit 0 (lowest)
    else if (!name.compare("HLT_Photon26_Photon18_v")                                                               ) HLT_ |= (bit<< 1);  // bit 1
    else if (!name.compare("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v")                          ) HLT_ |= (bit<< 2);  // ...
    else if (!name.compare("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v")                                  ) HLT_ |= (bit<< 3);
    else if (!name.compare("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v")                                      ) HLT_ |= (bit<< 4);
    else if (!name.compare("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v")                                  ) HLT_ |= (bit<< 5);
    else if (!name.compare("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v")                                          ) HLT_ |= (bit<< 6);
    else if (!name.compare("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v")      ) HLT_ |= (bit<< 7);
    else if (!name.compare("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v")             ) HLT_ |= (bit<< 8);
    else if (!name.compare("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon10_R9Id85_OR_CaloId10_Iso50_Mass80_v")      ) HLT_ |= (bit<< 9);
    else if (!name.compare("HLT_Photon36_Photon22_v")                                                               ) HLT_ |= (bit<<10);
    else if (!name.compare("HLT_PAPhoton15_TightCaloIdVL_v")                                                        ) HLT_ |= (bit<<11);
    else if (!name.compare("HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_v")                                        ) HLT_ |= (bit<<12);
    else if (!name.compare("HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned_v")                 ) HLT_ |= (bit<<13);
    else if (!name.compare("HLT_Photon30_v")                                                                        ) HLT_ |= (bit<<14);
    else if (!name.compare("HLT_Photon30_CaloIdVL_v")                                                               ) HLT_ |= (bit<<15);
    else if (!name.compare("HLT_Photon50_CaloIdVL_IsoL_v")                                                          ) HLT_ |= (bit<<16);
    else if (!name.compare("HLT_Photon50_CaloIdVL_v")                                                               ) HLT_ |= (bit<<17);
    else if (!name.compare("HLT_Photon75_CaloIdVL_v")                                                               ) HLT_ |= (bit<<18);
    else if (!name.compare("HLT_Photon90_CaloIdVL_v")                                                               ) HLT_ |= (bit<<19);
    else if (!name.compare("HLT_Photon135_v")                                                                       ) HLT_ |= (bit<<20);
    else if (!name.compare("HLT_Photon150_v")                                                                       ) HLT_ |= (bit<<21);
    else if (!name.compare("HLT_Photon160_v")                                                                       ) HLT_ |= (bit<<22);
    else if (!name.compare("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v")                                 ) HLT_ |= (bit<<23);
    else if (!name.compare("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v")                                         ) HLT_ |= (bit<<24);
    else if (!name.compare("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v")                                         ) HLT_ |= (bit<<25);
    else if (!name.compare("HLT_Photon36_R9Id85_Photon22_R9Id85_v")                                                 ) HLT_ |= (bit<<26);
    else if (!name.compare("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) HLT_ |= (bit<<27);
    else if (!name.compare("HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele15_CaloIdT_CaloIsoVL_trackless_v")      ) HLT_ |= (bit<<28);
    else if (!name.compare("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v")                               ) HLT_ |= (bit<<29);
    else if (!name.compare("HLT_Ele27_WP80_v")                                                                      ) HLT_ |= (bit<<30);
    else if (!name.compare("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v")) HLT_ |= (bit<<31);
    else if (!name.compare("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v")                             ) HLT_ |= (bit<<32);
    else if (!name.compare("HLT_Mu22_TkMu22_v")                                                                     ) HLT_ |= (bit<<33);
    else if (!name.compare("HLT_Mu17_TkMu8_v")                                                                      ) HLT_ |= (bit<<34);
    else if (!name.compare("HLT_Mu17_Mu8_v")                                                                        ) HLT_ |= (bit<<35);
    else if (!name.compare("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")                                    ) HLT_ |= (bit<<36);
    else if (!name.compare("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")                                    ) HLT_ |= (bit<<37);
    else if (!name.compare("HLT_Mu22_TkMu8_v")                                                                      ) HLT_ |= (bit<<38);
    else if (!name.compare("HLT_Mu22_Mu8_v")                                                                        ) HLT_ |= (bit<<39);
    else if (!name.compare("HLT_Mu22_Photon22_CaloIdL_v")                                                           ) HLT_ |= (bit<<40);
    else if (!name.compare("HLT_IsoMu24_eta2p1_v")                                                                  ) HLT_ |= (bit<<41);
    else if (!name.compare("HLT_IsoMu24_v")                                                                         ) HLT_ |= (bit<<42);
    else if (!name.compare("HLT_PAMu12_v")                                                                          ) HLT_ |= (bit<<43);
    else if (!name.compare("HLT_PAL1DoubleMuOpen_v")                                                                ) HLT_ |= (bit<<44);
    else if (!name.compare("HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF_v")                                              ) HLT_ |= (bit<<45);
    else if (!name.compare("HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v")                                              ) HLT_ |= (bit<<46);
    else if (!name.compare("HLT_PFJet40_v")                                                                         ) HLT_ |= (bit<<47);
    else if (!name.compare("HLT_PFJet80_v")                                                                         ) HLT_ |= (bit<<48);
    else if (!name.compare("HLT_PFJet140_v")                                                                        ) HLT_ |= (bit<<49);
    else if (!name.compare("HLT_PFJet200_v")                                                                        ) HLT_ |= (bit<<50);
    else if (!name.compare("HLT_PFJet260_v")                                                                        ) HLT_ |= (bit<<51);
    else if (!name.compare("HLT_Physics_v1")                                                                        ) HLT_ |= (bit<<52);
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
