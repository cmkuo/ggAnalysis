#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          npfPho_;
vector<float>  pfPhoEt_;
vector<float>  pfPhoEta_;
vector<float>  pfPhoPhi_;
vector<float>  pfPhoChIso_;
vector<float>  pfPhoPhoIso_;
vector<float>  pfPhoNeuIso_;
vector<int>    pfPhoEleIndex_;
vector<int>    pfPhoPhoIndex_;

void ggNtuplizer::branchesPFPhotons(TTree* tree) {
  
  tree->Branch("npfPho",        &npfPho_);
  tree->Branch("pfPhoEt",       &pfPhoEt_);
  tree->Branch("pfPhoEta",      &pfPhoEta_);
  tree->Branch("pfPhoPhi",      &pfPhoPhi_);
  tree->Branch("pfPhoChIso",    &pfPhoChIso_);
  tree->Branch("pfPhoPhoIso",   &pfPhoPhoIso_);
  tree->Branch("pfPhoNeuIso",   &pfPhoNeuIso_);
  tree->Branch("pfPhoEleIndex", &pfPhoEleIndex_);
  tree->Branch("pfPhoPhoIndex", &pfPhoPhoIndex_);
}

void ggNtuplizer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  pfPhoEt_      .clear();
  pfPhoEta_     .clear();
  pfPhoPhi_     .clear();
  pfPhoChIso_   .clear();
  pfPhoPhoIso_  .clear();
  pfPhoNeuIso_  .clear();
  pfPhoEleIndex_.clear();
  pfPhoPhoIndex_.clear();

  npfPho_ = 0;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);
  
  edm::Handle<pat::PackedCandidateCollection> cands;
  e.getByToken(pckPFCandidateCollection_, cands);

  for (unsigned int i = 0; i< cands->size(); ++i) {
    const pat::PackedCandidate &pf = (*cands)[i];

    if (pf.pdgId() == 22 && pf.pt() > 1.5 && fabs(pf.eta()) < 2.5) {

      int nEle_     = 0;
      int eleIndex_ = -1;
      for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {
	for (const edm::Ref<pat::PackedCandidateCollection> & refEle : iEle->associatedPackedPFCandidates()) {
	  if (reco::deltaR(pf.p4(), refEle->p4()) < 0.0001) eleIndex_ = nEle_;
	}
	nEle_++;
      }
      
      int nPho_     = 0;
      int phoIndex_ = -1;
      for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {
	for (const edm::Ref<pat::PackedCandidateCollection> & refPho : iPho->associatedPackedPFCandidates()) {
	  if (reco::deltaR(pf.p4(), refPho->p4()) < 0.0001) phoIndex_ = nPho_;
	}
	nPho_++;
      }

      float chiso_  = 0;
      float phoiso_ = 0;
      float neuiso_ = 0;
      float dR_     = 0;

      for (unsigned int j = 0; j< cands->size(); ++j) {
	const pat::PackedCandidate &pff = (*cands)[j];

	if (abs(pff.pdgId()) == 211 && pff.pt() > 0.2) {
	  dR_ = reco::deltaR(pf.p4(), pff.p4());
	  if (dR_ < 0.3 && dR_ > 0.0001) chiso_ += pff.pt(); 	    
	}

	if (abs(pff.pdgId()) == 22  && pff.pt() > 0.5) {
	  dR_ = reco::deltaR(pf.p4(), pff.p4());
	  if (dR_ < 0.3 && dR_ > 0.01) phoiso_ += pff.pt(); 	    
	}

	if (abs(pff.pdgId()) == 130 && pff.pt() > 0.5) {
	  dR_ = reco::deltaR(pf.p4(), pff.p4());
	  if (dR_ < 0.3 && dR_ > 0.01) neuiso_ += pff.pt(); 	    
	}

      }

      pfPhoEt_      .push_back(pf.pt());
      pfPhoEta_     .push_back(pf.eta());
      pfPhoPhi_     .push_back(pf.phi());
      pfPhoChIso_   .push_back(chiso_);
      pfPhoPhoIso_  .push_back(phoiso_);
      pfPhoNeuIso_  .push_back(neuiso_);
      pfPhoEleIndex_.push_back(eleIndex_);
      pfPhoPhoIndex_.push_back(phoIndex_);

      npfPho_++;
    }

  }

  
}

