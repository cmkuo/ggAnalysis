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

void ggNtuplizer::branchesPFPhotons(TTree* tree) {
  
  tree->Branch("npfPho",      &npfPho_);
  tree->Branch("pfPhoEt",     &pfPhoEt_);
  tree->Branch("pfPhoEta",    &pfPhoEta_);
  tree->Branch("pfPhoPhi",    &pfPhoPhi_);
  tree->Branch("pfPhoChIso",  &pfPhoChIso_);
  tree->Branch("pfPhoPhoIso", &pfPhoPhoIso_);
  tree->Branch("pfPhoNeuIso", &pfPhoNeuIso_);
}

void ggNtuplizer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  pfPhoEt_    .clear();
  pfPhoEta_   .clear();
  pfPhoPhi_   .clear();
  pfPhoChIso_ .clear();
  pfPhoPhoIso_.clear();
  pfPhoNeuIso_.clear();

  npfPho_ = 0;

  edm::Handle<pat::PackedCandidateCollection> cands;
  e.getByToken(pckPFCandidateCollection_, cands);

  for (unsigned int i = 0; i< cands->size(); ++i) {
    const pat::PackedCandidate &pf = (*cands)[i];

    if (pf.pdgId() == 22 && pf.pt() > 1.5 && fabs(pf.eta()) < 2.5) {

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

      pfPhoEt_    .push_back(pf.pt());
      pfPhoEta_   .push_back(pf.eta());
      pfPhoPhi_   .push_back(pf.phi());
      pfPhoChIso_ .push_back(chiso_);
      pfPhoPhoIso_.push_back(phoiso_);
      pfPhoNeuIso_.push_back(neuiso_);

      npfPho_++;
    }

  }

  
}

