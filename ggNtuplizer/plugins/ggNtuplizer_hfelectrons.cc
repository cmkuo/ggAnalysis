#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include <algorithm>

using namespace std;
using namespace reco;

// simple compare 
bool pairCompare(const std::pair<double, edm::Ptr<reco::Candidate> >& firstItem,
		 const std::pair<double, edm::Ptr<reco::Candidate> >& secondItem) {
  return firstItem.first > secondItem.first;
}

// pf particles
Int_t npfHF_;
vector<float> pfHFEn_;
vector<float> pfHFECALEn_;
vector<float> pfHFHCALEn_;
vector<float> pfHFPt_;
vector<float> pfHFPhi_;
vector<float> pfHFEta_;
vector<float> pfHFIso_;

void ggNtuplizer::branchesHFElectrons(TTree* tree) {

  tree->Branch("npfHF",      &npfHF_);
  tree->Branch("pfHFEn",     &pfHFEn_);
  tree->Branch("pfHFECALEn", &pfHFECALEn_);
  tree->Branch("pfHFHCALEn", &pfHFHCALEn_);
  tree->Branch("pfHFPt",     &pfHFPt_);
  tree->Branch("pfHFEta",    &pfHFEta_);  
  tree->Branch("pfHFPhi",    &pfHFPhi_);  
  tree->Branch("pfHFIso",    &pfHFIso_);
}

void ggNtuplizer::fillHFElectrons(const edm::Event &e) {
    
  // cleanup from previous execution
  pfHFEn_    .clear();
  pfHFECALEn_.clear();
  pfHFHCALEn_.clear();
  pfHFPt_    .clear();
  pfHFPhi_   .clear();
  pfHFEta_   .clear();
  pfHFIso_   .clear();
  npfHF_ = 0;

  // =================================================================
  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsAK4Label_, jetHandle);
  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
    return;
  }

  // Loop over jets
  for (edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
    if ( fabs(iJet->eta()) < 2.9 ) continue;
    TLorentzVector pfJet;
    pfJet.SetPtEtaPhiM(iJet->pt(), iJet->eta(), iJet->phi(), 0);

    vector<edm::Ptr<reco::Candidate> > constituentsUnsorted = iJet->getJetConstituents();
    vector<edm::Ptr<reco::Candidate> > constituents;
    constituents.clear();
    vector<std::pair<double, edm::Ptr<reco::Candidate> > > pfEn;
    pfEn.clear();
    for(auto& cand : constituentsUnsorted ) {
      std::pair<double, edm::Ptr<reco::Candidate> > newPair = make_pair(cand->pt(), cand);
      pfEn.push_back(newPair);
    }
    sort(pfEn.begin(), pfEn.end(), pairCompare);
    for(auto& cand : pfEn ) 
      constituents.push_back(cand.second);

    //int counter = 0;
    vector<edm::Ptr<pat::PackedCandidate> > pfs;
    pfs.clear();
    for(auto& cand : constituents ) {
      edm::Ptr<pat::PackedCandidate> pkcand(cand);
      pfs.push_back(pkcand);
      //cout << ++counter << " " << pkcand->pt() << " " << pkcand->eta() << " " << pkcand->hcalFraction() << endl;
    }

    // seed is the most energetic cluster in the jet
    //if ( pfs[0]->ecalEnergy() == 0 && pfs[0]->hcalEnergy() == 0 ) continue; // skip the spike

    TLorentzVector core;
    core.SetPtEtaPhiM(pfs[0]->pt(), pfs[0]->eta(), pfs[0]->phi(), 0);    
    TLorentzVector isolation;
    isolation.SetPtEtaPhiM(0,0,0,0);
    double ecalEnergy = (1.0 - pfs[0]->hcalFraction())*pfs[0]->energy();
    double hcalEnergy = pfs[0]->hcalFraction()*pfs[0]->energy();

    for(int i = 1; i != (int)pfs.size(); ++i) {
      pat::PackedCandidate pfCandidate =*(pfs[i]);
      TLorentzVector pf;
      pf.SetPtEtaPhiM(pfCandidate.pt(),
		      pfCandidate.eta(),
		      pfCandidate.phi(),
		      0);
      double dr = pf.DeltaR(core);
      if ( dr < 0.15 ) {
	core += pf;
	ecalEnergy += pfCandidate.energy()*(1.0 - pfCandidate.hcalFraction());
	hcalEnergy += pfCandidate.energy()*pfCandidate.hcalFraction();
      }
      if ( dr >= 0.15 && dr < 0.3 ) {
	isolation += pf;
      }
    } // Looped over constituents

    // Selection criteria
    //if ( spike ) continue; // should not have any PFCandidates identidied as spike
    //if ( core.Pt() < 5.0 ) continue; // sufficiently high pT
    if ( pfs[0]->pt() < 5.0 ) continue; // require the seed PFCandidate to be energetic enough
    
    if ( ecalEnergy == 0 ) continue; // must have ecal energy to be called an electron/photon in HF
    if ( hcalEnergy/ecalEnergy > 2.25 ) continue; // equivalent to ecal/hcal < 0.4 requirement
    if ( isolation.E()/core.E() > 0.55 ) continue;
    
    ++npfHF_;
    pfHFEn_.push_back(core.E());
    pfHFECALEn_.push_back(ecalEnergy);
    pfHFHCALEn_.push_back(hcalEnergy);
    pfHFPt_.push_back(core.Pt());
    pfHFPhi_.push_back(core.Phi());
    pfHFEta_.push_back(core.Eta());
    pfHFIso_.push_back(isolation.E());
  }

}
