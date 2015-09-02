#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include <algorithm>

using namespace std;
using namespace reco;

// simple compare 
bool pairCompare(const std::pair<double, reco::PFCandidate>& firstItem,
		 const std::pair<double, reco::PFCandidate>& secondItem) {
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

  edm::Handle<reco::PFCandidateCollection> pfCandidateHandle;
  e.getByToken(pfAllParticles_, pfCandidateHandle);
  if (!pfCandidateHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no reco::PFCandidateCollection in event";
    return;
  }

  // =================================================================
  // PF HF Electron candidates
  vector<reco::PFCandidate> hfObjectsUnsorted;
  hfObjectsUnsorted.clear();
  for(unsigned int i=0; (unsigned int)i!=pfCandidateHandle->size(); ++i) {
    reco::PFCandidate pfcand = pfCandidateHandle->at(i);
    reco::PFCandidate::ParticleType candidateType = pfcand.particleId();
    if ( candidateType != reco::PFCandidate::h_HF &&
	 candidateType != reco::PFCandidate::egamma_HF ) continue;
    hfObjectsUnsorted.push_back(pfcand);
  }
  // 
  // sort them in energy
  vector<std::pair<double, reco::PFCandidate> > hfObjectPairs;
  for(unsigned int i = 0; (unsigned int)i!=hfObjectsUnsorted.size(); ++i) {
    std::pair<double, reco::PFCandidate> newCluster =
      make_pair(hfObjectsUnsorted[i].pt(), hfObjectsUnsorted[i]);
    hfObjectPairs.push_back(newCluster);
  }
  std::sort(hfObjectPairs.begin(), hfObjectPairs.end(), pairCompare);
  vector<reco::PFCandidate> hfObjects;
  hfObjects.clear();
  for(unsigned int i = 0; (unsigned int)i!=hfObjectPairs.size(); ++i) 
    hfObjects.push_back(hfObjectPairs[i].second);

  // Seeding is here ====================================================
  bool seeding = true;
  if ( seeding ) {
    vector<int> seeds;
    seeds.clear();
    for(unsigned int ipf = 0; ipf != (unsigned int)hfObjects.size(); ++ipf) {
      reco::PFCandidate object = hfObjects[ipf];
      // remove PFCandidates that are marked as spikes
      if ( object.ecalEnergy() == 0 && object.hcalEnergy() == 0 ) continue;
      
      TLorentzVector pf;
      pf.SetPtEtaPhiM(object.pt(), object.eta(), object.phi(), 0);
      if ( object.pt() < 5 ) continue; // require an energetic seed
      bool rejectSeed = false;
      for(int i = 0; i != (int)seeds.size(); ++i) {
	reco::PFCandidate oldSeed = hfObjects[seeds[i]];
	TLorentzVector pf2;
	pf2.SetPtEtaPhiM(oldSeed.pt(), 
			 oldSeed.eta(), 
			 oldSeed.phi(), 
			 0);
	if ( pf2.DeltaR(pf) < 0.3 ) { // no other seeds within 0.3 in DeltaR
	  rejectSeed = true;
	  break;
	}
      }
      if (!rejectSeed) {
	//cout << "This seed is selected!" << endl;
	seeds.push_back(ipf);
      }
    } // loop over seeds
    
    // Loop over seeds
    for(unsigned int iseed = 0; iseed != (unsigned int)seeds.size(); ++iseed) {
      TLorentzVector electronHFCandidate;
      int i = seeds[iseed];

      TLorentzVector seed;
      seed.SetPtEtaPhiM(hfObjects[i].pt(),
			hfObjects[i].eta(),
			hfObjects[i].phi(),
			0);
      electronHFCandidate = seed;
      
      double energyValue = hfObjects[i].energy();
      double isoValue = 0;
      double ecalEnergy = hfObjects[i].ecalEnergy();
      double hcalEnergy = hfObjects[i].hcalEnergy();
      // loop over the rest of PF objects
      for(int ipf = 0; ipf != (int)hfObjects.size(); ++ipf) {
	if ( ipf == i ) continue; // do not include the seed itself
	reco::PFCandidate pf = hfObjects[ipf];
	TLorentzVector cand;
	cand.SetPtEtaPhiM(pf.pt(), pf.eta(), pf.phi(), 0);
	if ( cand.DeltaR(seed) < 0.15 ) {
	  energyValue += cand.E();
	  ecalEnergy += pf.ecalEnergy();
	  hcalEnergy += pf.hcalEnergy();
	  electronHFCandidate += cand;
	}
	else if ( cand.DeltaR(seed) < 0.3 ) {
	  if ( pf.ecalEnergy() == 0 &&
	       pf.hcalEnergy() == 0 ) isoValue += 10000.0;
	  isoValue += cand.E();
	}
      } 
      
      // Now, selection criteria
      bool passEoH = false;
      bool passIso = false;
      if ( hcalEnergy == 0 || ecalEnergy/hcalEnergy > 0.4 ) passEoH = true;
      if ( isoValue/energyValue < 0.4 ) passIso = true;
      
      if ( !(passEoH && passIso) ) continue;
      
      // object passed criteria, filling the tree
      pfHFEn_.push_back(electronHFCandidate.E());
      pfHFECALEn_.push_back(ecalEnergy);
      pfHFHCALEn_.push_back(hcalEnergy);
      pfHFPt_    .push_back(electronHFCandidate.Pt());
      pfHFPhi_   .push_back(electronHFCandidate.Phi());
      pfHFEta_   .push_back(electronHFCandidate.Eta());
      pfHFIso_   .push_back(isoValue/energyValue);
      ++npfHF_;
    }
  } else { // no seeding, save every single PFCandidate
    for(unsigned int ipf = 0; ipf != (unsigned int)hfObjects.size(); ++ipf) {
      reco::PFCandidate pf = hfObjects[ipf];
      pfHFEn_.push_back(pf.energy());
      pfHFECALEn_.push_back(pf.ecalEnergy());
      pfHFHCALEn_.push_back(pf.hcalEnergy());
      pfHFPt_    .push_back(pf.pt());
      pfHFPhi_   .push_back(pf.phi());
      pfHFEta_   .push_back(pf.eta());
      pfHFIso_   .push_back(-1);
      ++npfHF_;
    } // no seeding
  } // loop over HF PF objects    

}
