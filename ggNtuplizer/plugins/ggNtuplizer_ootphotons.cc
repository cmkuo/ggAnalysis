#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          nOOTPho_;
vector<float>  ootPhoE_;
vector<float>  ootPhoEt_;
vector<float>  ootPhoEta_;
vector<float>  ootPhoPhi_;
vector<float>  ootPhoSCEta_;
vector<float>  ootPhoSCPhi_;
vector<float>  ootPhoSCEtaWidth_;
vector<float>  ootPhoSCPhiWidth_;
vector<float>  ootPhoR9_;
vector<float>  ootPhoHoverE_;
vector<float>  ootPhoSigmaIEtaIEtaFull5x5_;
vector<float>  ootPhoR9Full5x5_;
vector<float>  ootPhoSeedTime_;
vector<float>  ootPhoSeedEnergy_;
vector<float>  ootPhoMIPTotEnergy_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesOOTPhotons(TTree* tree) {
  
  tree->Branch("nOOTPho",                     &nOOTPho_);
  tree->Branch("ootPhoE",                     &ootPhoE_);
  tree->Branch("ootPhoEt",                    &ootPhoEt_);
  tree->Branch("ootPhoEta",                   &ootPhoEta_);
  tree->Branch("ootPhoPhi",                   &ootPhoPhi_);
  tree->Branch("ootPhoSCEta",                 &ootPhoSCEta_);
  tree->Branch("ootPhoSCPhi",                 &ootPhoSCPhi_);
  tree->Branch("ootPhoSCEtaWidth",            &ootPhoSCEtaWidth_);
  tree->Branch("ootPhoSCPhiWidth",            &ootPhoSCPhiWidth_);
  tree->Branch("ootPhoR9",                    &ootPhoR9_);
  tree->Branch("ootPhoHoverE",                &ootPhoHoverE_);
  tree->Branch("ootPhoSigmaIEtaIEtaFull5x5",  &ootPhoSigmaIEtaIEtaFull5x5_);
  tree->Branch("ootPhoR9Full5x5",             &ootPhoR9Full5x5_);
  tree->Branch("ootPhoSeedTime",              &ootPhoSeedTime_);
  tree->Branch("ootPhoSeedEnergy",            &ootPhoSeedEnergy_);
  tree->Branch("ootPhoMIPTotEnergy",          &ootPhoMIPTotEnergy_);

}

void ggNtuplizer::fillOOTPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  ootPhoE_                   .clear();
  ootPhoEt_                  .clear();
  ootPhoEta_                 .clear();
  ootPhoPhi_                 .clear();
  ootPhoSCEta_               .clear();
  ootPhoSCPhi_               .clear();
  ootPhoSCEtaWidth_          .clear();
  ootPhoSCPhiWidth_          .clear();
  ootPhoR9_                  .clear();
  ootPhoHoverE_              .clear();
  ootPhoSigmaIEtaIEtaFull5x5_.clear();
  ootPhoR9Full5x5_           .clear();
  ootPhoSeedTime_            .clear();
  ootPhoSeedEnergy_          .clear();
  ootPhoMIPTotEnergy_        .clear();

  nOOTPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > ootPhotonHandle;
  e.getByToken(ootPhotonCollection_, ootPhotonHandle);

  if (!ootPhotonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no OOT pat::Photons in event";
    return;
  }

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Photon>::const_iterator iOOTPho = ootPhotonHandle->begin(); iOOTPho != ootPhotonHandle->end(); ++iOOTPho) {

    ootPhoE_                     .push_back(iOOTPho->energy());
    ootPhoEt_                    .push_back(iOOTPho->et());
    ootPhoEta_                   .push_back(iOOTPho->eta());
    ootPhoPhi_                   .push_back(iOOTPho->phi());
    ootPhoSCEta_                 .push_back((*iOOTPho).superCluster()->eta());
    ootPhoSCPhi_                 .push_back((*iOOTPho).superCluster()->phi());
    ootPhoSCEtaWidth_            .push_back((*iOOTPho).superCluster()->etaWidth());
    ootPhoSCPhiWidth_            .push_back((*iOOTPho).superCluster()->phiWidth());
    ootPhoR9_                    .push_back(iOOTPho->r9());
    ootPhoHoverE_                .push_back(iOOTPho->hadTowOverEm());

    DetId seed = (iOOTPho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());

    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      ootPhoSeedTime_  .push_back((*theSeedHit).time());
      ootPhoSeedEnergy_.push_back((*theSeedHit).energy());
    } else {
      ootPhoSeedTime_  .push_back(-99.);
      ootPhoSeedEnergy_.push_back(-99.);
    }
    
    ootPhoSigmaIEtaIEtaFull5x5_ .push_back(iOOTPho->full5x5_sigmaIetaIeta());
    ootPhoR9Full5x5_            .push_back(iOOTPho->full5x5_r9());
    ootPhoMIPTotEnergy_         .push_back(iOOTPho->mipTotEnergy());

    nOOTPho_++;
  }

}

void ggNtuplizer::cleanupOOTPhotons() {

}
