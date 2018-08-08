#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t             nMu_;
vector<float>     muPt_;
vector<float>     muEn_;
vector<float>     muEta_;
vector<float>     muPhi_;
vector<int>       muCharge_;
vector<int>       muType_;
vector<int>       muIDbit_;
vector<float>     muD0_;
vector<float>     muDz_;
vector<float>     muSIP_;
vector<float>     muChi2NDF_;
vector<float>     muInnerD0_;
vector<float>     muInnerDz_;
vector<int>       muTrkLayers_;
vector<int>       muPixelLayers_;
vector<int>       muPixelHits_;
vector<int>       muMuonHits_;
vector<int>       muStations_;
vector<int>       muMatches_;
vector<int>       muTrkQuality_;
vector<float>     muIsoTrk_;
vector<float>     muPFChIso_;
vector<float>     muPFPhoIso_;
vector<float>     muPFNeuIso_;
vector<float>     muPFPUIso_;
vector<ULong64_t> muFiredTrgs_;
vector<ULong64_t> muFiredL1Trgs_;
vector<float>     muInnervalidFraction_;
vector<float>     musegmentCompatibility_;
vector<float>     muchi2LocalPosition_;
vector<float>     mutrkKink_;
vector<float>     muBestTrkPtError_;
vector<float>     muBestTrkPt_;
vector<int>       muBestTrkType_;

void ggNtuplizer::branchesMuons(TTree* tree) {

  tree->Branch("nMu",           &nMu_);
  tree->Branch("muPt",          &muPt_);
  tree->Branch("muEn",          &muEn_);
  tree->Branch("muEta",         &muEta_);
  tree->Branch("muPhi",         &muPhi_);
  tree->Branch("muCharge",      &muCharge_);
  tree->Branch("muType",        &muType_);
  tree->Branch("muIDbit",       &muIDbit_);
  tree->Branch("muD0",          &muD0_);
  tree->Branch("muDz",          &muDz_);
  tree->Branch("muSIP",         &muSIP_);
  tree->Branch("muChi2NDF",     &muChi2NDF_);
  tree->Branch("muInnerD0",     &muInnerD0_);
  tree->Branch("muInnerDz",     &muInnerDz_);
  tree->Branch("muTrkLayers",   &muTrkLayers_);
  tree->Branch("muPixelLayers", &muPixelLayers_);
  tree->Branch("muPixelHits",   &muPixelHits_);
  tree->Branch("muMuonHits",    &muMuonHits_);
  tree->Branch("muStations",    &muStations_);
  tree->Branch("muMatches",     &muMatches_);
  tree->Branch("muTrkQuality",  &muTrkQuality_);
  tree->Branch("muIsoTrk",      &muIsoTrk_);
  tree->Branch("muPFChIso",     &muPFChIso_);
  tree->Branch("muPFPhoIso",    &muPFPhoIso_);
  tree->Branch("muPFNeuIso",    &muPFNeuIso_);
  tree->Branch("muPFPUIso",     &muPFPUIso_);
  tree->Branch("muFiredTrgs",   &muFiredTrgs_);
  tree->Branch("muFiredL1Trgs", &muFiredL1Trgs_);
  tree->Branch("muInnervalidFraction",   &muInnervalidFraction_);
  tree->Branch("musegmentCompatibility", &musegmentCompatibility_);
  tree->Branch("muchi2LocalPosition",    &muchi2LocalPosition_);
  tree->Branch("mutrkKink",              &mutrkKink_);
  tree->Branch("muBestTrkPtError",       &muBestTrkPtError_);
  tree->Branch("muBestTrkPt",            &muBestTrkPt_);
  tree->Branch("muBestTrkType",          &muBestTrkType_);
}

void ggNtuplizer::fillMuons(const edm::Event& e, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  muPt_                  .clear();
  muEn_                  .clear();
  muEta_                 .clear();
  muPhi_                 .clear();
  muCharge_              .clear();
  muType_                .clear();
  muIDbit_               .clear();
  muD0_                  .clear();
  muDz_                  .clear();
  muSIP_                 .clear();
  muChi2NDF_             .clear();
  muInnerD0_             .clear();
  muInnerDz_             .clear();
  muTrkLayers_           .clear();
  muPixelLayers_         .clear();
  muPixelHits_           .clear();
  muMuonHits_            .clear();
  muStations_            .clear();
  muMatches_             .clear();
  muTrkQuality_          .clear();
  muIsoTrk_              .clear();
  muPFChIso_             .clear();
  muPFPhoIso_            .clear();
  muPFNeuIso_            .clear();
  muPFPUIso_             .clear();
  muFiredTrgs_           .clear();
  muFiredL1Trgs_         .clear();
  muInnervalidFraction_  .clear();
  musegmentCompatibility_.clear();
  muchi2LocalPosition_   .clear();
  mutrkKink_             .clear();
  muBestTrkPtError_      .clear();
  muBestTrkPt_           .clear();
  muBestTrkType_         .clear();

  nMu_ = 0;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  if (!muonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
    return;
  }

  for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {

    if (iMu->pt() < 3) continue;
    if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;

    muPt_    .push_back(iMu->pt());
    muEn_    .push_back(iMu->energy());
    muEta_   .push_back(iMu->eta());
    muPhi_   .push_back(iMu->phi());
    muCharge_.push_back(iMu->charge());
    muType_  .push_back(iMu->type());
    muD0_    .push_back(iMu->muonBestTrack()->dxy(pv));
    muDz_    .push_back(iMu->muonBestTrack()->dz(pv));
    muSIP_   .push_back(fabs(iMu->dB(pat::Muon::PV3D))/iMu->edB(pat::Muon::PV3D));

    int tmpmuIDbit = 0;
    if (iMu->passed(reco::Muon::CutBasedIdLoose))        tmpmuIDbit += pow(2,  0);
    if (iMu->passed(reco::Muon::CutBasedIdMedium))       tmpmuIDbit += pow(2,  1);
    if (iMu->passed(reco::Muon::CutBasedIdMediumPrompt)) tmpmuIDbit += pow(2,  2);
    if (iMu->passed(reco::Muon::CutBasedIdTight))        tmpmuIDbit += pow(2,  3);
    if (iMu->passed(reco::Muon::CutBasedIdGlobalHighPt)) tmpmuIDbit += pow(2,  4);
    if (iMu->passed(reco::Muon::CutBasedIdTrkHighPt))    tmpmuIDbit += pow(2,  5);
    if (iMu->passed(reco::Muon::PFIsoVeryLoose))         tmpmuIDbit += pow(2,  6);
    if (iMu->passed(reco::Muon::PFIsoLoose))             tmpmuIDbit += pow(2,  7);
    if (iMu->passed(reco::Muon::PFIsoMedium))            tmpmuIDbit += pow(2,  8);
    if (iMu->passed(reco::Muon::PFIsoTight))             tmpmuIDbit += pow(2,  9);
    if (iMu->passed(reco::Muon::PFIsoVeryTight))         tmpmuIDbit += pow(2, 10);
    if (iMu->passed(reco::Muon::TkIsoLoose))             tmpmuIDbit += pow(2, 11);
    if (iMu->passed(reco::Muon::TkIsoTight))             tmpmuIDbit += pow(2, 12);
    if (iMu->passed(reco::Muon::SoftCutBasedId))         tmpmuIDbit += pow(2, 13);
    if (iMu->passed(reco::Muon::SoftMvaId))              tmpmuIDbit += pow(2, 14);
    if (iMu->passed(reco::Muon::MvaLoose))               tmpmuIDbit += pow(2, 15);
    if (iMu->passed(reco::Muon::MvaMedium))              tmpmuIDbit += pow(2, 16);
    if (iMu->passed(reco::Muon::MvaTight))               tmpmuIDbit += pow(2, 17);
    if (iMu->passed(reco::Muon::MiniIsoLoose))           tmpmuIDbit += pow(2, 18);
    if (iMu->passed(reco::Muon::MiniIsoMedium))          tmpmuIDbit += pow(2, 19);
    if (iMu->passed(reco::Muon::MiniIsoTight))           tmpmuIDbit += pow(2, 20);
    if (iMu->passed(reco::Muon::MiniIsoVeryTight))       tmpmuIDbit += pow(2, 21);
    muIDbit_.push_back(tmpmuIDbit);

    muFiredTrgs_  .push_back(matchMuonTriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()));
    muFiredL1Trgs_.push_back(matchL1TriggerFilters(iMu->pt(), iMu->eta(), iMu->phi()));

    muBestTrkPtError_        .push_back(iMu->muonBestTrack()->ptError());
    muBestTrkPt_             .push_back(iMu->muonBestTrack()->pt());
    muBestTrkType_           .push_back(iMu->muonBestTrackType());
    musegmentCompatibility_  .push_back(iMu->segmentCompatibility());
    muchi2LocalPosition_     .push_back(iMu->combinedQuality().chi2LocalPosition);
    mutrkKink_               .push_back(iMu->combinedQuality().trkKink);

    const reco::TrackRef glbmu = iMu->globalTrack();
    const reco::TrackRef innmu = iMu->innerTrack();

    if (glbmu.isNull()) {
      muChi2NDF_ .push_back(-99.);
      muMuonHits_.push_back(-99);
    } else {
      muChi2NDF_.push_back(glbmu->normalizedChi2());
      muMuonHits_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
    }

    if (innmu.isNull()) {
      muInnerD0_     .push_back(-99.);
      muInnerDz_     .push_back(-99.);
      muTrkLayers_   .push_back(-99);
      muPixelLayers_ .push_back(-99);
      muPixelHits_   .push_back(-99);
      muTrkQuality_  .push_back(-99);

      muInnervalidFraction_ .push_back(-99);
    } else {
      muInnerD0_     .push_back(innmu->dxy(pv));
      muInnerDz_     .push_back(innmu->dz(pv));
      muTrkLayers_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
      muPixelLayers_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
      muPixelHits_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
      muTrkQuality_  .push_back(innmu->quality(reco::TrackBase::highPurity));

      muInnervalidFraction_ .push_back(innmu->validFraction());
    }

    muStations_   .push_back(iMu->numberOfMatchedStations());
    muMatches_    .push_back(iMu->numberOfMatches());
    muIsoTrk_     .push_back(iMu->isolationR03().sumPt);
    muPFChIso_    .push_back(iMu->pfIsolationR04().sumChargedHadronPt);
    muPFPhoIso_   .push_back(iMu->pfIsolationR04().sumPhotonEt);
    muPFNeuIso_   .push_back(iMu->pfIsolationR04().sumNeutralHadronEt);
    muPFPUIso_    .push_back(iMu->pfIsolationR04().sumPUPt);

    nMu_++;
  }

}
